# =============================================================================
# Astronomical Image Processing
# =============================================================================

from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

# =============================================================================
# Section 5.2: Reading the Data
# =============================================================================

# import data
hdulist = fits.open("mosaic.fits")
data = hdulist[0].data
header = hdulist[0].header

# crop data
#data = data[1000:4000, 1000:4000]

# print header and data
#print(header)
#print(data)

# =============================================================================
# Section 5.2.6: The Statistics of the Image
# =============================================================================

# ✓ COMPLETED:
# - Masked data to isolate background pixels (3300-3550 range)
# - Created histogram of pixel values
# - Fitted Gaussian function to histogram
# - Excluded outlier bin from fit for better accuracy
# - Extracted background level (mu) and noise (sigma)
# - Plotted histogram with fitted Gaussian curve

def fit_gaussian():
    # mask data
    xmin, xmax = 3300, 3550
    masked_data = data[(data >= xmin) & (data <= xmax)]

    # histogram
    counts, bins = np.histogram(masked_data, bins=(xmax - xmin))
    bin_centres = 0.5 * (bins[:-1] + bins[1:])

    # Gaussian model
    def gaussian(x, A, mu, sigma):
        return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

    # initial guesses
    A0 = counts.max()
    mu0 = bin_centres[np.argmax(counts)]
    sigma0 = np.std(masked_data)

    # find the "wrong" bin (the highest value of the histogram)
    bad = np.argmax(counts)

    # exclude it from the fit
    ok = np.ones_like(counts, dtype=bool)
    ok[bad] = False

    popt, pcov = curve_fit(gaussian,bin_centres[ok], counts[ok], p0=[A0, mu0, sigma0])
    A, mu, sigma = popt

    plt.hist(masked_data, bins=(xmax - xmin))
    plt.plot(bin_centres, gaussian(bin_centres, A, mu, sigma))

    print(f"Background level (mu) = {mu:.1f}, Noise (sigma) = {sigma:.1f}")

    plt.xlabel("Pixel values")
    plt.ylabel("Number of pixels")

    return mu, sigma

# =============================================================================
# Section 5.3: Source Detection
# =============================================================================

# ✓ COMPLETED:
# - Set detection threshold at 5 sigma above background
# - Created boolean mask image to flag processed areas
# - Iteratively found highest pixel value in unmasked regions
# - Stored source positions (x, y) and peak values
# - Masked circular aperture (12 pixel diameter) around each detected source
# - Detected all sources above threshold
# - Stored results in sources list: (x, y, peak_value)

def detect_sources():
    mu, sigma = fit_gaussian()

    # Detection threshold: 5 sigma above background
    detection_threshold = mu + 5 * sigma

    # Create mask image to track processed pixels
    mask = np.zeros(data.shape, dtype=bool)

    # List to store detected sources (x, y, peak_value)
    sources = []

    # Aperture radius for masking (12 pixel diameter)
    aperture_radius = 6

    print(f"Starting source detection with threshold = {detection_threshold:.1f}")
    print(f"Background level (mu) = {mu:.1f}, Noise (sigma) = {sigma:.1f}")

    # Iteratively find sources
    iteration = 0
    while True:
        # Create a copy of data with masked regions set to very low value
        masked_data = np.copy(data).astype(float)
        masked_data[mask] = -np.inf

        # Find highest pixel value
        max_value = np.max(masked_data)

        # Check if above detection threshold
        if max_value < detection_threshold:
            break

        # Find position of maximum
        max_index = np.argmax(masked_data)
        y, x = np.unravel_index(max_index, masked_data.shape)

        # Store source
        sources.append((x, y, max_value))

        # Create circular mask around source, dont understand how it works
        yy, xx = np.ogrid[:data.shape[0], :data.shape[1]]
        distance = np.sqrt((xx - x)**2 + (yy - y)**2)
        mask[distance <= aperture_radius] = True

        iteration += 1
        if iteration % 100 == 0:
            print(f"  Detected {iteration} sources...")

    print(f"\nDetected {len(sources)} sources above threshold")
    print(f"Brightest source at (x={sources[0][0]}, y={sources[0][1]}) with value {sources[0][2]:.1f}")

    return sources

# =============================================================================
# Section 5.4: Source Photometry
# =============================================================================

# TODO: Implement aperture photometry
# - Count pixels within fixed aperture (diameter ~3" = 12 pixels)
# - Use annular reference aperture for local background
# - Subtract background contribution from aperture flux
# - Calculate flux for each detected source

def aperture_photometry(sources, mu):
    aperture_radius = 6    # 12 pixel diameter ≈ 3 arcsec at 0.258"/pixel
    annulus_inner = 8
    annulus_outer = 15

    height, width = data.shape

    # Build a global source mask to exclude sources from background annuli
    source_mask = np.zeros(data.shape, dtype=bool)
    for (sx, sy, _) in sources:
        ymin = max(0, sy - aperture_radius)
        ymax = min(height, sy + aperture_radius + 1)
        xmin = max(0, sx - aperture_radius)
        xmax = min(width, sx + aperture_radius + 1)
        yy, xx = np.ogrid[ymin:ymax, xmin:xmax]
        dist = np.sqrt((xx - sx)**2 + (yy - sy)**2)
        source_mask[ymin:ymax, xmin:xmax][dist <= aperture_radius] = True

    results = []

    for i, (x, y, peak) in enumerate(sources):
        # Work on a local cutout around the source for efficiency
        ymin = max(0, y - annulus_outer)
        ymax = min(height, y + annulus_outer + 1)
        xmin = max(0, x - annulus_outer)
        xmax = min(width, x + annulus_outer + 1)

        cutout = data[ymin:ymax, xmin:xmax]
        local_source_mask = source_mask[ymin:ymax, xmin:xmax]

        # Distance from source centre for each pixel in the cutout
        yy, xx = np.ogrid[ymin:ymax, xmin:xmax]
        dist = np.sqrt((xx - x)**2 + (yy - y)**2)

        # Source aperture: sum all pixels within aperture_radius
        in_aperture = dist <= aperture_radius
        aperture_sum = np.sum(cutout[in_aperture])
        n_aperture = np.count_nonzero(in_aperture)

        # Background annulus: exclude pixels belonging to any detected source
        in_annulus = (dist >= annulus_inner) & (dist <= annulus_outer)
        clean_annulus = in_annulus & ~local_source_mask

        if np.count_nonzero(clean_annulus) > 0:
            bg_per_pixel = np.mean(cutout[clean_annulus])
        else:
            bg_per_pixel = mu  # fallback to global background

        # Net flux = aperture counts - background contribution
        bg_total = bg_per_pixel * n_aperture
        net_flux = aperture_sum - bg_total

        results.append((x, y, peak, net_flux, bg_per_pixel))

        if (i + 1) % 100 == 0:
            print(f"  Photometry for {i + 1}/{len(sources)} sources...")

    print(f"\nCompleted aperture photometry for {len(results)} sources")
    print(f"Aperture radius: {aperture_radius} px, Background annulus: {annulus_inner}-{annulus_outer} px")
    

# =============================================================================
# Section 5.5: Calibrating the Fluxes
# =============================================================================

# TODO: Convert counts to magnitudes
# - Read MAGZPT and MAGZRR from FITS header
# - Convert counts to instrumental magnitudes: mag_i = -2.5 * log10(counts)
# - Apply zero point: m = MAGZPT + mag_i
# - Calculate errors including MAGZRR



# =============================================================================
# Section 5.6: Producing the Catalogue
# =============================================================================

# TODO: Create and save source catalogue
# - Store (x, y) positions of sources
# - Store total counts and background
# - Store calibrated magnitudes and errors
# - Save to ASCII file for further analysis



# =============================================================================
# Section 5.7: Analyzing the Data
# =============================================================================

# TODO: Create number count plot
# - Calculate N(< m) vs m (cumulative number counts)
# - Plot log(N(m)) vs magnitude with error bars
# - Compare to theoretical relation: log N(m) = 0.6m + constant
# - Compare with published results (e.g., Yasuda et al. 2001)

# =============================================================================
# Running Code
# =============================================================================

#fit_gaussian()

sources, mu, sigma = detect_sources()
results = aperture_photometry(sources, mu)

plt.show()