# =============================================================================
# Astronomical Image Processing
# =============================================================================

from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from astropy.visualization import ZScaleInterval, LogStretch, HistEqStretch, ImageNormalize

# =============================================================================
# Section 5.2: Reading the Data
# =============================================================================

# import data
hdulist = fits.open("mosaic.fits")
data = hdulist[0].data
header = hdulist[0].header

# Define crop region (y_start:y_end, x_start:x_end)
crop_region = (slice(1050, 1700), slice(390, 1300))

# crop data
data = data[crop_region]

# print header and data
#print(header)
#print(data)

# =============================================================================
# Section 5.2.6: The Statistics of the Image
# =============================================================================

# COMPLETED:
# - Masked data to isolate background pixels (3300-3550 range)
# - Created histogram of pixel values
# - Fitted Gaussian function to histogram
# - Excluded outlier bin from fit for better accuracy
# - Extracted background level (mu) and noise (sigma)
# - Plotted histogram with fitted Gaussian curve

def fit_gaussian():
    xmin, xmax = 3300, 3550
    masked_data = data[(data >= xmin) & (data <= xmax)]

    # Integer bin centres
    bin_centres = np.arange(xmin, xmax + 1)
    bin_edges = np.arange(xmin - 0.5, xmax + 1.5, 1)

    counts, bins = np.histogram(masked_data, bins=bin_edges)

    def gaussian(x, A, mu, sigma):
        return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

    A0 = counts.max()
    mu0 = bin_centres[np.argmax(counts)]
    sigma0 = np.std(masked_data)

    popt, pcov = curve_fit(gaussian, bin_centres, counts, p0=[A0, mu0, sigma0])
    A, mu, sigma = popt

    plt.hist(masked_data, bins=bin_edges)
    plt.plot(bin_centres, gaussian(bin_centres, A, mu, sigma))

    print(f"Background level (mu) = {mu:.1f}, Noise (sigma) = {sigma:.1f}")

    plt.xlabel("Pixel values")
    plt.ylabel("Number of pixels")

    return mu, sigma


# =============================================================================
# Section 5.3: Source Detection
# =============================================================================

# COMPLETED:
# - Seed pixel must be above 4 sigma to be a candidate
# - Expanding rings: grow radius from 1 outward, compute mean of each ring
# - Stop when ring mean falls below 3 sigma — record radius
# - Reject sources with radius < 2 (hot pixels), no cap on radius
# - Mask each detected source using its measured radius
# - Stored results in sources list: (x, y, peak_value, radius)

def detect_sources():
    mu, sigma = fit_gaussian()

    # Thresholds
    seed_threshold = mu + 4 * sigma       # initial pixel must exceed this
    ring_threshold = mu + 3 * sigma       # expanding disc average must exceed this
    min_radius = 2                         # reject single hot pixels

    # Create mask image to track processed pixels
    mask = np.zeros(data.shape, dtype=bool)

    # List to store detected sources (x, y, peak_value, radius)
    sources = []

    print(f"Seed threshold (4σ) = {seed_threshold:.1f}")
    print(f"Ring threshold (3σ) = {ring_threshold:.1f}")
    print(f"Background level (mu) = {mu:.1f}, Noise (sigma) = {sigma:.1f}")

    # Iteratively find sources
    iteration = 0
    while True:
        # Create a copy of data with masked regions set to very low value
        masked_data = np.copy(data).astype(float)
        masked_data[mask] = -np.inf

        # Find highest pixel value
        max_value = np.max(masked_data)

        # Check if above seed threshold
        if max_value < seed_threshold:
            break

        # Find position of maximum
        max_index = np.argmax(masked_data)
        y, x = np.unravel_index(max_index, masked_data.shape)

        # Expand outward from the peak pixel in rings
        # At each radius r, compute the mean of pixels in the ring (r-1 < dist <= r)
        # Stop when the ring average falls below the threshold
        found_radius = 0
        max_expand = max(100)
        for r in range(1, max_expand + 1):
            # Define local cutout bounds around the source
            ymin = max(0, y - r)
            ymax = y + r + 1
            xmin = max(0, x - r)
            xmax = x + r + 1

            cutout = data[ymin:ymax, xmin:xmax]
            yy, xx = np.ogrid[ymin:ymax, xmin:xmax] # list of xy values within rectangle around peak point
            dist = np.sqrt((xx - x)**2 + (yy - y)**2)

            in_ring = (dist > r - 1) & (dist <= r)
            if np.count_nonzero(in_ring) == 0:
                break
            ring_mean = np.mean(cutout[in_ring])

            if ring_mean < ring_threshold:
                # Ring average dropped below threshold — stop expanding
                found_radius = r - 1
                break
        else:
            # Reached image boundary without dropping below threshold
            found_radius = r

        if found_radius < min_radius:
            # Reject: too small (likely a hot pixel)
            mask[y, x] = True
            continue

        # Store source with its measured radius
        sources.append((x, y, max_value, found_radius))

        # Mask the detected object using its measured radius
        ymin = max(0, y - found_radius)
        ymax = min(height, y + found_radius + 1)
        xmin = max(0, x - found_radius)
        xmax = min(width, x + found_radius + 1)
        yy, xx = np.ogrid[ymin:ymax, xmin:xmax]
        dist = np.sqrt((xx - x)**2 + (yy - y)**2)
        mask[ymin:ymax, xmin:xmax][dist <= found_radius] = True

        iteration += 1
        if iteration % 100 == 0:
            print(f"  Detected {iteration} sources...")

    print(f"\nDetected {len(sources)} sources above threshold")
    if sources:
        print(f"Brightest source at (x={sources[0][0]}, y={sources[0][1]}) with value {sources[0][2]:.1f}, radius={sources[0][3]}px")

    return sources, mu, sigma

# =============================================================================
# Section 5.4: Source Photometry
# =============================================================================

# COMPLETED:
# - Used each source's measured radius as the aperture
# - Background annulus scales with source radius (gap=2px, width=7px)
# - Subtracted background contribution from aperture flux
# - Calculated flux for each detected source

def aperture_photometry():
    sources, mu, sigma = detect_sources()

    # Annulus offset relative to each source's measured radius
    annulus_gap = 2       # gap between aperture edge and annulus start
    annulus_width = 7     # width of the background annulus

    height, width = data.shape

    # Build a global source mask to exclude sources from background annuli
    source_mask = np.zeros(data.shape, dtype=bool)
    for (sx, sy, _, sr) in sources:
        ymin = max(0, sy - sr)
        ymax = min(height, sy + sr + 1)
        xmin = max(0, sx - sr)
        xmax = min(width, sx + sr + 1)
        yy, xx = np.ogrid[ymin:ymax, xmin:xmax]
        dist = np.sqrt((xx - sx)**2 + (yy - sy)**2)
        source_mask[ymin:ymax, xmin:xmax][dist <= sr] = True

    results = []

    for i, (x, y, peak, radius) in enumerate(sources):
        annulus_inner = radius + annulus_gap
        annulus_outer = annulus_inner + annulus_width

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

        # Source aperture: sum all pixels within this source's measured radius
        in_aperture = dist <= radius
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

        results.append((x, y, peak, net_flux, bg_per_pixel, radius))

        if (i + 1) % 100 == 0:
            print(f"  Photometry for {i + 1}/{len(sources)} sources...")

    print(f"\nCompleted aperture photometry for {len(results)} sources")
    print(results[0])  # print first source's photometry results for verification

    return results

# =============================================================================
# Section 5.5: Calibrating the Fluxes
# =============================================================================

# COMPLETED:
# - Read MAGZPT and MAGZRR from FITS header
# - Convert net flux (counts) to calibrated magnitudes: m = MAGZPT - 2.5 * log10(counts)
# - Skip sources with non-positive net flux (cannot take log)
# - Calculate magnitude error from MAGZRR

def calibrate_fluxes():
    results = aperture_photometry()

    # Read zero point and its error from FITS header
    magzpt = header['MAGZPT']
    magzrr = header['MAGZRR']

    print(f"Zero point: MAGZPT = {magzpt}, MAGZRR = {magzrr}")

    calibrated = []

    for (x, y, _, net_flux, _, _) in results:
        if net_flux < 1:
            continue  # skip sources with negligible flux (< 1 count)

        # m = ZP_inst + mag_i = ZP_inst - 2.5 * log10(counts)
        mag = magzpt - 2.5 * np.log10(net_flux)
        mag_err = magzrr

        calibrated.append((x, y, net_flux, mag, mag_err))

    print(f"Calibrated {len(calibrated)} sources (skipped {len(results) - len(calibrated)} with non-positive flux)")

    return calibrated

# =============================================================================
# Section 5.6: Producing the Catalogue
# =============================================================================

# COMPLETED:
# - Combine photometry results with calibrated magnitudes
# - Save ASCII catalogue with: x, y, peak, net_flux, bg_per_pixel, magnitude, mag_error
# - Use numpy.savetxt for clean formatted output

def produce_catalogue(results, calibrated):
    # Look up peak value and background from photometry results
    result_lookup = {(x, y): (peak, bg_per_pixel) for (x, y, peak, _, bg_per_pixel, _) in results}

    rows = []
    for (x, y, net_flux, mag, mag_err) in calibrated:
        peak, bg_per_pixel = result_lookup[(x, y)]
        rows.append([x, y, peak, net_flux, bg_per_pixel, mag, mag_err])

    catalogue = np.array(rows)

    np.savetxt('catalogue.csv', catalogue,
               header='x  y  peak_value  net_flux  bg_per_pixel  magnitude  mag_error',
               fmt=['%d', '%d', '%.2f', '%.2f', '%.2f', '%.4f', '%.4f'])

    print(f"Catalogue saved: {len(catalogue)} sources to catalogue.csv")

# =============================================================================
# Section 5.7: Analyzing the Data
# =============================================================================

# COMPLETED:
# - Calculate N(< m) vs m (cumulative number counts per deg²)
# - Plot log(N) vs magnitude with Poisson error bars
# - Overlay theoretical relation: log N(m) = 0.6m + constant

def number_counts(calibrated):
    # Extract magnitudes
    mags = np.array([mag for (_, _, _, mag, _) in calibrated])

    # Image area in square degrees
    pixel_scale = 0.258  # arcsec per pixel
    height, width = data.shape
    area_deg2 = height * width * pixel_scale**2 / 3600**2

    print(f"Image area: {height} x {width} pixels = {area_deg2:.6f} deg²")

    # Magnitude bin edges (0.5 mag steps)
    mag_bins = np.arange(np.floor(mags.min()), np.ceil(mags.max()) + 0.5, 0.5)

    # Cumulative counts: N(< m) for each bin edge
    N_cumulative = np.array([np.sum(mags <= m) for m in mag_bins])

    # Only keep bins with at least 1 source
    valid = N_cumulative > 0
    mag_plot = mag_bins[valid]
    N_raw = N_cumulative[valid]

    # Normalize to per square degree
    N_per_deg2 = N_raw / area_deg2

    # Poisson error bars propagated to log10: sigma_log10 = 1 / (sqrt(N) * ln(10))
    log10_N = np.log10(N_per_deg2)
    log10_err = 1.0 / (np.sqrt(N_raw) * np.log(10))

    # Plot number counts
    plt.figure()
    plt.errorbar(mag_plot, log10_N, yerr=log10_err, fmt='ko', markersize=4, label='Observed')

    # Fit a straight line to the data: log10(N) = gradient * m + intercept
    gradient, intercept = np.polyfit(mag_plot, log10_N, 1)
    print(f"Measured gradient: {gradient:.4f} (theoretical: 0.6)")

    # Theoretical relation: log10(N) = 0.6m + constant
    mid = len(mag_plot) // 2
    constant = log10_N[mid] - 0.6 * mag_plot[mid]
    mag_theory = np.linspace(mag_plot[0], mag_plot[-1], 100)
    plt.plot(mag_theory, 0.6 * mag_theory + constant, 'r--', label='Theory: 0.6m + const')
    plt.plot(mag_theory, gradient * mag_theory + intercept, 'b--', label=f'Fit: {gradient:.4f}m + const')

    plt.xlabel('Magnitude')
    plt.ylabel('log$_{10}$(N(< m)) [per deg$^2$]')
    plt.legend()
    plt.title('Cumulative Number Counts')

    print(f"Magnitude range: {mags.min():.2f} to {mags.max():.2f}")
    print(f"Total sources: {len(mags)}")

# =============================================================================
# Source Detection Visualization
# =============================================================================

# COMPLETED:
# - Create a 2D image marking detected source positions
# - Overlay detection map on original image for visual comparison
# - Show which sources were detected and their spatial distribution

def visualise_sources(sources, calibrated, label_map):
    """
    Create detection map and overlay with original image

    Parameters:
    sources: list of (x, y, peak_value, radius, n_pixels) tuples
    calibrated: list of calibrated source data
    label_map: 2D array where each pixel is labelled with its source index (-1 = background)
    """
    # Create figure with subplots
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Plot 1: Original image with histogram equalization
    norm = ImageNormalize(data, stretch=HistEqStretch(data))
    axes[0].imshow(data, cmap='gray', origin='lower', norm=norm)
    axes[0].set_title('Original Image')
    axes[0].set_xlabel('X (pixels)')
    axes[0].set_ylabel('Y (pixels)')

    # Plot 2: Flood-fill detection map — show actual detected regions (binary)
    detection_map = np.zeros_like(data, dtype=float)
    detection_map[label_map >= 0] = 1.0
    axes[1].imshow(detection_map, cmap='gray', origin='lower', vmin=0, vmax=1)
    axes[1].set_title(f'Detected Sources (n={len(sources)})')
    axes[1].set_xlabel('X (pixels)')
    axes[1].set_ylabel('Y (pixels)')

    # Plot 3: Original image with point markers at source centres
    axes[2].imshow(data, cmap='gray', origin='lower', norm=norm)
    if len(sources) > 0:
        src_x = np.array([s[0] for s in sources])
        src_y = np.array([s[1] for s in sources])
        axes[2].plot(src_x, src_y, 'r+', markersize=5, markeredgewidth=0.5, label=f'Detected ({len(sources)})')
        axes[2].legend()
    axes[2].set_title('Overlay: Original + Detections')
    axes[2].set_xlabel('X (pixels)')
    axes[2].set_ylabel('Y (pixels)')

    plt.tight_layout()

    print(f"\nVisualisation created: {len(sources)} sources detected")
    print(f"Calibrated sources: {len(calibrated)}")

# =============================================================================
# Magnitude vs Source Size
# =============================================================================

def magnitude_vs_size(sources, calibrated):
    # Build lookup from (x, y) -> radius from sources
    radius_lookup = {(x, y): r for (x, y, _, r, _) in sources}

    mags = []
    radii = []
    for (x, y, _, mag, _) in calibrated:
        if (x, y) in radius_lookup:
            mags.append(mag)
            radii.append(radius_lookup[(x, y)])

    mags = np.array(mags)
    radii = np.array(radii)

    plt.figure()
    plt.scatter(radii, mags, s=5, alpha=0.5)
    plt.xlabel('Source Radius (pixels)')
    plt.ylabel('Magnitude')
    plt.title('Magnitude vs Source Size')
    plt.gca().invert_yaxis()  # brighter objects (lower mag) at the top

    print(f"Plotted {len(mags)} sources: radius range {radii.min()}-{radii.max()} px")

# =============================================================================
# Running Code
# =============================================================================

#fit_gaussian()
sources, mu, sigma, label_map = detect_sources()
results = aperture_photometry()
calibrated = calibrate_fluxes()
produce_catalogue(results, calibrated)
number_counts(calibrated)
visualise_sources(sources, calibrated, label_map)
magnitude_vs_size(sources, calibrated)

plt.show()