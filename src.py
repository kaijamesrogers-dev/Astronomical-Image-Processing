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

hdulist = fits.open("mosaic.fits")
data = hdulist[0].data
header = hdulist[0].header

#print(header)
#print(data)

# =============================================================================
# Section 5.2.6: The Statistics of the Image
# =============================================================================

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

popt, pcov = curve_fit(gaussian, bin_centres, counts, p0=[A0, mu0, sigma0])
A, mu, sigma = popt

plt.hist(masked_data, bins=(xmax - xmin))
plt.plot(bin_centres, gaussian(bin_centres, A, mu, sigma))

plt.xlabel("Pixel values")
plt.ylabel("Number of pixels")
plt.show()

# =============================================================================
# Section 5.3: Source Detection
# =============================================================================

# TODO: Implement source detection
# - Find highest pixel value in image
# - Extract photometry for brightest source
# - Use mask image to flag processed areas
# - Iterate to find all sources above detection threshold



# =============================================================================
# Section 5.4: Source Photometry
# =============================================================================

# TODO: Implement aperture photometry
# - Count pixels within fixed aperture (diameter ~3" = 12 pixels)
# - Use annular reference aperture for local background
# - Subtract background contribution from aperture flux
# - Calculate flux for each detected source



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