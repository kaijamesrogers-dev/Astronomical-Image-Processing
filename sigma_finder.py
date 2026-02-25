from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

# =============================================================================
# User settings
# =============================================================================

SEED_SIGMA_MULTIPLIER = 3.0
RING_SIGMA_MULTIPLIER = 0.5
FIT_RANGE = (3300, 3550)

# =============================================================================
# Section 5.2: Reading the Data
# =============================================================================

with fits.open("mosaic.fits") as hdulist:
    data = hdulist[0].data

# Define crop region (y_start:y_end, x_start:x_end)
crop_region = (slice(1050, 1700), slice(390, 1300))

# Crop data
data = data[crop_region]


# =============================================================================
# Gaussian fitting
# =============================================================================


def fit_gaussian(data_array, title):
    xmin, xmax = FIT_RANGE
    flat = np.asarray(data_array).ravel()
    in_range = flat[(flat >= xmin) & (flat <= xmax)]

    if in_range.size == 0:
        raise ValueError(f"No pixels found in fit range [{xmin}, {xmax}].")

    bin_centres = np.arange(xmin, xmax + 1)
    bin_edges = np.arange(xmin - 0.5, xmax + 1.5, 1)
    counts, _ = np.histogram(in_range, bins=bin_edges)

    def gaussian(x, A, mu, sigma):
        return A * np.exp(-(x - mu) ** 2 / (2 * sigma**2))

    A0 = max(counts.max(), 1)
    mu0 = bin_centres[np.argmax(counts)]
    sigma0 = max(np.std(in_range), 1e-6)

    popt, _ = curve_fit(gaussian, bin_centres, counts, p0=[A0, mu0, sigma0])
    A, mu, sigma = popt

    model = gaussian(bin_centres, A, mu, sigma)
    residuals = counts - model
    dof = max(len(counts) - 3, 1)
    reduced_chi2 = float(np.sum((residuals**2) / np.maximum(model, 1.0)) / dof)

    plt.figure()
    plt.hist(in_range, bins=bin_edges, alpha=0.7, label="Pixels")
    plt.plot(bin_centres, model, "r-", label="Gaussian fit")
    plt.title(title)
    plt.xlabel("Pixel values")
    plt.ylabel("Number of pixels")
    plt.legend()

    print(f"{title}: mu={mu:.3f}, sigma={sigma:.3f}, reduced_chi2={reduced_chi2:.4f}")

    return float(mu), float(abs(sigma)), reduced_chi2


# =============================================================================
# Source detection
# =============================================================================


def detect_sources(image_data, mu, sigma, seed_sigma_mult, ring_sigma_mult):
    seed_threshold = mu + seed_sigma_mult * sigma
    ring_threshold = mu + ring_sigma_mult * sigma
    min_radius = 1

    mask = np.zeros(image_data.shape, dtype=bool)
    sources = []

    print(f"Seed threshold (mu + {seed_sigma_mult}*sigma) = {seed_threshold:.2f}")
    print(f"Ring threshold (mu + {ring_sigma_mult}*sigma) = {ring_threshold:.2f}")

    height, width = image_data.shape
    iteration = 0

    while True:
        masked_data = np.copy(image_data).astype(float)
        masked_data[mask] = -np.inf

        max_value = np.max(masked_data)
        if max_value < seed_threshold:
            break

        y, x = np.unravel_index(np.argmax(masked_data), masked_data.shape)

        found_radius = 0
        for r in range(1, 101):
            ymin = max(0, y - r)
            ymax = min(height, y + r + 1)
            xmin = max(0, x - r)
            xmax = min(width, x + r + 1)

            cutout = masked_data[ymin:ymax, xmin:xmax]
            yy, xx = np.ogrid[ymin:ymax, xmin:xmax]
            dist = np.sqrt((xx - x) ** 2 + (yy - y) ** 2)

            in_ring = (dist > r - 1) & (dist <= r)
            if np.count_nonzero(in_ring) == 0:
                break

            ring_median = np.median(cutout[in_ring])
            if ring_median < ring_threshold:
                found_radius = r - 1
                break
        else:
            found_radius = r

        if found_radius <= min_radius:
            mask[y, x] = True
            continue

        radius = int(1.5 * found_radius)
        if radius <= 6:
            radius = 6

        sources.append((x, y, max_value, radius))

        ymin2 = max(0, y - radius)
        ymax2 = min(height, y + radius + 1)
        xmin2 = max(0, x - radius)
        xmax2 = min(width, x + radius + 1)

        yy2, xx2 = np.ogrid[0:(ymax2 - ymin2), 0:(xmax2 - xmin2)]
        dist2 = np.sqrt((xx2 - (x - xmin2)) ** 2 + (yy2 - (y - ymin2)) ** 2)
        mask[ymin2:ymax2, xmin2:xmax2][dist2 <= radius] = True

        iteration += 1
        if iteration % 100 == 0:
            print(f"  Detected {iteration} sources...")

    print(f"Detected {len(sources)} sources above threshold")
    print(f"Masked fraction of image = {100 * np.mean(mask):.2f}%")

    background_only = image_data[~mask]
    return sources, mask, background_only


# =============================================================================
# Run
# =============================================================================

mu1, sigma1, chi2_1 = fit_gaussian(data, "Gaussian Fit 1: Before Source Removal")

sources, source_mask, background_only = detect_sources(
    data,
    mu1,
    sigma1,
    seed_sigma_mult=SEED_SIGMA_MULTIPLIER,
    ring_sigma_mult=RING_SIGMA_MULTIPLIER,
)

mu2, sigma2, chi2_2 = fit_gaussian(background_only, "Gaussian Fit 2: After Source Removal")

print("\nFit comparison:")
print(f"Before removal -> mu={mu1:.3f}, sigma={sigma1:.3f}, reduced_chi2={chi2_1:.4f}")
print(f"After removal  -> mu={mu2:.3f}, sigma={sigma2:.3f}, reduced_chi2={chi2_2:.4f}")

plt.show()
