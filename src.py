from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from math import erfc, sqrt

with fits.open("mosaic.fits") as hdulist:
    image_data = hdulist[0].data.astype(float)
    header = hdulist[0].header

# Header data
print(header)

# FITS data summary
print(f"Image shape: {image_data.shape}")
print(f"Image dtype: {image_data.dtype}")

# Section 5.2.6 - The Statistics of the Image
# Create a histogram of pixel values to analyze the image statistics

pixels = image_data[np.isfinite(image_data)].ravel()
num_pixels = pixels.size

# Suppress very bright pixels (stars/blooming) to reveal the background Gaussian
high_clip = np.percentile(pixels, 99.5)
pixels_suppressed = pixels[pixels <= high_clip]

# Background/noise estimate from suppressed distribution
background_mean = np.mean(pixels_suppressed)
background_sigma = np.std(pixels_suppressed)

# Create figure with histograms
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Full histogram (log y-axis to show long positive tail)
ax1.hist(pixels, bins=300, color='blue', edgecolor='black')
ax1.set_xlabel('Pixel Value (counts)')
ax1.set_ylabel('Frequency')
ax1.set_title('Histogram of All Pixel Values')
ax1.set_yscale('log')

# Histogram with suppressed high values (zoomed in on the main Gaussian distribution)
# Clip extreme values to better visualize the background noise distribution
ax2.hist(pixels_suppressed, bins=250, color='green', edgecolor='black')
ax2.set_xlabel('Pixel Value (counts)')
ax2.set_ylabel('Frequency')
ax2.set_title(f'Histogram ({high_clip:.0f} counts clip, 99.5th percentile)')

# Mark mean and +3 sigma level used for source detection discussions
ax2.axvline(background_mean, color='black', linestyle='--', linewidth=1.5, label='Background mean')
ax2.axvline(
    background_mean + 3 * background_sigma,
    color='red',
    linestyle='--',
    linewidth=1.5,
    label='+3 sigma threshold',
)
ax2.legend()

plt.tight_layout()
plt.savefig('histogram.png', dpi=100)
if "agg" in plt.get_backend().lower():
    plt.close(fig)
else:
    plt.show()

# Print image statistics
print("\n--- Image Statistics (Section 5.2.6) ---")
print(f"Finite pixels analysed: {num_pixels:,}")
print(f"Background mean (clipped): {background_mean:.2f} counts")
print(f"Background noise sigma (clipped): {background_sigma:.2f} counts")
print(f"Median pixel value: {np.median(pixels):.2f} counts")
print(f"Min / Max pixel values: {np.min(pixels):.2f} / {np.max(pixels):.2f} counts")
print(f"Bright-pixel suppression level (99.5th percentile): {high_clip:.2f} counts")

print("\nSuggested one-sided Gaussian detection thresholds:")
for nsigma in (3, 4, 5):
    # Probability that pure noise exceeds +N sigma in one pixel
    p_tail = 0.5 * erfc(nsigma / sqrt(2.0))
    expected_false = num_pixels * p_tail
    threshold_counts = background_mean + nsigma * background_sigma
    print(
        f"{nsigma} sigma -> threshold {threshold_counts:.2f} counts, "
        f"expected noise pixels above threshold: {expected_false:.2f}"
    )
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from math import erfc, sqrt

with fits.open("mosaic.fits") as hdulist:
    image_data = hdulist[0].data.astype(float)
    header = hdulist[0].header

# Header data
print(header)

# FITS data summary
print(f"Image shape: {image_data.shape}")
print(f"Image dtype: {image_data.dtype}")

# Section 5.2.6 - The Statistics of the Image
# Create a histogram of pixel values to analyze the image statistics

pixels = image_data[np.isfinite(image_data)].ravel()
num_pixels = pixels.size

# Suppress very bright pixels (stars/blooming) to reveal the background Gaussian
high_clip = np.percentile(pixels, 99.5)
pixels_suppressed = pixels[pixels <= high_clip]

# Background/noise estimate from suppressed distribution
background_mean = np.mean(pixels_suppressed)
background_sigma = np.std(pixels_suppressed)

# Create figure with histograms
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Full histogram (log y-axis to show long positive tail)
ax1.hist(pixels, bins=300, color='blue', edgecolor='black')
ax1.set_xlabel('Pixel Value (counts)')
ax1.set_ylabel('Frequency')
ax1.set_title('Histogram of All Pixel Values')
ax1.set_yscale('log')

# Histogram with suppressed high values (zoomed in on the main Gaussian distribution)
# Clip extreme values to better visualize the background noise distribution
ax2.hist(pixels_suppressed, bins=250, color='green', edgecolor='black')
ax2.set_xlabel('Pixel Value (counts)')
ax2.set_ylabel('Frequency')
ax2.set_title(f'Histogram ({high_clip:.0f} counts clip, 99.5th percentile)')

# Mark mean and +3 sigma level used for source detection discussions
ax2.axvline(background_mean, color='black', linestyle='--', linewidth=1.5, label='Background mean')
ax2.axvline(
    background_mean + 3 * background_sigma,
    color='red',
    linestyle='--',
    linewidth=1.5,
    label='+3 sigma threshold',
)
ax2.legend()

plt.tight_layout()
plt.savefig('histogram.png', dpi=100)
if "agg" in plt.get_backend().lower():
    plt.close(fig)
else:
    plt.show()

# Print image statistics
print("\n--- Image Statistics (Section 5.2.6) ---")
print(f"Finite pixels analysed: {num_pixels:,}")
print(f"Background mean (clipped): {background_mean:.2f} counts")
print(f"Background noise sigma (clipped): {background_sigma:.2f} counts")
print(f"Median pixel value: {np.median(pixels):.2f} counts")
print(f"Min / Max pixel values: {np.min(pixels):.2f} / {np.max(pixels):.2f} counts")
print(f"Bright-pixel suppression level (99.5th percentile): {high_clip:.2f} counts")

print("\nSuggested one-sided Gaussian detection thresholds:")
for nsigma in (3, 4, 5):
    # Probability that pure noise exceeds +N sigma in one pixel
    p_tail = 0.5 * erfc(nsigma / sqrt(2.0))
    expected_false = num_pixels * p_tail
    threshold_counts = background_mean + nsigma * background_sigma
    print(
        f"{nsigma} sigma -> threshold {threshold_counts:.2f} counts, "
        f"expected noise pixels above threshold: {expected_false:.2f}"
    )
