from astropy.io import fits
import matplotlib.pyplot as plt

# import data
hdulist = fits.open("mosaic.fits")
data = hdulist[0].data
header = hdulist[0].header

# print header and data
#print(header)
#print(data)

# mask data
xmin, xmax = 3000, 4000
masked_data = data[(data >= xmin) & (data <= xmax)]

# plot masked data
plt.hist(masked_data, bins=100)
plt.xlabel("Pixel value")
plt.ylabel("Number of pixels")
plt.show()
