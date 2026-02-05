from astropy.io import fits

hdulist = fits.open("mosaic.fits")

#header data
print(hdulist[0].header)

#FITS data
print(hdulist[0].data)
