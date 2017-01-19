from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
import numpy as np
from matplotlib import pyplot as plt


c1 = SkyCoord(10, 20, unit='deg')  # Defaults to ICRS



c2 = SkyCoord([1, 2, 3], [-30, 45, 8], frame='icrs', unit='deg')

sc = SkyCoord("1h12m43.2s", "+1d12m43s", frame=Galactic, unit='deg')



import pyfits

from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table

file = get_pkg_data_filename('im1.fits')

hdu = fits.open(file)[0]
wcs = WCS(hdu.header)

#sc = SkyCoord(ra, dec, frame='icrs')

fig = plt.figure()
fig.add_subplot(111, projection=wcs)
plt.imshow(hdu.data, origin='lower', cmap=plt.cm.viridis)
plt.xlabel('RA')
plt.ylabel('Dec')

#plt.show()

data = fits.getdata('im1.fits')
psc = Table(data)

hdulist_im = fits.open('im1.fits')

# Extract image and header
image = hdulist_im[0].data
header = hdulist_im[0].header

# Instantiate WCS object
w = WCS(header)

