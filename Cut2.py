import numpy as np
from astropy.modeling.models import Gaussian2D
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy import wcs

import pyfits


f = pyfits.open("im1.fits")
h, data = f[0].header, f[0].data
#w = wcs.WCS(h)

import matplotlib.pyplot as plt

plt.grid(color='white', ls='solid')
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.imshow(data, origin='lower')

plt.show()

from astropy.nddata import Cutout2D
from astropy import units as u
from astropy import coordinates


#Position - center of the image, size (nx, ny)
position = (697, 608)
size = (150, 1300)     # of pixels


cutout = Cutout2D(data, position, size)
plt.imshow(cutout.data, origin='lower')

plt.show()


plt.imshow(data, origin='lower')



cutout.plot_on_original(color='white', alpha = 1, linewidth = 2)


plt.show()
