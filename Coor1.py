import aplpy

import numpy as np 
import matplotlib
import matplotlib.pyplot as plt

from astropy import fits
#from astropy import Skycoord
from astropy import wcs
from astropy import coordinates
from astropy import units as u

hdu_list = fits.open('4imAF.fits')
hdu_list.info()

image_data = hdu_list[0].image_data

header = fits.getheader('4imAF.fits')


plt.imshow(image_data, cmap = 'gray')
plt.colorbar()

plt.show()