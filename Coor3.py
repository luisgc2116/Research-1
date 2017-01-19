import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar
from astropy.io import fits

from math import exp
import math 
from matplotlib import pyplot
from scipy.interpolate import UnivariateSpline

from astropy import wcs
from astropy import coordinates
from astropy import units as u

import pywcsgrid2



import aplpy
#from astropy import Skycoord
from astropy import coordinates
from astropy import units as u

hdu_list = fits.open('4imAF.fits')
hdu_list.info()

image_data = hdu_list[0].image_data

header = fits.getheader('4imAF.fits')


plt.imshow(image_data, cmap = 'gray')
plt.colorbar()

plt.show()

#x_pix, y_pix = fig.world2pixel(45.3332, 22.12)
