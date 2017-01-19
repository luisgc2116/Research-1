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


hdu = fits.open("4imAF.fits")

# Data to plot is in rotated frame (2)
h2 = hdu[0].header
w2 = wcs.WCS(h2)


import numpy as np 
import pywcs 
import pyfits 
from numpy import sin, cos, radians 

def read_fits(name, hdu):
	hdulist = pyfits.open(name) 
	img_hdu = hdulist[hdu] 
	wcs = pywcs.WCS(img_hdu.header) 

	return img_hdu, hdulist, wcs 
#DEBUG = False 
def rotate(degs): 
	rads = radians(degs) 
	s = sin(rads) 
	c = cos(rads) 
	
	return np.array([[c, -s], [s, c]]) 

def write_fits(hdulist, name, clobber=True, checksum=True):
	hdulist.writeto(name, clobber=clobber, checksum=checksum)

hdulist = pyfits.open('4imAF.fits')
read_fits("4imAF.fits", 0)
rotate(34)
#write_fits(hdulist, 'result1.fits')


from skimage import data
from skimage.transform import rotate
from skimage.util import img_as_ubyte


image = data.camera()
rotate(image, 34, resize = True)

img = img_as_ubyte('4imAF.fits')

rotated_img = img_as_ubyte(transform.rotate(img, 34, resize=True))

plt.show()















