#!/usr/bin/env python
import numpy as np
import pylab as P

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar
from astropy.io import fits
from math import exp
from scipy.interpolate import UnivariateSpline
import math 
from astropy.modeling import models, fitting
from matplotlib import pyplot


#Opens file 
Wise1 = fits.open('4imDHN.fits')# im1\Cut2image4.fits
im1 = Wise1[0].data

#unravels data into a one-dimensional array to be passed onto the function P.hist for a histogram
b = im1.ravel()

b[np.isnan(b)] = 0
b[np.isinf(b)] = 0

#mu, sigma = 200, 25
#x = mu + sigma*P.randn(10000)

# the histogram of the data with histtype='step'
n, bins, patches = P.hist(b, bins = range(300,500), normed=1)# histtype='stepfilled'
P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

# add a line showing the expected distribution
#y = P.normpdf( bins, mu, sigma)
#l = P.plot(bins, y, 'k--', linewidth=1.5)

def f(x, p1, p2, p3):
	 return p1*np.exp((-np.absolute(((x-p2)/p3))**2)*np.log(2))


p0 = (80, 400, 200)
popt, pcov = curve_fit(f, bins, n, p0=p0)


P.show()
