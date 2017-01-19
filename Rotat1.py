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


hdu = fits.open("im1.fits")

h2 = hdu[0].header
w2 = wcs.WCS(h2)


import aplpy
fig = aplpy.FITSFigure(hdu[0])

#fig.show_markers(gal, layer='marker_set_1', edgecolor='red',facecolor='none', marker='o', s=10, alpha=0.5)

#fig.set_xaxis_coord_type('latitude')

fig.show_colorscale(cmap='gist_heat')
fig.frame.set_linewidth(1)
fig.frame.set_color('black')
#fig.set_theme('publication')

fig.add_colorbar()
fig.colorbar.show()
fig.colorbar.set_location('right')
#fig.colorbar.set_axis_label_rotation(120)

fig.add_grid()
fig.grid.show()
fig.set_system_latex(True)



plt.show()


#x_pix, y_pix = fig.world2pixel(45.3332, 22.12)
