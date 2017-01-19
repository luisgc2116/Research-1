import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar
from astropy.io import fits
import pywcsgrid2

from math import exp
import math 
from matplotlib import pyplot
from scipy.interpolate import UnivariateSpline

from astropy import wcs
from astropy import coordinates
from astropy import units as u


hdu = fits.open("4imAF.fits")

# Data to plot is in rotated frame (2)
h2 = hdu[0].header
w2 = wcs.WCS(h2)

# To make the rotation easier, first shift the WCS reference pixel to
# the center of the grid
nx, ny = h2["naxis1"], h2["naxis2"]
i0, j0 = (float(nx) + 1.0)/2, (float(ny) + 1.0)/2
[ra0, dec0], = w2.wcs_pix2world([[i0, j0]], 1)
h2.update(crpix1=i0, crpix2=j0, crval1=ra0, crval2=dec0)





# First plot the image in the original (non-rotated) frame

ax2 = pywcsgrid2.axes
arcsec = 1./3600
#ax2.imshow(hdu.data, origin='lower', vmin=6, vmax=16, cmap=plt.cm.gray_r)
#ax2["fk5"].plot([ra0, ra0],[dec0 + 10*arcsec, dec0 + 15*arcsec], "r")
#ax2.grid()
#ax2.figure.savefig("pywcsgrid2-rotate-test-pixframe.pdf", dpi=150)






# Now construct a header that is in the equatorial frame (1)

h1 = h2.copy()
h1.update(
    cd1_1=-np.hypot(h2["CD1_1"], h2["CD1_2"]), 
    cd1_2=0.0, 
    cd2_1=0.0, 
    cd2_2=np.hypot(h2["CD2_1"], h2["CD2_2"]), 
    orientat=0.0)
# Finally plot the image in the new frame
ax = pywcsgrid2.axes(wcs=h1, aspect=1)
ax.set_xlim(-nx/4, 5*nx/4)
ax.set_ylim(-ny/4, 5*ny/4)
ax["fk5"].plot([ra0, ra0], 
               [dec0 + 10*arcsec, dec0 + 15*arcsec], "r")
ax[h2].imshow_affine(hdu.data, origin='lower', 
                     vmin=6, vmax=16, cmap=plt.cm.gray_r)
ax.grid()
ax.figure.savefig("pywcsgrid2-rotate-test-eqframe.pdf", dpi=150)

#hdu = fits.PrimaryHDU(im1)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('1_3imAF_1.fits')


