import numpy as np
import astropy.wcs as wcs
from astropy.io import fits
import pywcsgrid2
import matplotlib.pyplot as plt



#opening the fits file
hdu = fits.open("4imAF.fits")
im1 = hdu[0].data

h2 = hdu[0].header
w2 = wcs.WCS(h2)

nx, ny = h2["naxis1"], h2["naxis2"]
i0, j0 = (float(nx) + 1.0)/2, (float(ny) + 1.0)/2
[ra0, dec0], = w2.wcs_pix2world([[i0, j0]], 1)
h2.update(crpix1=i0, crpix2=j0, crval1=ra0, crval2=dec0)


# First plot the image in the original (rotated) frame
ax2 = pywcsgrid2.axes(header=h2)
arcsec = 1./3600
ax2.imshow(hdu.data, origin='lower', 
           vmin=6, vmax=16, cmap=plt.cm.gray_r)
ax2["fk5"].plot([ra0, ra0], 
                [dec0 + 10*arcsec, dec0 + 15*arcsec], "r")
ax2.grid()
ax2.figure.savefig("pywcsgrid2-rotate-test-pixframe.pdf", dpi=150)





#Rotation
'''from astropy.modeling.models import Rotation2D

SkyRotation = Rotation2D.rename('SkyRotation')

'''
'''
I don't know how to put the inputs; in this case, we have the 1. rotation angle and 2. the input fits file.
This is a very trivial question, I know, but for whatever reason it's proving to be a lot more difficult than
what it's supposed to be.

'''


#This is the kind of attempt I did. Here the angle is 36 degrees and im1 is the array.
#rotation = SkyRotation('36')
#rotation.im1



