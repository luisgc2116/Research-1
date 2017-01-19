import numpy as np
from astropy.modeling.models import Gaussian2D


y, x = np.mgrid[0:500, 0:500]
data = Gaussian2D(1, 50, 100, 10, 5, theta=0.5)(x, y)

import matplotlib.pyplot as plt
plt.imshow(data, origin='lower')

plt.show()

from astropy.nddata import Cutout2D
from astropy import units as u


#Position - center of the image, size (nx, ny)
position = (49.7, 100.1)

#size = (40, 50)     # of pixels
#size = (40*u.pixel, 50*u.pixel)
#size = 40
size = 40 * u.pixel

cutout = Cutout2D(data, position, size)

plt.imshow(cutout.data, origin='lower')

plt.show()

plt.imshow(data, origin='lower')
cutout.plot_on_original(color='white')

plt.show()

'''
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


pos_galactic = '13h11m29.96s -01d19m18.7s'
position = SkyCoord(pos_galactic, frame='icrs')

wcs = WCS(naxis=2)
rho = np.pi / 3.
scale = 0.05 / 3600.

wcs.wcs.cd = [[scale*np.cos(rho), -scale*np.sin(rho)],
              [scale*np.sin(rho), scale*np.cos(rho)]]
wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
wcs.wcs.crval = [position.ra.to(u.deg).value,
                 position.dec.to(u.deg).value]
wcs.wcs.crpix = [50, 100]


cutout = Cutout2D(data, position, (30, 40), wcs=wcs)
plt.imshow(cutout.data, origin='lower')  

plt.show()

'''













