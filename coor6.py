import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

filename = get_pkg_data_filename('im1.fits')

hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)

ax = plt.subplot(projection=wcs)

ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')

ax.coords[0].set_axislabel('Galactic Longitude')
ax.coords[1].set_axislabel('Galactic Latitude')
ax.coords.grid(color='yellow', ls='solid', alpha = 1)

overlay = ax.get_coords_overlay('fk5')
overlay.grid(color='white', ls='solid', alpha = 1)
overlay[0].set_axislabel('Right Ascension (J2000)')
overlay[1].set_axislabel('Declination (J2000)')

plt.subplot(projection=wcs)
plt.imshow(hdu.data, origin='lower')
#plt.grid(color='white', ls='solid')
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')

plt.show()

'''
f = aplpy.FITSFigure('im1.fits')
f.show_colorscale(vmin=0)
f.show_contour('im1.fits',colors='white')
'''

'''
from astropy.wcs import WCS
from astropy.io import fits

data = fits.getdata('im1.fits')
wcs = WCS('im1.fits')


fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, projection=wcs)
ax.get_transform("galactic")
ax.imshow(data, vmin=0, vmax=1.e-4, origin='lower')
ax.grid(color='white', alpha=1, ls='solid')
ax.set_xlabel("Galactic Longitude")
ax.set_ylabel("Galactic Latitude")
overlay = ax.get_coords_overlay('fk5')
plt.show()

'''
'''

fig1 = plt.figure(figsize=(8,8))
ax = fig1.add_axes([0.25, 0.25, 0.6, 0.6], projection=wcs)
ax.imshow(hdu, vmin=0, vmax=1.e-4, origin='lower', cmap=plt.cm.gist_heat)

overlay = ax.get_coords_overlay('fk5')

ax.coords['glon'].set_ticks(color='white')
ax.coords['glat'].set_ticks(color='white')

ax.coords['glon'].set_axislabel('Galactic Longitude')
ax.coords['glat'].set_axislabel('Galactic Latitude')

ax.coords.grid(color='yellow', linestyle='solid', alpha=1)

overlay['ra'].set_ticks(color='white')
overlay['dec'].set_ticks(color='white')

overlay['ra'].set_axislabel('Right Ascension')
overlay['dec'].set_axislabel('Declination')

overlay.grid(color='white', linestyle='solid', alpha=1)

plt.show()
'''