import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

# Read in file
hdulist = fits.open('im1.fits')

# Extract image and header
image = hdulist[0].data
header = hdulist[0].header

# Instantiate WCS object
w = WCS(header)

# Plot the image
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.imshow(image, cmap=plt.cm.gist_heat,
          origin='lower', vmin=0, vmax=1000.)

# Loop over lines of longitude
for lon in np.linspace(-180., 180., 13):
    grid_lon = np.repeat(lon, 100)
    grid_lat = np.linspace(-90., 90., 100)
    px, py = w.wcs_world2pix(grid_lon, grid_lat, 1)
    ax.plot(px, py, color='white', alpha=0.5)

# Loop over lines of latitude
for lat in np.linspace(-60., 60., 5):
    grid_lon = np.linspace(-180., 180., 100)
    grid_lat = np.repeat(lat, 100)
    px, py = w.wcs_world2pix(grid_lon, grid_lat, 1)
    ax.plot(px, py, color='white', alpha=0.5)

ax.set_xlim(0, image.shape[1])
ax.set_ylim(0, image.shape[0])
ax.set_xticklabels('')
ax.set_yticklabels('')
fig.savefig('wcs_extra.png', bbox_inches='tight')

plt.show()