
from astropy import wcs
from astropy import coordinates
from astropy import units as u

mywcs = wcs.WCS(header)

pix_coords = [(20,30), (50,10), (500,200)]
world_coords = mywcs.wcs_pix2world(pix_coords, 0)
print(world_coords)


my_coords = coordinates.SkyCoord(world_coords*u.deg, frame='fk5')
print(my_coords)
print(my_coords.to_string(style='hmsdms'))

