import pyfits
import matplotlib.pyplot as plt


f = pyfits.open("im1.fits")
h, d = f[0].header, f[0].data



import pywcsgrid2


ax = pywcsgrid2.axes(header=h)

ax.set_display_coord_system("gal")
ax.set_ticklabel_type("absdeg", "absdeg")

axis = ax['gal'].new_floating_axis(1, 0.)
ax.axis['b=0'] = axis

plt.imshow(d, origin="lower")
ax.grid(alpha = 1, linestyle = 'solid', color = "w")


plt.show()



