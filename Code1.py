import pyfits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from astropy.io import fits

#load the first image
data1 = pyfits.getdata('4imAF.fits') 

#set a noise floor for the image. Poor SNR pixels will give you bad temperature measurements.
floor1 = 0.045

#set low SNR pixels to nans
av1 = np.where(data1<floor1)
data1[av1] = np.nan

#load the second image
data2 = pyfits.getdata('3imAF.fits') 

#set a noise floor for the image. Poor SNR pixels will give you bad temperature measurements.
floor2 = 0.070

#set low SNR pixels to nans
av2 = np.where(data2 < floor2)
data2[av2] = np.nan

#Color temperatures

#wavelength of first image in microns
l1 = 12e-4
#extinction factor for the first image (if you haven't already applied an extinction correction)
ext1 = 2.58

#wavelength of second image in microns
l2 = 22e-4
#extinction factor for the second image
ext2 = 1.54

#set constants (in cgs)
c = 3e10
h = 6.626e-27
c = 2.9979e10
kb = 1.38e-16

#convert wavelengths of images to frequencies
v1 = c / l1
v2 = c / l2

#define function for b_nu that we will solve for T
func = lambda T : (v1/v2)**5*(np.expm1((h*v2)/(kb*T)))/(np.expm1((h*v1)/(kb*T)))-I_1/I_2

#provide an initial guess for T for the solver to use. 
Tinit = 500.0

#get the shape of the images so we know how big to make the for loop
shp = np.shape(data1)

#create a new empty array to store the temperature values
Tim = np.zeros((shp[0],shp[1]))

for i  in range (0,shp[0]):
    for j in range (0,shp[1]):

        #store the pixel values from the images
        I_1 = data1[i,j]
        I_2 = data2[i,j]
        
        #apply the extinction factors
        I_1 = I_1- 312
        I_2 = I_2 -1204
    
        #use the solver to find the temperatures
        TS = fsolve(func,Tinit)
    
        #store the values in the array
        Tim[i,j] = TS

#Where the solver returned Tinit, turn into a nan. 
Tim[Tim == Tinit] = np.nan

#plot the color temperature map
plt.figure()
plt.title('Color-Tempature Map')
plt.imshow(Tim[::-1])
#plt.clim(lowval, highval)
plt.colorbar()
plt.show()



#fits.writeto('out.fits', data, header)
