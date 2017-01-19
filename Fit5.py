import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar
from astropy.io import fits
from math import exp
from scipy.interpolate import UnivariateSpline
import math 

Wise1 = fits.open('4imAF_5.fits')# im1\Cut2image4.fits
im1 = Wise1[0].data

b = im1.ravel()

b[np.isnan(b)] = 0
b[np.isinf(b)] = 0


#------------------PART I: Histogram----------------------------------

#Creates the Histogram (hist, bin_edges)
#             bin_edges is the x-axis
#             hist is the y-axis
hist, bin_edges = np.histogram(b, bins = range(260,339))#(260,350))#(1000,1220))

#plots the histogram
plt.bar(bin_edges[:-1], hist, width = 1)  
plt.xlim(min(bin_edges), max(bin_edges))


#For the curve fitting, we need to use the centers of each bin for a more accurate fit, which we define here as bin_centres
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2


hist2, bin_edges2 = np.histogram(b, bins = range(260,600))#(260,350))#(1000,1220))

#plots the histogram
plt.bar(bin_edges2[:-1], hist2, width = 1)  
plt.xlim(min(bin_edges2), max(bin_edges2))


#For the curve fitting, we need to use the centers of each bin for a more accurate fit, which we define here as bin_centres
bin_centres2 = (bin_edges2[:-1] + bin_edges2[1:])/2




#------------------PART II: Fitting Function--------------------------

#Defining the function that will fit the histogram data
def f(x, p1, p2, p3):
    
    return p1*np.exp((-np.absolute(((x-p2)/p3))**2)*np.log(2))
    #return p3*(p1/((x-p2)**2 + (p1/2)**2)) 
    #return p3*(p1/((x-p2)**2 + (p1/2)**2)) 
    #return p1*exp(-(np.log(x - p2)/p3)**(3 + np.tanh(4*np.sqrt((x-p2)**2)/p3 - 3))*np.log(2))
    #return p1*exp(-(np.log(x - p2)/p3)**(3 + np.tanh(4*(x-p2)*np.sign(x-p2))/p3 - 3))*np.log(2))
          #p1*exp(-(|x - p2|/p3)**(3+tanh(4*|x-p2|/p3 - 3))*ln2) #original function 



# Making the actual fit by the function curvefit(f, x, y, guesspoints)

#guesspoints
p0 = (80, 400, 200)

#popt, pcov = curve_fit(f, xdata, hist)
popt, pcov = curve_fit(f, bin_centres, hist, p0=p0)

#popt is a 1D array
#pcov is a 2D array


#------------------PART III: Finding the Highest point--------------------------

# This part of the code simply finds the peak of the fitted function
fm = lambda x: -f(x, *popt)
r = minimize_scalar(fm, bounds=(1, 5))

print "maximum:", r["x"], f(r["x"], *popt)  #example ~ maximum: 2.99846874275 18.3928199902


#------------------PART IV: Plotting Everything--------------------------

#Creates 'blue' outline of the histogram data
#plt.plot(bin_centres, hist, label='Test data')

#The fitted function is outlined 'green' 
##hist, bin_edges = np.histogram(b, bins = range(863,894))#(260,350))#(1000,1220))
hist_fit = f(bin_centres, *popt)
fit1 = plt.plot(bin_centres, hist_fit, label='Fitted data')
plt.setp(fit1, color='g')

#Graphically marks the highest point with a 'black' dot
x_curve = np.linspace(1, 5, 100)
plt.plot(x_curve, f(x_curve, *popt))
plt.plot(r['x'], f(r['x'], *popt), 'ko')

#Prints axis labels
plt.title("Fitted Histogram 4imAF_5")
plt.xlabel("Pixel Intensity (DN)")
plt.ylabel("Number of Pixels with Intensity Value")

#plt.imshow(im1, origin='lower') #2D plot -> not useful at the moment

#Command to show all plots
plt.show()