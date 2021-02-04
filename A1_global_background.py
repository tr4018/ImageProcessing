#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: tomandrews_h
"""
"""
This file contains code to plot a histogram for the upper bound of the image
corresponding to oversaturation that needs to be masked

The second part of this file determines the global background noise by fitting
a gassian to a histogram spread of the noise.

"""
from astropy.io import fits
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from photutils import CircularAperture
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
                               


#%%%
"""
Plots histogram for upper bound of image
"""
hdulist = fits.open('A1_mosaic.fits')
x_min = 263
y_min = 243

x_max = 2350
y_max = 4400
data_list = []

for y in range(y_min, y_max):
    for x in range(x_min,x_max):
        if hdulist[0].data[y][x] > 25000:  #filters data for pixels with count above 25000
            data_list.append(hdulist[0].data[y][x])
        else:
            continue

fig = plt.figure()
ax = fig.add_subplot(111) 
frequency,count, a = ax.hist(data_list, bins = 500, color = 'steelblue')  #plot histogram
 
ax.set_xlabel('Pixel Count', size='15')
ax.set_ylabel('Frequency', size='15')
ax.set_title('Image Upper Bound', fontweight = 'bold', size = '15')

ax.xaxis.grid(True, ls = '--', which= 'minor')
ax.yaxis.grid(True, ls = '--', which= 'minor')
ax.xaxis.grid(True, ls = '-', which= 'major')
ax.yaxis.grid(True, ls = '-', which= 'major')
ax.xaxis.set_major_locator(MultipleLocator(2000))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
ax.xaxis.set_minor_locator(MultipleLocator(1000))
ax.yaxis.set_major_locator(MultipleLocator(200))
ax.yaxis.set_minor_locator(MultipleLocator(100))
       
fig.savefig('UpperBoundtest')
plt.show()
#%%
"""
Identify brightest stars/galaxies
"""

def max_brightness(image_input, mask_image_ones):
    image_high = image_input * mask_image_ones  #multiply mask image by CCD Image
    ind = np.unravel_index(np.argmax(image_high, axis=None), image_input.shape)   #returns location of brightest pixel
    return ind
        
"""
Move 0s to matrix mask_image_ones to mask off particular region
"""
def clean_image_new(y1, y2, x1, x2):    #input boundaries of region to mask
    for y in range(y1, y2):
        for x in range(x1,x2):
            if image[y][x] >= 0:
                mask_image_ones[y][x] = 0
            else:
                continue
#%%
"""
Clean up messy/ bleeding objects with updated approach
"""

image = hdulist[0].data
mask_image_ones = sp.ones(image.shape, dtype = int)       #create a mask image of all 1s

clean_image_new(0, y_min, 0, 2570)  #remove boundaries from image
clean_image_new(y_max, 4611, 0, 2570)
clean_image_new(0, 4611, 0, x_min)
clean_image_new(0, 4611, x_max, 2570)

"""
creates circular mask object of gievn radius and masks this from image
requires the global mask image as an input to remove then add the mask object to this image
"""
def circle_mask(position_centre, radius, mask_image, image_input):   
    aperture = CircularAperture(position_centre, r = radius)  #creates an aperture about the input pixel location. 
    masks = aperture.to_mask(method ='center') #creates mask obejct
    
    mask = masks
    big_circle = mask.to_image(shape = image_input.shape)  #creates mask matrix of same shape as input image

    central_circle_mask = 1 - big_circle  #invert matrix
    
    mask_image_ones = mask_image*central_circle_mask  #add mask to global mask
    return mask_image_ones


mask_image_ones = circle_mask(position_centre = (1434, 3199), radius = 290., mask_image = mask_image_ones,image_input = image)



for y in range(106, 450):   #remove gaussians at bottom of image
    for x in range(1014, 1692):
        if hdulist[0].data[y][x] > 6000:
            mask_image_ones[y][x] = 0
        else:
            continue

for y in range(449, 489):     #remove gaussians at bottom of image
    for x in range(1370,1470):
        if hdulist[0].data[y][x] > 6000:
            mask_image_ones[y][x] = 0
        else:
            continue
        
clean_image_new(2221,2363, 860, 950)   #manual masking of image
clean_image_new(3199, 3421, 724, 836)
clean_image_new(0, 4571, 1422, 1454)
clean_image_new(2700,2846, 930, 1024)
clean_image_new(4052, 4140, 513, 609)
clean_image_new(1398, 1454, 2056, 2138)
clean_image_new(3996, 4070, 1429,1495)
clean_image_new(3704,3806, 2094, 2175)
clean_image_new(2276,2338, 2097, 2165)
mask_image_ones = circle_mask(position_centre = (2466,3410), radius = 34., mask_image = mask_image_ones, image_input = image)

#%%
"""
Plots histogram fior background noise of entire image
"""
hdulist = fits.open('A1_mosaic.fits')

data_list = []

for y in range(y_min, y_max):
    for x in range(x_min,x_max):
        if hdulist[0].data[y][x]*mask_image_ones[y][x] < 3500 and hdulist[0].data[y][x]*mask_image_ones[y][x] > 0: #ignore data that isn't in the region of the background noise
            data_list.append(hdulist[0].data[y][x])
        else:
            continue
#%%

fig = plt.figure()
ax = fig.add_subplot(111) 
frequency,count, a = ax.hist(data_list, bins = 155, color = 'steelblue')
 
bin_middle = []
for i in range(1,len(count)):
    x = (count[i]+count[i-1])/2 # calculates centre of each bin
    bin_middle.append(x) #appendsto a list that can be used to estimate the parameters of the gaussian fit
 
def gaus(x,b,c,d):
    return d*sp.exp(-(x-b)**2/(2*c**2))

popt,pcov = curve_fit(gaus,bin_middle,frequency,p0=[3500, 10, 3000000])   #returns gaus paramaters as defined in function gaus

g = sp.arange(min(bin_middle),3500,0.5)

ax.plot(g,gaus(g, *popt), label="Gaussian Fit", color='darkred')  #plots gauss
ax.vlines(popt[0], 0, popt[2], linestyles='dashed', label= '{} = 3418.96{}0.02'.format(u'\u03bc',u'\u00b1'), color='darkslategray')
ax.hlines(gaus(popt[0]+popt[1],*popt),popt[0], popt[0]-popt[1], linestyles='solid' , color = 'darkslategray')
ax.set_xlabel('Pixel Count', size='15')
ax.set_ylabel('Frequency', size='15')
ax.set_title('Histogram of Background Noise', fontweight = 'bold', size = '15')
ax.text(popt[0]+abs(popt[1])+2, gaus(popt[0]+popt[1],*popt)-10000,'{}=11.39{}0.02'.format(u'\u03C3',u'\u00b1'), fontsize = 16, color = 'darkslategray')
ax.xaxis.grid(True, ls = '--', which= 'minor')
ax.yaxis.grid(True, ls = '--', which= 'minor')
ax.xaxis.grid(True, ls = '-', which= 'major')
ax.yaxis.grid(True, ls = '-', which= 'major')
ax.xaxis.set_major_locator(MultipleLocator(20))

ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(50000))
ax.yaxis.set_minor_locator(MultipleLocator(25000))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))

ax.arrow(popt[0]-popt[1]-8,gaus(popt[0]+popt[1],*popt),5,0,fc="darkslategray", ec="darkslategray",
head_width=5000, head_length=1, overhang = 3)
ax.arrow(popt[0]+8,gaus(popt[0]+popt[1],*popt),-5,0,fc="darkslategray", ec="darkslategray",
head_width=5000, head_length=1, overhang = 3)
plt.legend(fontsize = 'medium')

print(popt[0], popt[1], popt[2])  
print(sp.sqrt(sp.diag(pcov)))        
plt.savefig('Histogram of Background Noisetest')
plt.show()