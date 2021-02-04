#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 14:52:46 2021

@author: tomandrews_h
"""

"""
This is a test file used to test the object_find function on a synthetic galaxy and test on small boxes of the image.
When testing on small boxes on the image, it is noted that you have to re run the corresponding cell each time. 
The main purpose of this file is to produce a plot to determine the effective radius of a synthetic galaxy, hence the graph settings
in object_find are tailored to the synthetic galaxy, however they can be easily adjusted to produce a plot of any galaxy if you wish 
to input your own box to test the function on.

"""

from astropy.io import fits
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from astropy.table import Table
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
#%%
params = {
        'axes.labelsize': 15,
        'font.size': 17,
        'font.family': 'sans-serif',
        'font.style': 'normal',
        'font.serif': 'Arial',
        'legend.fontsize': 10,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'figure.figsize':[11,7]}

plt.rcParams.update(params)
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
"""
Determines source magnitude from total pixel count of object
"""
def source_magnitude(pc):
    return -2.5*sp.log10(pc) + 25.3


#%%
"""
Produces aperture around input pixel location and performs adaptaive aperture method
This is a test version of obejct_find and so produces a plot after each iteration, this function is not desigend to be produecd over the whole image
This function detrermies effective radius of a galaxy using a tanh fit, this is not used in the main function as 'curve_fit' is time consuming. 
"""
def object_find(x_centre, y_centre, mask_image_ones):
    radius_list = sp.arange(1,51,1)     #initial radius list to iterate through (aticapte radius's won't be higher than 12)
    radius_plot = []              #list to append radius values              
    sum_list = []               # list to append cumulative sum
    bkg_list =[]                #list to append background values
    err_fs_list =[]    #list to append error in final sum
    bkg_mean_error_list =[]  #list to append bkg mean error
    bkg_diff = -10              #initial value for bkg_diff
    for radius in radius_list:
        if abs(bkg_diff) < 1 or (bkg_diff*(-1))<= 0:        #criteria for adaptive aperture
            break
        else:
            position_centre = (x_centre,y_centre)           #location of centre of aperture
            aperture = CircularAperture(position_centre, r = radius)    #create aperture object
            #aperture.plot(color='white', lw=2)  #overplot the aperture circle onto a colormap plot (need to produce an plot before running function)
            
            aperture2 = CircularAperture(position_centre, r = radius+1)     #create second aperture object to caclulate adjusted annulus area, radius is same size as annulus
            
            annulus_aperture = CircularAnnulus(position_centre, r_in = radius, r_out = radius+1)    #create annulus aperture
            #annulus_aperture.plot(color='red', lw=2) #overplot the annulus circle onto a colormap plo
            
            apers = [aperture, annulus_aperture]
            mask_image_zeros = 1- mask_image_ones               
            a = 1 - mask_image_ones            #invert mask image matrix for photoutils
            mask_global = a.astype(bool)
            
     
            phot_table = aperture_photometry(cropped_image, apers,mask = mask_global)  #counts flux in aperture and annulus, stores in phot_table
            mask1 = aperture.to_mask(method = 'center')     #creates new mask object for aperture
            mask = mask1
            image1 = mask.to_image(shape = cropped_image.shape)#creates corresponding mask matrix of same shape as our image
            
            aperture_area = aperture.area - np.sum(image1*np.copy(mask_image_zeros)) #adjusted apertrue area, removing area of aperture already masked
           
            
            mask2 = aperture2.to_mask(method = 'center')
            mask = mask2
            image2 = mask.to_image(shape= cropped_image.shape)
            annulus_area = (aperture2.area - np.sum(image2*np.copy(mask_image_zeros)))- aperture_area #adjusted annulus, removing area of annulus already masked
           
            
            bkg_mean = phot_table['aperture_sum_1']/annulus_area
            bkg_list.append(bkg_mean)  #append bkg mean for each aperture to list
            bkg_sum = bkg_mean*aperture_area #calculate backgroud contribution for aperture
            
            err_bkg_mean = sp.sqrt(bkg_mean/3.1)    #error in the mean value of bkg
            bkg_mean_error_list.append(sum(err_bkg_mean))   
            err_bkg = sp.sqrt(bkg_sum/3.1)  #error in bkg given by poisson statistics (divided by gain = 3.1)
            
            final_sum = phot_table['aperture_sum_0'] - bkg_sum  #calculate total contribution with bkg removed
            
            err_fs = sp.sqrt((phot_table['aperture_sum_0'])/3.1 + err_bkg**2)  #poission error in final sum, added in quadrature with error from background(gain = 3.1)
            err_fs_list.append(sum(err_fs))  #err_fs is in form of single list, so sum just takes absolute value and appends to main list
            phot_table['residual_aperture_sum'] = final_sum  #update phototable with aperture sum
            sum_list.append(final_sum)
            radius_plot.append(radius)
            
            if radius == 1:
                continue
        
            else:
                bkg_diff = bkg_list[radius-1] -bkg_list[radius-2]
                continue
                
            continue
    
    masks = aperture.to_mask(method ='center')  #creates mask of final aperture
    mask = masks
    image_l = mask.to_image(shape = cropped_image.shape)  
    b  = 1 - image_l   #invert mask matrix
    
    mask_image_ones = mask_image_ones*b  #add new masked object to global mask
    
    
    def tanhh(x, a, b,c):  #define tanh to fit to cumulative frequency plot
        return c*sp.tanh((x-a)/b)
    
    if len(radius_plot) >3:     # avoids catalogueing of high background noise
        popt,pcov = curve_fit(tanhh,sp.ravel(radius_plot),sp.ravel(sum_list),p0=[6, 0.5,5000])  #returns tanh fit parameters 
        r_e = sp.sqrt(sp.diag(pcov))   #error covariance matrix
        radius = np.arctanh(0.5)*abs(popt[1])+abs(popt[0])  #interpolate tanh function to find effective radius
        radius_err = sp.sqrt((np.arctanh(0.5)*r_e[1])**2 + r_e[0]**2) #standard error propogation to find error in radius
        
    else:
        radius = None
        
    
    """
    produces plots within object_find
    following graph is designed for synethetic 2d gaussian, however parameters can be adjustedfor any galaxy
    """
    f1, ax1 = plt.subplots()
    plt.title('Galaxy plot', fontweight = 'bold')
    ax1.set_xlabel('Pixel Radius')
    ax1.set_ylabel('Cumulative Pixel Count')
    ax1.errorbar(radius_plot,sum_list,yerr = err_fs_list, marker ='x',  label = 'Cumulative Frequency', linestyle = ' ')  #errorbars for final sum
    #print(bkg_list, 'bkg list')
    
    x_values = sp.arange(1,max(radius_plot),0.1)
    
    ax1.plot(x_values,tanhh(x_values, *popt), label="Tanh Fit", color='darkslategrey')  #plot tanh fit
    ax1.vlines(radius, 0, tanhh(radius, *popt), linestyles='dashed', label='Effective radius = {:.2f}{}{:.2f}'.format(radius,u'\u00b1', radius_err), color='darkred')
    ax1.hlines(abs(popt[2]/2),0, radius, linestyles='dashed', label='Half of Total Flux = {:.2e}'.format(popt[2]/2), color='darkcyan' )
    
    plt.xlim([0,14])
    plt.ylim([0,600000])
    ax3 = ax1.twinx()
    plt.ylim([3000,13000])
    ax3.set_ylabel('Background mean')
    ax3.errorbar(radius_plot,bkg_list,yerr = sp.array(bkg_mean_error_list),marker = 'x', color = 'red', label = 'Sky Background', linestyle = ' ') #errorbars for bkg
    ax1.xaxis.grid(True, ls = '--', which= 'minor')
    ax1.yaxis.grid(True, ls = '--', which= 'minor')
    ax1.xaxis.grid(True, ls = '-', which= 'major')
    ax1.yaxis.grid(True, ls = '-', which= 'major')
    ax3.xaxis.grid(True, ls = '-', which= 'major')
    ax3.yaxis.grid(True, ls = '-', which= 'major')
    ax3.set_yticks(np.linspace(ax3.get_yticks()[0], ax3.get_yticks()[-1], len(ax1.get_yticks())))
    ax1.xaxis.set_major_locator(MultipleLocator(2))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax1.xaxis.set_minor_locator(MultipleLocator(1))
    ax1.yaxis.set_major_locator(MultipleLocator(100000))
    ax1.yaxis.set_minor_locator(MultipleLocator(50000))
    ax3.yaxis.set_major_locator(MultipleLocator(2000))
    ax3.yaxis.set_minor_locator(MultipleLocator(1000))
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax3.get_legend_handles_labels()
    
    f1.legend(handles1+handles2, labels1+labels2, fancybox = True, bbox_to_anchor=(0.9, 0.5),loc = 'center right')
    f1.savefig('galaxy_plot3')
    
    return radius, final_sum, mask_image_ones
#%%

""" Make a 2d gaussian plot

    size is the length of a side of the square
    sigma is sd of gaussian
    
    central_maximum is a scaling so intensity is
    approximately equal to that of a galaxy
    """
    
def twodGaussian(size, sigma ,central_max, center):

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return central_max*np.exp(-((x-x0)**2 + (y-y0)**2) / sigma**2) + 3418

cropped_image = twodGaussian(size=30, sigma = 4, central_max = 10000, center = None)  #creates image of synthetic data
mask_image_ones = sp.ones(cropped_image.shape, dtype = int)  #create mask image

x_min = 0 #boundaries of 2d gaussian image
y_min = 0
x_max = 30
y_max = 30
#%%
"""
create columns for phototable t
"""
x= []
y = []
pxc =[]
err_bkg =[]
eff_r = []
error = []
mag = []
#%%
"""

"""
plt.imshow(cropped_image, interpolation='nearest', origin = 'lower')  #produces color map of gaussian
y_centre, x_centre = max_brightness(cropped_image, mask_image_ones)
#x_centre = 353
#y_centre = 3037
radius, final_sum, mask_image_ones = object_find(x_centre, y_centre, mask_image_ones)
if radius != None:     #skips if invalid object
    x.append(x_centre)   #appends data to phototable
    y.append(y_centre)
    pxc.append(final_sum)
    eff_r.append(radius)
    error.append(sp.sqrt(final_sum))
    mag.append(source_magnitude(final_sum))
    t = Table([x,y,pxc, eff_r, error, mag], names = ('x center', 'y center', 'pixel count', 'effective radius', 'error', 'source magnitude'))
    y_centre, x_centre = max_brightness(cropped_image, mask_image_ones)   #finds next brightest image pixel
    
print(y_centre,x_centre)
f1 = plt.figure()
f1.add_subplot(111)
plt.title('Masked Image')
plt.spy(mask_image_ones, origin = 'lower')  #prdouces final mask for gaussian 
#%%
"""
Set boundaries to test object_find on specific boxes of data
"""

hdulist = fits.open('A1_mosaic.fits')

#box containing two close galaxies
x_min = 389
y_min = 1550

x_max = 451
y_max = 1612

#%%
"""
Clean up messy/ bleeding objects with updated approach
"""

image = hdulist[0].data
mask_image_ones = sp.ones(image.shape, dtype = int)       #create a mask image of all 1s

clean_image_new(0, y_min, 0, 2570)  #remove boundaries from image, it is sim
clean_image_new(y_max, 4611, 0, 2570)
clean_image_new(0, 4611, 0, x_min)
clean_image_new(0, 4611, x_max, 2570)

"""
creates circular mask object and masks this from image
"""
def circle_mask(position_centre, radius, mask_image, image_input):  
    aperture = CircularAperture(position_centre, r = radius)
    masks = aperture.to_mask(method ='center')
    
    mask = masks
    big_circle = mask.to_image(shape = image_input.shape) 

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
crop image after manually masking
"""

cropped_image_x1 = sp.delete(hdulist[0].data, slice(x_max,2571),1)
cropped_image_x2 = sp.delete(cropped_image_x1, slice(0,x_min),1)

cropped_image_y1 = sp.delete(cropped_image_x2, slice(y_max, 4612), 0)
cropped_image = sp.delete(cropped_image_y1, slice(0,y_min), 0)

masked_image_x1 = sp.delete(mask_image_ones, slice(x_max,2571),1)
masked_image_x2 = sp.delete(masked_image_x1, slice(0,x_min),1)

masked_image_y1 = sp.delete(masked_image_x2, slice(y_max, 4612), 0)
mask_image_ones = sp.delete(masked_image_y1, slice(0,y_min), 0)

plt.title('Masked Image')
plt.spy(mask_image_ones, origin = 'lower')  #displays the initial mask image before beginning identifying objects

"""
create columns for table t
"""
x= []
y = []
pxc =[]
err_bkg =[]
eff_r = []
error = []
mag = []

f2 = plt.figure()
ax21 = f2.add_subplot(111)
ax21.imshow(cropped_image, interpolation='nearest', origin = 'lower')  #color plot of cropped ccd image

y_centre, x_centre = max_brightness(cropped_image, mask_image_ones)
#%%
"""
this cell can be run for however many times you want, it is  manual version of the while loop in the main file
for example with the current image boundaries(below) this cell can be run twice to show the photometry method working 
x_min = 389
y_min = 1550

x_max = 451
y_max = 1612

note that the graph settings for plots produced in object_find are tailored to the 2d gaussian, so you will need to adjust the ssettings to get more apealing plots
"""

radius, final_sum, mask_image_ones = object_find(x_centre, y_centre, mask_image_ones)
if radius != None:
    x.append(x_centre)
    y.append(y_centre)
    pxc.append(final_sum)
    eff_r.append(radius)
    error.append(sp.sqrt(final_sum))
    mag.append(source_magnitude(final_sum))
    t = Table([x,y,pxc, eff_r, error, mag], names = ('x center', 'y center', 'pixel count', 'effective radius', 'error', 'source magnitude'))
    y_centre, x_centre = max_brightness(cropped_image, mask_image_ones)
    

#%%    
"""
simple plot of final mask image, interesting to compare with imshow plot 
"""
f3 = plt.figure()
ax3 = f3.add_subplot(111)
ax3.spy(mask_image_ones, origin = 'lower')
ax3.set_title('Adjacent Objects', fontweight = 'bold')  
