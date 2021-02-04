#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: tomandrews_h

"""
"""
This is the main program fully developed after testing and calculation of global abckground noise (3464)
The code will run over the image and identify obekcts in order of brightness, after each valid object is identfied and catlogued, 
it is then masked  so as not to be rediscovered
All analysis is done in the code and plots of the galaxy number count will be produced and saved at the end
Each cell can be run sequentially(or run code all at once). The loop over the image takes approx 3 hrs to run.
Data will save directly into your repository and then read from your repository before producing plots.

"""
from astropy.io import fits
import scipy as sp
import matplotlib.pyplot as plt
import numpy as np
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from astropy.table import Table
import pandas as pd
from matplotlib.ticker import (MultipleLocator)
                               
#%%
"""
Mask edges to avoid edge effects
"""

hdulist = fits.open('A1_mosaic.fits')

x_min = 263
y_min = 243

x_max = 2350
y_max = 4400

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

"""
Produces aperture around input pixel location and performs adaptaive aperture method
"""
def object_find(x_centre, y_centre, mask_image_ones):
    radius_list = sp.arange(1,51,1)     #initial radius list to iterate through (aticapte radius's won't be higher than 12)
    radius_plot = []              #list to append radius values              
    sum_list = []               # list to append cumulative sum
    bkg_list =[]                #list to append background values
    bkg_diff = -10              #initial value for bkg_diff
    for radius in radius_list:
        if abs(bkg_diff) < 1 or (bkg_diff*(-1))<= 0:        #criteria for adaptive aperture
            break
        else:
            position_centre = (x_centre,y_centre)           #location of centre of aperture
            aperture = CircularAperture(position_centre, r = radius)    #create aperture object
            #aperture.plot(color='white', lw=2)
            
            aperture2 = CircularAperture(position_centre, r = radius+1)     #create second aperture object to caclulate adjusted annulus area, radius is same size as annulus
            
            annulus_aperture = CircularAnnulus(position_centre, r_in = radius, r_out = radius+1)    #create annulus aperture
            #annulus_aperture.plot(color='red', lw=2)
            
            apers = [aperture, annulus_aperture]
            mask_image_zeros = 1- mask_image_ones               
            a = 1 - mask_image_ones            #invert mask image matrix for photoutils
            mask_global = a.astype(bool)
            
     
            phot_table = aperture_photometry(image, apers,mask = mask_global)  #counts flux in aperture and annulus, stores in phot_table
            mask1 = aperture.to_mask(method = 'center')     #creates new mask object for aperture
            mask = mask1
            image1 = mask.to_image(shape = image.shape)#creates corresponding mask matrix of same shape as our image
            
            aperture_area = aperture.area - np.sum(image1*np.copy(mask_image_zeros)) #adjusted apertrue area, removing area of aperture already masked
           
            
            mask2 = aperture2.to_mask(method = 'center')
            mask = mask2
            image2 = mask.to_image(shape= image.shape)
            annulus_area = (aperture2.area - np.sum(image2*np.copy(mask_image_zeros)))- aperture_area #adjusted annulus, removing area of annulus already masked
            
            
            bkg_mean = phot_table['aperture_sum_1']/annulus_area 
            bkg_list.append(bkg_mean)   #append bkg mean for each aperture to list
            bkg_sum = bkg_mean*aperture_area  #caclulate backgroud contribution for aperture
            final_sum = phot_table['aperture_sum_0'] - bkg_sum      #calculate total contribution with bkg removed
            err_fs = sp.sqrt((final_sum/3.1)+(bkg_sum/3.1)) #poission error in final sum, added in quadrature with error from background(gain = 3.1)
            phot_table['residual_aperture_sum'] = final_sum  #update phototable with aperture sum
            sum_list.append(final_sum)
            radius_plot.append(radius)
           
            
            if radius == 1:   #only calculated difference in background between two annuli if radius greater than 1
                continue
            else:
                bkg_diff = bkg_list[radius-1] -bkg_list[radius-2]
                continue
                
            continue
           
    masks = aperture.to_mask(method ='center')  #creates mask of final aperture
    mask = masks
    image_l = mask.to_image(shape = image.shape)  
    b  = 1 - image_l   #invert mask matrix
    
    mask_image_ones = mask_image_ones*b  #add new masked object to global mask so identified object not rediscovered
    
    if len(radius_plot) <3:  #avoid cataloguing small objects
        radius = None
        
    return radius, final_sum, err_fs, mask_image_ones
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


plt.title('Mask Image Initial')
plt.spy(mask_image_ones, origin = 'lower')
plt.savefig('clean_image_test')
plt.show()
#%%
"""
Calculates centre and determines radius and brightness
This code runs over enire image until criteria for min brightness is met
"""

x= []   #creates columns for table t
y = []
pxc =[]  #pixel count
ap_r = []  #final apertrue radius
error = []  #error in magnitudew
mag = []        #magnitude

y_centre, x_centre = max_brightness(image, mask_image_ones)  #give initial brighetst pixel for object_find
i=1                     #initial count
while image[y_centre][x_centre] > 3464.6:   #lower bound set by global background
    radius, final_sum, err_fs, mask_image_ones = object_find(x_centre, y_centre, mask_image_ones)  #returns values from object_find
    if radius != None and final_sum >0:   #provided this is a valid object, stores relevant information
        print('object no. {i}  location {location}  pc {pc}'.format(i =i, location = (x_centre, y_centre), pc = image[y_centre][x_centre]))
        x.append(x_centre)
        y.append(y_centre)
        pxc.append(final_sum)
        ap_r.append(radius)
        error.append(sp.sqrt(((2.5/(sp.log(10)*final_sum))*err_fs)**2 + (0.02**2)))
        mag.append(source_magnitude(final_sum))
        t = Table([x,y,pxc, ap_r, error, mag], names = ('x center', 'y center', 'pixel count', 'effective radius', 'error', 'source magnitude'))
        y_centre, x_centre = max_brightness(image, mask_image_ones)  #returns next brightest pixel(now that previous object has been masked)
        i += 1
        continue
    else:
        y_centre, x_centre = max_brightness(image, mask_image_ones)  #if invalid object, move on 
        continue

f2 = plt.figure()
f2.add_subplot(111)
plt.title('Final Mask Image')
plt.spy(mask_image_ones, origin = 'lower')
plt.savefig('final_mask_test')
#%%%
from matplotlib.lines import Line2D

def logN(set_mag, mag):  #determines number of galaxies with mag less than a set magnitude e.g if set_mag = 13, then will return the numebr of galaxies with mag less than 13
    count = 0
    for m in mag:
        if m < set_mag:
            count += 1
            print(m, 'm')
            continue
        else:
            continue
    
    return count

data = []

for mag_input in range(0, int(max(mag))+2):   #range of mag_input defiend by range of magnitude collected
    count =  logN(mag_input, mag)           #for each mag_input, the number of magnitudes less than mag_input will be recorded
    if count == 0:
        continue
    else:
        print(count, 'count')
        data.append([mag_input, sp.log10(count)])  #record the mag_input and the log of the count

m_list =[]
logN_list =[]
m_listall = []
logN_listall =[]
for m, logN in data:   #isolates data for magnitudes less than 17, easier when plotting
    if m<=17:   
        m_list.append(m)
        logN_list.append(logN)
        m_listall.append(m)
        logN_listall.append(logN)
        continue
    else:
        m_listall.append(m)
        logN_listall.append(logN)
        continue
#%%
"""
store collected data
"""
df_phototable = pd.DataFrame({'x': x,
                   'y': y,
                   'pxc':   pxc,
                   'eff_r': eff_r,
                   'error':error,
                   'mag':mag})
df_magdata = pd.DataFrame({'m_list':m_list,
                   'logN_list': logN_list})
df_magdata_all = pd.DataFrame({'m_listall':m_listall,
                   'logN_listall':logN_listall})
    
df_phototable.to_csv('phototable_df5.csv', index=False)
df_magdata.to_csv('magdata_df5.csv', index = False)
df_magdata_all.to_csv('magdata_all_df5.csv', index = False)
df_readme = pd.DataFrame({'data set 5': 'min pc 3464.6 with only saturated objects masked; see image "saturated mask locations" with all data',
                          'x':  (x_min,x_max),
                          'y':  (y_min,y_max),
                          })
df_readme.to_csv('readme_df5.csv', index = False)
#%%
"""
Read data files
"""

magnitude_data = pd.read_csv('magdata_df5.csv')
m_list = magnitude_data['m_list']
logN_list = magnitude_data['logN_list']

magnitude_data_all = pd.read_csv('magdata_all_df5.csv')
m_listall = magnitude_data_all['m_listall']
logN_listall = magnitude_data_all['logN_listall']
#%%
"""
Produce galaxy magnitude count plot
"""
f3 = plt.figure()
ax3 =f3.add_subplot(111)
yerr = 1/(sp.log(10)*sp.sqrt(10**(logN_list))) #error in LogN (Huang et al, 1946)
fit, cov = np.polyfit(m_list, logN_list, 1,w = 1/yerr, cov = True)  #returns fit parameters from data
poly = sp.poly1d(fit)

yerr2 = 1/(sp.log(10)*sp.sqrt(10**(logN_listall[7:])))  #error in LogN (Huang et al, 1946)
fit2, cov2 = np.polyfit(m_listall[7:], logN_listall[7:], 1,w = 1/yerr2, cov = True)  #returns fit parameters from data
poly2 = sp.poly1d(fit2)

#produce two serpeate linear fits for the faint and bright region
ax3.plot(m_listall[:9],poly(m_listall[:9]), color = 'darkred', linestyle = '--',label = 'bright gradient: {fit:.3f} +/- {error:.3f}'.format(fit = fit[0], error =sp.sqrt(cov[0][0])))
ax3.plot(m_listall[7:],poly2(m_listall[7:]), color = 'lightsalmon', linestyle = '--',label = 'faint gradient: {fit:.3f} +/- {error:.3f}'.format(fit = fit2[0], error =sp.sqrt(cov2[0][0])))

ax3.errorbar(m_listall, logN_listall, yerr = 1/(sp.log(10)*sp.sqrt(10**(logN_listall))), color = 'darkslategray',marker = 'x', linestyle ='None')
ax3.set_xlabel('Magnitude/ mag')
ax3.set_ylabel('Log [N(m)/ $mag^{-1}$$deg^{-2}$]')
ax3.set_title ('Galaxy Magnitude Plot', fontweight = 'bold')

ax3.xaxis.grid(True, ls = '--', which= 'minor')
ax3.yaxis.grid(True, ls = '--', which= 'minor')
ax3.xaxis.grid(True, ls = '-', which= 'major')
ax3.yaxis.grid(True, ls = '-', which= 'major')

ax3.xaxis.set_major_locator(MultipleLocator(2.5))
ax3.xaxis.set_minor_locator(MultipleLocator(1.25))
ax3.yaxis.set_major_locator(MultipleLocator(0.5))
ax3.yaxis.set_minor_locator(MultipleLocator(0.25))
plt.ylim([0.75,3.625])
plt.xlim([9,23.5])
print('gradient: ', fit[0], ' +/-', sp.sqrt(cov[0][0]))
print('y_intercept: ', fit[1], ' +/-', sp.sqrt(cov[1][1]))
ax3.legend(loc = 'lower right')
f3.savefig('logN_plot_df5_3')
plt.show()
#%%
"""
crop image after analysis done
"""

cropped_image_x1 = sp.delete(hdulist[0].data, slice(x_max,2571),1)
cropped_image_x2 = sp.delete(cropped_image_x1, slice(0,x_min),1)

cropped_image_y1 = sp.delete(cropped_image_x2, slice(y_max, 4612), 0)
cropped_image = sp.delete(cropped_image_y1, slice(0,y_min), 0)

masked_image_x1 = sp.delete(mask_image_ones, slice(x_max,2571),1)
masked_image_x2 = sp.delete(masked_image_x1, slice(0,x_min),1)

masked_image_y1 = sp.delete(masked_image_x2, slice(y_max, 4612), 0)
mask_image_ones_cropped = sp.delete(masked_image_y1, slice(0,y_min), 0)
