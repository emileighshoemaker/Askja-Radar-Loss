#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:39:16 2023

@author: Emileigh
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from readgssi import dzt
import os,sys
import scipy
from scipy import signal
from scipy import stats
from scipy.optimize import curve_fit
from textwrap import wrap
from tkinter import *
import warnings


# Data directories and locations, User Inputs
##################  DIRECTORIES AND USER INPUTS  ##################
in_dir = ""
out_dir = ""

# input the range of traces you wish to average the over
slice_left = 0
slice_right = 1400

# Time Zero Sample number
tz_sample = 91
# Depth indices (samples) to conduct the linear least-squares fit over
top = 276
bottom = 509

########## CONSTANTS ###########
freq = 200e6 # GPR Center Frequency in Hz
# Layer permittivity values are later used in a three-layer depth correction model
eps_1 = 4.37 # layer 1 real dielectric permittivity
eps_2 = 3.20 # layer 2 real dielectric permittivity
eps_3 = 9.00 # layer 3 real dielectric permittivity
eps = [eps_1, eps_2, eps_3]
mean_eps_r = np.mean(eps)  # real dielectric permittivity
c = 2.99e8 # speed of light in m/s for converting twtt to depth
w_freq = 2*np.pi*freq #angular frequency
#wavespeed = c/np.sqrt(eps_r)
lambda_m = c/(freq * np.sqrt(mean_eps_r)) #calculating wavelength in the media given mean permittivity of all layers in 3 layer model
lambda_free = (c/freq)
###############################

# Input the file to analyze
files = glob.glob(in_dir + "/ICE192.prj/csv/" + "*031_P_I11_proc_amp.csv") # input the radargram number
print("FILE EXISTS:", os.path.exists(files[0])) # Check that this file exists

# Read in CSV files created by readgssi within bounds of traces selected by the user
# readgssi is available here: https://github.com/iannesbitt/readgssi
for file in files:
    amp_in = np.loadtxt(file, delimiter= ",", dtype = None) #dtype = np.int32).astype(np.float)

    amp_slice = amp_in[:, slice_left:slice_right]
    amp_slice[amp_slice == 0] = np.nan
    
    # Normalize to the surface reflection or noise by extracting those picks in RAGU
    # Also read in RAGU pick CSV
    # Radar Analysis Graphical Utility available here: https://github.com/btobers/RAGU
    fname = file.split("/")[-1]
    print(fname) #check the file name
    pfix = file.rstrip("/Iceland2019_GPR/ICE192.prj/csv/" + fname) #remove csv directory
    pk_file = pfix + "/Iceland2019_2021_Picks/" + fname.rstrip("_P_I11_proc_amp.csv") + "1_pk_es.csv" #change to directory where normalization picks are stored, add proper extension
    
    # Read in amplitude values to normalize with surf_amp will either be surface amplitude array or noise amplitude array from picks in NOSEpick
    norm_amp = np.loadtxt(pk_file, delimiter = ",", skiprows = 1, usecols = 7, dtype = None) # skip the header and make sure to select proper column
    norm_amp_slice = norm_amp[slice_left:slice_right]
    #norm_amp_slice[np.where(norm_amp_slice == 0)] = np.nan
    #print (surf_amp)
    print ('read successful')
    
    # Read in twtt values from pick csv of layers for use later to calculate accurate reflector depths using calculated permittivity
    pcsv = pd.read_csv(pk_file)
    pcols = pcsv.columns
    tph_base_twtt = pcsv['TephraBase_IceTop__twtt']
    ice_base_twtt = pcsv['IceBase__twtt']
    # take average twtt of those reflectors
    avg_tph_twtt = np.mean(tph_base_twtt[slice_left:slice_right])
    avg_ice_base_twtt = np.mean(ice_base_twtt[slice_left:slice_right])
    
    # find the value in the time zeroed twtt array that is closest to the average value for a particular reflector
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]
    
    
    # divide by the surface reflection to calculate relative amplitudes
    # convert to dB after, divide amplitudes first    
    div_amps = np.zeros(np.shape(amp_slice))
    div_amps_gain = np.zeros(np.shape(amp_slice))
    for i in range(amp_slice.shape[1]):
        div_amps[:,i] = amp_slice[:,i] / norm_amp_slice[i]
        
    # Convert to dB, plot horizontally averaged power with depth       
    #avg_amps = np.nanmean(div_amps**2, axis = 1)
    pow_dB = 20*np.log10((div_amps))
    horizontal_pow_avg = (np.nanmean(pow_dB, axis = 1))
    horizontal_pow_avg = np.nan_to_num(horizontal_pow_avg, nan = 0)
    

    # Extract two way travel time from radargram via readgssi (make sure you only use time-zeroed processed rgram - no other processing performed)
    header, data, gps = dzt.readdzt(pfix + "/Iceland2019_GPR/ICE192.prj/" + fname.rstrip("_P_I11_proc_amp.csv") + "1.DZT")
    twtt_full = np.arange(900)* header['ns_per_zsample']   # array of two way travel times extracted from readgssi for the radargram
    twtt_surf = twtt_full[tz_sample:900] - twtt_full[tz_sample]
    #depth = ((c*twtt_surf)/(2 * np.sqrt(eps_r)))
    
    # find the index in the time zeroed array corresponding to the average value of the tephra base reflector
    v1 = find_nearest(twtt_surf, avg_tph_twtt)
    t1, = np.where(twtt_surf == v1)
    # new array of depths from surface to tephra base
    z1 = (c * twtt_surf[0: t1[0]]) / (2 * np.sqrt(eps_1))
    
    # find the index in the time zeroed array corresponding to the average value of the ice base reflector
    v2 = find_nearest(twtt_surf, avg_ice_base_twtt)
    t2, = np.where(twtt_surf == v2)
    # new array of depths from the tephra base to the ice base
    z2 = ((c * twtt_surf[0 : t2[0] - t1[0]]) / (2 * np.sqrt(eps_2)) ) + np.max(z1)
    
    # new array of depths from the ice base to the end of the twtt array
    t3, = np.where(twtt_surf == np.max(twtt_surf))
    z3 = ((c * twtt_surf[0: t3[0] + 1 - t2[0]]) / (2 * np.sqrt(eps_3))) + np.max(z2) 
    
    z = np.concatenate((z1, z2, z3), axis = 0) # add the three new depth corrected sections together
    y = z
    x = horizontal_pow_avg[tz_sample:900]
    
    # Geometric Spreading Correction, loss tangent/Q factor, and Plotting 
##########################################################################################
    #R = x # set R equal depth 
    geom_spread2 = y**-2

    geom_spread2_dB = 20*np.log10(geom_spread2)
    
    corrected_loss = x - geom_spread2_dB

    slope, intercept, r_value, p_value, std_err = stats.linregress(y[top:bottom], corrected_loss[top:bottom])
    line = slope*y[top:bottom] + intercept
    print("slope = ", slope, "p-value", p_value, "r-value", r_value, "r-squared", r_value**2)

    tand = (np.abs(slope/8.686) * lambda_m)/np.pi #relationship between loss tangent and inverse of Q factor (have to take out of dB to isolate the spatial attenuation coefficient alpha)
    print("tand =", tand)
    
    Q = 1/(tand)
    print("Q =", Q)
    

 ########################################  ERROR CALCULATION   #################################################
    
    
    slopelist = [] # create list to store slopes for future calculations
    interceptlist = [] # create a list for the intercepts as well
    for i in range(slice_left, slice_right, 10):
    
        amp_slice2 = amp_in[:, i-5:i+5] # add five on either side to make ten because Python
        amp_slice2[amp_slice2 == 0] = np.nan
        
        # Normalize to the surface reflection or noise by extracting those picks in RAGU
        # Also read in RAGU pick CSV
        # Radar Analysis Graphical Utility available here: https://github.com/btobers/RAGU
        fname = file.split("/")[-1]
        #print(fname) #check the file name
        pfix = file.rstrip("/Iceland2019_GPR/ICE192.prj/csv/" + fname) #remove csv directory
        pk_file = pfix + "/Iceland2019_2021_Picks/" + fname.rstrip("_P_I11_proc_amp.csv") + "1_pk_es.csv" #change to directory where normalization picks are stored, add proper extension
        
        # Read in amplitude values to normalize with surf_amp will either be surface amplitude array or noise amplitude array from picks in NOSEpick
        norm_amp2 = np.loadtxt(pk_file, delimiter = ",", skiprows = 1, usecols = 7, dtype = None) # skip the header and make sure to select proper column
        norm_amp_slice2 = norm_amp2[i-5:i+5]
        
        # divide by the surface reflection or noise floor in order to get powers relative to a known value
        # convert to dB after, divide amplitudes first    
        div_amps2 = np.zeros(np.shape(amp_slice2))
        for i in range(amp_slice2.shape[1]):
            div_amps2[:,i] = amp_slice2[:,i] / norm_amp_slice2[i]

            
        # Convert to dB, plot horizontally averaged power with depth       
        pow_dB2 = 10*np.log10((div_amps2**2))
        horizontal_pow_avg2 = (np.nanmean(pow_dB2, axis = 1))
        horizontal_pow_avg2 = np.nan_to_num(horizontal_pow_avg2, nan = 0)
        
    
        # Extract two way travel time from radargram via readgssi (make sure you only use time-zeroed processed rgram - no other processing performed)
        header, data, gps = dzt.readdzt(pfix + "op/Iceland2019_GPR/ICE192.prj/" + fname.rstrip("_P_I11_proc_amp.csv") + "1.DZT")
        twtt_full2 = np.arange(900)* header['ns_per_zsample']   # array of two way travel times extracted from readgssi for the radargram
        twtt_surf2 = twtt_full2[tz_sample:900] - twtt_full2[tz_sample]
        
        ############## 3 layer depth correction ##################################################################
        # find the index in the time zeroed array corresponding to the average value of the tephra base reflector
        v12 = find_nearest(twtt_surf2, avg_tph_twtt)
        t12, = np.where(twtt_surf2 == v12)
        # new array of depths from surface to tephra base
        z12 = (c * twtt_surf2[0: t12[0]]) / (2 * np.sqrt(eps_1))
        
        # find the index in the time zeroed array corresponding to the average value of the ice base reflector
        v22 = find_nearest(twtt_surf2, avg_ice_base_twtt)
        t22, = np.where(twtt_surf2 == v22)
        # new array of depths from the tephra base to the ice base
        z22 = ((c * twtt_surf2[0: t22[0] - t12[0]]) / (2 * np.sqrt(eps_2))) + np.max(z12)
        
        # new array of depths from the ice base to the end of the twtt array
        t32, = np.where(twtt_surf2 == np.max(twtt_surf2))
        z32 = ((c * twtt_surf2[0: t32[0] + 1 - t22[0]]) / (2 * np.sqrt(eps_3))) + np.max(z22) #assuming permittivity of 9 messes things up right now
        
        z2 = np.concatenate((z12, z22, z32), axis = 0) # add the three new depth corrected sections together
        
        y2 = z2
        x2 = horizontal_pow_avg2[tz_sample:900]
        
        # Geometric Spreading Correction, loss tangent/Q factor
    ##########################################################################################
        #R = x # set R equal depth 
        geom_spread22 = y2**-2
    
        geom_spread22_dB = 20*np.log10(geom_spread22)
        
        corrected_loss2 = x2 - geom_spread22_dB
    
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(y2[top:bottom], corrected_loss2[top:bottom])
        line2 = slope2 * (y2[top:bottom]) + intercept2
        # Taking reciprocal of slope to get dB/m
        #print("slope = ", slope2**-1, "intercept =", intercept2, "p-value", p_value2, "r-value", r_value2, "r-squared", r_value2**2)
        slopelist.append(slope2) # append slopes from each loop
        interceptlist.append(intercept2)
        minslope = np.min(slopelist)
        maxslope = np.max(slopelist)
        indmin = slopelist.index(np.min(slopelist)) # grab the index corresponding to the min slope
        indmax = slopelist.index(np.max(slopelist)) # grab the index corresponding to the max slope
        minline = minslope * (y2[top:bottom]) + interceptlist[indmin] # define the line associated with the min slope
        maxline = maxslope * (y2[top:bottom]) + interceptlist[indmax] # define the line associated with the max slope
        
    
    
    slopevar = np.var(slopelist) # calculate the variance in calculated slopes
    print("variance in slope =", slopevar)
    
    
##################### END ERROR CALCULATION SECTION ##########################
    
#####################   PLOTTING  ############################################ 
    fig, ax = plt.subplots()
    ax.plot(x, y, c = "k", label = "Average Trace")
    ax.plot(corrected_loss, y, 'grey', label = "Corrected Losses")
    ax.plot(line, y[top:bottom], "b--", label = "Total Loss Fit")
    plt.suptitle("2019 Radargram 31, 200 MHz , Caldera Site 3") 
            
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = True
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
    ax.set_xlabel("Signal (dB)")
    plt.ylabel ("Depth (m)")
    
    plt.gca().invert_yaxis()
    plt.xlim(xmin = -60, xmax = 0)
    plt.ylim(ymin = 6.0, ymax = 0)
    plt.legend(loc = "lower right")
    plt.show()
    #plt.savefig(out_dir + "Site1_200mhz" + fname.rstrip(".csv") + ".jpg", dpi = 300) 
    
    warnings.simplefilter(action = 'ignore', category = FutureWarning)
    warnings.simplefilter(action = 'ignore', category = RuntimeWarning)
    
    
    