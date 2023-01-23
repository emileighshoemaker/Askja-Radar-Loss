# Askja-Radar-Loss
Code used in Shoemaker et al. to determine the loss rate of averages of radargram traces collected using a Geophysical Survey System Inc. ground-penetrating radar.

# Dependencies
Our attenuation code is dependent upon two other Python software packages that must be installed prior to running the script.
1. Readgssi - This sofrware is necessary in order to open and read in Geophysical Survey Systems Inc. (GSSI) ground-penetrating radar files. This software package is available for download at: https://github.com/iannesbitt/readgssi.
2. RAGU - The Radar Analysis Graphical Utility (RAGU) was used to pick reflectors in radargrams collected by a GSSI system. Radargram picks were exported as csv files using this software. These csv files are later read in by the attenuation script. RAGU is available for download at: https://github.com/btobers/RAGU. 

# User Inputs
Prior to running the script, Users must input the following at the top of the script:

```
# Data directories for the input data and where to save plots
in_dir =
out_dir = 

# Input the range of traces you wish to average over
slice_left = 
slice_right = 

# Time zero sample number
tz_sample = 

# Depth indicies (samples) to conduct the linear least-squares fit over
top = 
bottom = 

# GPR Center Frequency in Hz
freq = 

# Layer permittivity values are later used in a three-layer depth correction model
# Layer 1 real dielectric permittivity
eps_1 = 

# Layer 2 real dielectric permittivity
eps_2 = 

# Layer 3 real dielectric permittivity
eps_3 = 
```

# Outputs
The script outputs a series of quantities and plots the average trace, corrected trace, and linear least-squares fit to the range of depths specified by inputs from the User.
