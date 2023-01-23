# Askja-Radar-Loss
Code used in Shoemaker et al. to determine the loss rate of averages of radargram traces collected using a Geophysical Survey System Inc. ground-penetrating radar.

# Dependencies
Our attenuation code is dependent upon two other Python software packages that must be installed prior to running the script.
1. Readgssi - This sofrware is necessary in order to open and read in Geophysical Survey Systems Inc. (GSSI) ground-penetrating radar files. This software package is available for download at: https://github.com/iannesbitt/readgssi.
2. RAGU - The Radar Analysis Graphical Utility (RAGU) was used to pick reflectors in radargrams collected by a GSSI system. Radargram picks were exported as csv files using this software. These csv files are later read in by the attenuation script. RAGU is available for download at: https://github.com/btobers/RAGU. 

# User Inputs
#### Prior to running the script Users must input the following:
'''
in_dir =
'''

# Outputs
