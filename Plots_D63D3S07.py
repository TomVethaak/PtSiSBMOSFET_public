# -*- coding: utf-8 -*-
from Standardized_plots import *
from DataDirectoryLibrary import *

# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

DataDirectory   = DataDirectoryD63D3S07
IMG_base        = 'img/D63D3S07'
SweepType       = "TVd"    # VgVd, HVd, TVd
# The X axis can be automatically extracted from VgVd file names
# For H or T sweeps all we get in the file names is #0001, #0002 etc
SweepStart      = 1.05      # The H or T value of the first file (#0001)
SweepStep       = -0.05     # The H or T step in Labview units

# Global variables are carried by the dict Global
PhaseShift      = 0       # degrees
Rfilters        = 41240/2     # Ohm
Cfilters        = 11.3e-9   # Farad
freq            = 140
Global          = {
    # Interpolation overhead: how many extra points are we going to sample?
    'InterpolationOverhead'     : 10,
    
    # At higher lock-in frequencies we get a shift in phase.
    'PhaseShift'                : PhaseShift,   # degrees
    'Rfilters'                  : Rfilters,     # Ohm
    'Cfilters'                  : Cfilters,     # Farad
    'freq'                      : freq,         # Hertz
    'Zfilters'                  : Rfilters/(1+1j*Rfilters*Cfilters*2*np.pi*freq),     # Complex number
    'alpha'                     : 0.88,
    'VdOffset'                  : -0.03,
    
    # Extracted gap in mV
    'DeltaPtSi'                 : 0.31
    }

# Color map
RootDirectory   = GetRootDirectory( User )
rcParams['image.cmap']='viridis'

#%% IMPORT
# This cell is for importing the data, it takes the longest.

if User == "Francois":
    # Francois stores the data on his C drive
    RootDirectory   = "C:/Users/FL134692/Documents/ANR - Projets/2019 SUNISIDEUP/Data/C2N_PtSi/"
elif User == "Tom":
    # Tom stores the data on his Z drive
    RootDirectory   = "C:/Users/vethaak/Documents/PhD/Data/Data sorted by lot/W12 Laurie Calvet's wafer/Christophe's fridge/"
else:
    print("Error: User not found")
    exit

# Start by listing the files in the directory
DataFileList    = np.array( [ entry for entry in os.listdir( RootDirectory + DataDirectory )[1*bool(SweepType=="TVd" or SweepType=="HVd")::] ] )
DataArray       = np.array( [ np.loadtxt( RootDirectory + DataDirectory + Filename ) for Filename in DataFileList ])

# We are dealing with two naming systems: either Vg=-X.12345V_Y.dat, or Vg=-X.1234V_.dat.
# The most robust way of extracting the Vg value seems to find any expression between "Vg=" and "V"
# Let's use regular expressions: 
# 'DataFileList' is an array with strings.
# For each element 'filename' we use the regular expression "(?<=Vg=).*(?=V)" to find all numbers that are sandwiched between "Vg=" and "V"
# We then take the first (and only) element of that list, index [0], and convert it from a string to a float (a number that uses two bytes)
if SweepType == "VgVd":
    XValues     = np.array( [ float( re.findall( "(?<=Vg=).*(?=V)" , Filename )[0] ) for Filename in DataFileList ] )
elif SweepType == "TVd" or "HVd":
    XValues     = np.array( [ SweepStart + XIterator * SweepStep for XIterator in range( 0 , len( DataArray ) ) ] )
else:
    print("Error: SweepType not found")
    exit

# If the measurement started at larger Vg values, and then moved to smaller (or more negative) values, invert the VgValues and DataArray arrays
# imshow requires the arrays to be sorted from left to right.
if XValues[0]>XValues[-1]:
    XValues     = XValues[::-1]
    DataArray   = DataArray[::-1]

#%% Select data
# Do not interpolate over all data to prevent memory issues
DataXRange  = np.array([-10,0])     # V     
DataYRange  = np.array([-10,10])    # mV

#%% Plot a few selected H sweep curves in one plot

XCurveSelection = [0.05,    0.6,    0.8,    0.9,    1.0]
VdOffsetPerX    = [0.089,   0.092,  0.091,  0.089,  0.089]
PlotXRange      = np.array([-0.11,0.11])
PlotYRange      = np.array([0,225])

Plot_Gdiff_T_manual( Global, IMG_base+'_Gdiff_Vg-4.0V_T_large.png', DataArray, XValues, XCurveSelection, PlotXRange, PlotYRange, VdOffsetPerX )