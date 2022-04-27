# -*- coding: utf-8 -*-
from Standardized_plots import *
from DataDirectoryLibrary import *

# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

DataDirectory   = DataDirectoryD63D3S08
IMG_base        = 'img/D63D3S08'
SweepType       = "HVd"    # VgVd, HVd, TVd
# The X axis can be automatically extracted from VgVd file names
# For H or T sweeps all we get in the file names is #0001, #0002 etc
SweepStart      = 0         # The H or T value of the first file (#0001)
SweepStep       = 0.005     # The H or T step in Labview units

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

#%% Initialize 
RootDirectory   = GetRootDirectory( User )
rcParams['image.cmap']='viridis'                # Color map

#%% IMPORT

DataFileList    = GetDataFileList( RootDirectory + DataDirectory, SweepType )
DataArray       = GetDataArray( DataFileList, RootDirectory + DataDirectory )
XValues         = GetXValues( SweepType, SweepStart, SweepStep, DataFileList, DataArray )

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

XCurveSelection = [0,    0.01,    0.02,    0.03,   0.04]
VdOffsetPerX    = [0.095,  0.095,  0.09,  0.09,  0.09]
PlotXRange      = np.array([-0.11,0.11])
PlotYRange      = np.array([0,200])

Plot_Gdiff_H_Manual( Global, IMG_base+'_Gdiff_Vg-4.5V_H_large.png', DataArray, XValues, XCurveSelection, PlotXRange, PlotYRange, VdOffsetPerX )