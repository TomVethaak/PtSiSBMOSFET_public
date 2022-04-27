# -*- coding: utf-8 -*-
from Standardized_plots import *
from DataDirectoryLibrary import *

# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

# Global variables are carried by the dict Global
PhaseShift      = 0       # degrees
Rfilters        = 41240/2     # Ohm
Cfilters        = 11.3e-9   # Farad
freq            = 140
Global          = {
    # Interpolation overhead: how many extra points are we going to sample?
    'InterpolationOverhead'     : 10,
    
    # At higher lock-in frequencies we get a shift in phase.
    'PhaseShift'                : PhaseShift,       # degrees
    'Rfilters'                  : Rfilters,     # Ohm
    'Cfilters'                  : Cfilters,   # Farad
    'freq'                      : freq,        # Hertz
    'Zfilters'                  : Rfilters/(1+1j*Rfilters*Cfilters*2*np.pi*freq),     # Complex number
    'alpha'                     : 0.88,
    'VdOffset'                  : -0.03,
    
    # Extracted gap in mV
    'DeltaPtSi'                 : 0.31
    }

#%% Initialize 
RootDirectory   = GetRootDirectory( User )      
rcParams['image.cmap']='viridis'                # Color map

#%% -------- %%#
#   HVd scan   #
#   --------   #

DataDirectory   = DataDirectoryD63D3S08
IMG_base        = 'img/D63D3S07_08'
SweepType       = "HVd"    # VgVd, HVd, TVd
# The X axis can be automatically extracted from VgVd file names
# For H or T sweeps all we get in the file names is #0001, #0002 etc
SweepStart      = 0         # The H or T value of the first file (#0001)
SweepStep       = 0.005     # The H or T step in Labview units

#%% IMPORT

# Start by listing the files in the directory
HDataFileList   = GetDataFileList( RootDirectory + DataDirectory, SweepType )
HDataArray      = GetDataArray( HDataFileList, RootDirectory + DataDirectory )
HValues         = GetXValues( SweepType, SweepStart, SweepStep, HDataFileList, HDataArray )

# If the measurement started at larger Vg values, and then moved to smaller (or more negative) values, invert the VgValues and DataArray arrays
# imshow requires the arrays to be sorted from left to right.
if HValues[0]>HValues[-1]:
    HValues     = HValues[::-1]
    HDataArray  = HDataArray[::-1]

#%% Plot a few selected H sweep curves in one plot

HSelection      = [0,    0.01,    0.02,    0.03,   0.04]
VdOffsetPerH    = [0.095,  0.095,  0.09,  0.09,  0.09]
PlotXRange      = np.array([-0.11,0.11])
HPlotYRange      = np.array([0,200])

#%% -------- %%#
#   TVd scan   #
#   --------   #

DataDirectory   = DataDirectoryD63D3S07
SweepType       = "TVd"    # VgVd, HVd, TVd
# The X axis can be automatically extracted from VgVd file names
# For H or T sweeps all we get in the file names is #0001, #0002 etc
SweepStart      = 1.05      # The H or T value of the first file (#0001)
SweepStep       = -0.05     # The H or T step in Labview units

#%% IMPORT

# Start by listing the files in the directory
TDataFileList   = GetDataFileList( RootDirectory + DataDirectory, SweepType )
TDataArray      = GetDataArray( TDataFileList, RootDirectory + DataDirectory )
TValues         = GetXValues( SweepType, SweepStart, SweepStep, TDataFileList, TDataArray )

# If the measurement started at larger Vg values, and then moved to smaller (or more negative) values, invert the VgValues and DataArray arrays
# imshow requires the arrays to be sorted from left to right.
if TValues[0]>TValues[-1]:
    TValues     = TValues[::-1]
    TDataArray  = TDataArray[::-1]

#%% Plot a few selected H sweep curves in one plot

TSelection      = [0.05,    0.6,    0.8,    0.9,    1.0]
VdOffsetPerT    = [0.089,   0.092,  0.091,  0.089,  0.089]
TPlotYRange     = np.array([0,225])

#%% ComboPlot
ComboPlot_Gdiff_H_T_manual( Global, IMG_base+'_Gdiff_Vg-4.0V_T_Vg-4.5V_H_large.png', HDataArray, HValues, HSelection, PlotXRange, HPlotYRange, VdOffsetPerH,
                                                                                     TDataArray, TValues, TSelection,             TPlotYRange, VdOffsetPerT)