# -*- coding: utf-8 -*-
from Standardized_plots import *
from DataDirectoryLibrary import *

# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

DataDirectory   = DataDirectoryD73D4S02
IMG_base        = 'img/D73D4S02'
SweepType       = "VgVd"    # VgVd, HVd, TVd
# The X axis can be automatically extracted from VgVd file names
# For H or T sweeps all we get in the file names is #0001, #0002 etc
SweepStart      = 1         # The H or T value of the first file (#0001)
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

#%% Initialize 
RootDirectory   = GetRootDirectory( User )

#%% IMPORT
DataFileList    = GetDataFileList( RootDirectory + DataDirectory, SweepType )
DataArray       = GetDataArray( DataFileList, RootDirectory + DataDirectory )
XValues         = GetXValues( SweepType, SweepStart, SweepStep, DataFileList, DataArray )

# If the measurement started at larger Vg values, and then moved to smaller (or more negative) values, invert the VgValues and DataArray arrays
# imshow requires the arrays to be sorted from left to right.
if XValues[0]>XValues[-1]:
    XValues     = XValues[::-1]
    DataArray   = DataArray[::-1]

#%% Get data
# Do not interpolate over all data to prevent memory issues
DataXRange  = np.array([-10,0])     # V     
DataYRange  = np.array([-10,10])    # mV

#%% Combo plot 2D Gdiff, slices

VgValues    = [-4.4, -5.0, -5.6, -5.7]
XRange0     = [np.min(DataArray[:,0,5]), np.max(DataArray[:,0,5])]
YRange0     = [-0.4,0.4]
XRange10    = [-0.4,0.4]
XRange11    = [-0.15,0.15]
YRange10    = [0, 1.4E2]
YRange11    = [0, 0.8E2]
YOffset     = -0.14

# 2D plot
Plot_2D_Gdiff( Global, IMG_base + '_GdiffDUT.png', DataArray, VgValues, DataXRange, DataYRange, XRange0,      YRange0,    YOffset, 'SETUP')
# Combo 2D and 1D plots
ComboPlot_Gdiff( Global, IMG_base + '_Vd_combo_plot2.png', DataArray, VgValues, DataXRange, DataYRange, XRange0, XRange10, XRange11,
                YRange0, YRange10, YRange11, YOffset, 'SETUP')

#%% Plot only Gdiff at 4 Vg values

VgValues    = [-4.4,    -5.0,   -5.6,   -5.7]
XShifts     = [0.028,   0.030,  0.028,  0.035]
XRange      = [-0.39,0.39]
YRange      = [0, 1.3E2]

Plot_Gdiff_Vg_Manual( Global, IMG_base + '_Gdiff.png', DataArray, VgValues, XRange, YRange, XShifts )

#%% Combo plot 2D Gdiff and A, Gamma

XRange      = [ min(DataArray[:,0,5]), max(DataArray[:,0,5]) ]
YRange      = [-0.5,0.5]
YOffset     = -0.14

ComboPlot_Gdiff_Gamma_A( Global, IMG_base + '_Gdiff_interface_combo_plot.png',
                        DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, 'SETUP')

#%% Plot I(V)

VgValues    = [-5.3, -5.5, -5.7]
XRange      = [0,1.2]
YRange      = [0,1.2]

Plot_Isn( Global, IMG_base+'_I.png', DataArray, VgValues, XRange, YRange )