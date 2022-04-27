# -*- coding: utf-8 -*-
from Standardized_plots import *
from DataDirectoryLibrary import *

# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

DataDirectory   = DataDirectoryD61D4S01
IMG_base        = 'img/D61D4S01'
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

# Color map
RootDirectory   = GetRootDirectory( User )
rcParams['image.cmap']='viridis'

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

#%% Combo plot 2D Gdiff, slices

VgValues    = [-2.99,  -3.60, -4.14]
DataXRange  = [-10,10]
DataYRange  = [-10,10]
XRange0     = [np.min(DataArray[:,0,5]), np.max(DataArray[:,0,5])]
YRange0     = [-0.5,0.5]
XRange10    = [-0.5,0.5]
XRange11    = [-0.15,0.15]
YRange10    = [0, 2.7E2]
YRange11    = [0, 2.5E1]
YOffset     = -0.04

# 2D plot
Plot_2D_Gdiff( Global, IMG_base+'_GdiffDUT.png', DataArray, VgValues, DataXRange, DataYRange, XRange0, YRange0, YOffset, 'SETUP')
# Combo 2D and 1D plots
ComboPlot_Gdiff( Global, IMG_base+'_Vd_combo_plot.png', DataArray, VgValues, DataXRange, DataYRange, XRange0, XRange10, XRange11,
                YRange0, YRange10, YRange11, YOffset, 'SETUP')

#%% Combo plot 2D Gdiff and A, Gamma

XRange      = [ min(DataArray[:,0,5]), max(DataArray[:,0,5]) ]
YRange      = [-0.5,0.5]
YOffset     = -0.04

ComboPlot_Gdiff_Gamma_A( Global, IMG_base+'_Gdiff_interface_combo_plot.png',
                        DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, 'SETUP')

#%% Combo plot 2D Gdiff and A, Gamma (SMALL)

XRange      = [ min(DataArray[:,0,5]), max(DataArray[:,0,5])-0.1 ]
YRange      = [-0.5,0.5]
YOffset     = -0.04

ComboPlot_Gdiff_Gamma_A_small( Global, IMG_base+'_Gdiff_interface_combo_plot_small.png',
                        DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, 'SETUP')

#%% Extract "RN" from minimum

VgValues    = [-2.99,  -3.60, -4.14]
XRange      = [-0.5,0.5]
YRange      = [0,6E1]

Plot_Gdiff_Vg( Global, IMG_base+'_Gdiff.png', DataArray, VgValues, XRange, YRange )

#%% Plot I(V)

VgValues    = [-2.99,  -3.60, -4.14]
XRange      = [0,1.2]
YRange      = [0,1.2]

Plot_Isn( Global, IMG_base+'_I.png', DataArray, VgValues, XRange, YRange )

#%% Plot Rdiff and RS

XRange      = [ min(DataArray[:,0,5]), max(DataArray[:,0,5] ) ]
YRange      = [ 0, 1E6 ]
YRangeLog   = [ 1E4, 2E8 ]
#Plot_RNdiffRN( IMG_base+'_Rdiff_RN_Vg.png', DataArray, XRange, YRange )
LogPlot_RNdiffRN( Global, IMG_base+'_Rdiff_RS_Vg_log.png', DataArray, XRange, YRangeLog )

#%% Plot ratio RS/Rdiff

XRange      = [ min(DataArray[:,0,5]), max(DataArray[:,0,5] ) ]
YRange      = [ 0, 10 ]
YRangeLog   = [ 1E4, 2E8 ]

Plot_RNdiff_RN_Ratio( Global, IMG_base+'_Rdiff_RS_Vg_ratio.png', DataArray, XRange, YRange )