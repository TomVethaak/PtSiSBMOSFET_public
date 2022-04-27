# -*- coding: utf-8 -*-
from Standardized_plots import *
from DataDirectoryLibrary import *

# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

IMG_base        = 'img/D63D3S04_S11'
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

#%% D63D3S04 -----------------------------------------------------------------
DataDirectory   = DataDirectoryD63D3S04

#%% IMPORT
S04DataFileList    = GetDataFileList( RootDirectory + DataDirectory, SweepType )
S04DataArray       = GetDataArray( S04DataFileList, RootDirectory + DataDirectory )
S04XValues         = GetXValues( SweepType, SweepStart, SweepStep, S04DataFileList, S04DataArray )

# If the measurement started at larger Vg values, and then moved to smaller (or more negative) values, invert the VgValues and DataArray arrays
# imshow requires the arrays to be sorted from left to right.
if S04XValues[0] > S04XValues[-1]:
    S04XValues     = S04XValues[::-1]
    S04DataArray   = S04DataArray[::-1]

#%% D63D4S11 -----------------------------------------------------------------
DataDirectory   = DataDirectoryD63D3S11

#%% IMPORT
S11DataFileList    = GetDataFileList( RootDirectory + DataDirectory, SweepType )
S11DataArray       = GetDataArray( S11DataFileList, RootDirectory + DataDirectory )
S11XValues         = GetXValues( SweepType, SweepStart, SweepStep, S11DataFileList, S11DataArray )

# If the measurement started at larger Vg values, and then moved to smaller (or more negative) values, invert the VgValues and DataArray arrays
# imshow requires the arrays to be sorted from left to right.
if S11XValues[0] > S11XValues[-1]:
    S11XValues     = S11XValues[::-1]
    S11DataArray   = S11DataArray[::-1]

#%% Plot

S04DataFilter       = remove_filters( Global, S04DataArray )
S04DataSelection    = Select_Data_XRange( S04DataFilter, [ -3.0, -2.7 ] )
S04DataAverage      = S04DataSelection.mean( axis = 0 )
S04CenteredAverage  = Center_Curve( savgol_filter( ( S04DataAverage[:,0] , S04DataAverage[:,2] ) , 3, 1) )

S11DataFilter       = remove_filters( Global, S11DataArray )
S11DataSelection    = Select_Data_XRange( S11DataFilter, [ -3.0, -2.7 ] )
S11DataAverage      = S11DataSelection.mean( axis = 0 )
S11CenteredAverage  = Center_Curve( savgol_filter( ( S11DataAverage[:,0] , S11DataAverage[:,2] ) , 3, 1) )

XShifts             = [ 0.103, 0.045 ]
XRange              = [ -0.39, 0.39 ]
YRange              = [ 0, 5.6E1 ]
ListPlot_Gdiff_Vg_Manual( Global, IMG_base + '_combo_line_plot.png', [ S04CenteredAverage, S11CenteredAverage ], [ 'H = 0', 'H = 100 mT' ], XRange, YRange, XShifts )

# fig     = plt.figure(figsize=(5,3.7));
# plt.xlabel( '$eV_\mathrm{d}/\Delta$',       fontsize=11);
# plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
# plt.axvline(x = 0,  color = 'black',    linestyle = '--');
# plt.axvline(x = -Global['DeltaPtSi'], color = 'grey',     linestyle = '--');
# plt.axvline(x = Global['DeltaPtSi'],  color = 'grey',     linestyle = '--');
# plt.xlim( [ -0.39, 0.39 ] )
# plt.ylim( [ 0, 5.6E1 ] )
# plt.plot( S04CenteredAverage[0]+0.104, 1E6*S04CenteredAverage[1], label='H = 0')
# plt.plot( S11CenteredAverage[0]+0.045, 1E6*S11CenteredAverage[1], label='H = 100 mT' )
# handles, labels = plt.gca().get_legend_handles_labels()
# plt.legend( handles, labels, loc='upper center' )
# plt.savefig( IMGfilename, dpi=GlobalDPI, bbox_inches='tight' )
# plt.show()
plt.close()