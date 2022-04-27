from Standardized_plots import *

# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

# D63D3
DataDirectoryD63D3S04   = "2020-11 PtSi D63/D3/Base Temp/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-1V_-4.5V(0.01)/"
DataDirectoryD63D3S11   = "2020-11 PtSi D63/D3/Base Temp/H=100mT/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-4.5_-1(0.01)V/"
DataDirectoryD63D3S12   = "2020-11 PtSi D63/D3/Tsweep Vg=-4.5V H=100mT/"
DataDirectoryD63D3S13   = "2020-11 PtSi D63/D3/Tsweep Vg=-2.7V H=100mT/"
DataDirectoryD63D3S14   = "2020-11 PtSi D63/D3/Base Temp/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-4.9V_-4.0V(0.01)/"
DataDirectoryD63D3S15   = "2020-11 PtSi D63/D3/Base Temp/D3_Gdiff_Vd=+-5(0.01)mV_Vac=10uV_Vg=-4.5V_-2.93V(0.01)/"
DataDirectoryD63D3S16   = "2020-11 PtSi D63/D3/run2/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-2.7V_-3V(0.01)/"
DataDirectoryD63D3S07   = "2020-11 PtSi D63/D3/Tsweep Vg=-4.0V/"

# D63D4
DataDirectoryD63D4S01   = "2020-11 PtSi D63/D4/Base temp/D4_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-4.0V_-1.0V(0.01)/"
DataDirectoryD63D4S02   = "2020-11 PtSi D63/D4/Base temp/D4_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-2.2V_-1.0V(0.0002)/"
DataDirectoryD63D4S03   = "2020-11 PtSi D63/D4/Base temp/D4_Gdiff_Vd=+-5(0.01)mV_Vac=10uV_Vg=-4.0V_-1.0V(0.01)/"

# D61D4
DataDirectoryD61D4S01   = "2021-02 PtSi D61/D4/Base Temp/D4_Gdiff_Vd=+-5(0.01)mV_Vac=10uV_Vg=-2V_-4.15V(0.01)/"

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

# Color map
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

#%% Get data
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

Plot_Gdiff( Global, IMG_base+'_Gdiff.png', DataArray, VgValues, XRange, YRange )

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