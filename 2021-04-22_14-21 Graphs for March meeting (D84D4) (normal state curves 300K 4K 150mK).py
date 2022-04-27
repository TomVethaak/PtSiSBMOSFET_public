from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
rcParams.update({'figure.autolayout': True})
# Select the user that is running the code
User            = "Tom"
#User           = "Francois"

# D63D3
DataDirectoryD63D3S04   = "2020-11 PtSi D63/D3/Base Temp/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-1V_-4.5V(0.01)/"
DataDirectoryD63D3S11   = "2020-11 PtSi D63/D3/Base Temp/H=100mT/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-4.5_-1(0.01)V/"
DataDirectoryD63D3S12   = "2020-11 PtSi D63/D3/Tsweep Vg=-4.5V H=100mT/"
DataDirectoryD63D3S13   = "2020-11 PtSi D63/D3/Tsweep Vg=-2.7V H=100mT/"
DataDirectoryD63D3S14   = "2020-11 PtSi D63/D3/Base Temp/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-4.9V_-4.0V(0.01)/"
DataDirectoryD63D3S15   = "2020-11 PtSi D63/D3/Base Temp/D3_Gdiff_Vd=+-5(0.01)mV_Vac=10uV_Vg=-4.5V_-2.06V(0.01)/"
DataDirectoryD63D3S16   = "2020-11 PtSi D63/D3/run2/D3_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-2.7V_-3V(0.01)/"
DataDirectoryD63D3S07   = "2020-11 PtSi D63/D3/Tsweep Vg=-4.0V/"


# D63D4
DataDirectoryD63D4S01   = "2020-11 PtSi D63/D4/Base temp/D4_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-4.0V_-1.0V(0.01)/"
    # Diamonds:
DataDirectoryD63D4S02   = "2020-11 PtSi D63/D4/Base temp/D4_Gdiff_Vd=+-1(0.01)mV_Vac=10uV_Vg=-2.2V_-1.0V(0.0002)/"
DataDirectoryD63D4S03   = "2020-11 PtSi D63/D4/Base temp/D4_Gdiff_Vd=+-5(0.01)mV_Vac=10uV_Vg=-4.0V_-1.0V(0.01)/"

# D61D4
    # Wide scan showing both the coherence peaks and the ZBCP:
DataDirectoryD61D4S01   = "2021-02 PtSi D61/D4/Base Temp/D4_Gdiff_Vd=+-5(0.01)mV_Vac=10uV_Vg=-2V_-4.15V(0.01)/"
    # 80K
DataFileD61D4C01        = "2021-02 PtSi D61/D4_Gdiff_Vg=0_-2(0.01)V_Vac=10uV_Vd=+1mV_T80K.dat"
    # 80K bis
DataFileD61D4C02        = "2021-02 PtSi D61/D4_Gdiff_Vg=0_-2(0.01)V_Vac=10uV_Vd=+1mV_T80Kbis.dat"
    # 4K
DataFileD61D4C03        = "2021-02 PtSi D61/D4_Gdiff_Vg=0_-2(0.01)V_Vac=10uV_Vd=+1mV_T4K.dat"
    # 45mK
DataFileD61D4C04        = "2021-02 PtSi D61/D4_Gdiff_Vg=0_-2(0.01)V_Vac=10uV_Vd=+1mV_T45mK.dat"

# D84D4 from Bigoudene
DataFileD84D4C01        = "D84/300K/D4_300K_Vd=1mV_VSub=f_Vg=pm1.0(0.02)V_run2.dat"
DataFileD84D4C02        = "D84/4K/D4_Vd=1mV_VSub=f_Vg=pm3.5(0.05)V_run2.dat"
DataFileD84D4C03        = "D84/150mK/D4_Vd=0.12mV_VSub=f_Vg=0_-4(0.01)V_run1.dat"

DataFile        = DataFileD61D4C01
DataFileList    = np.array([DataFileD84D4C01,DataFileD84D4C02,DataFileD84D4C03])
DataDirectory   = DataDirectoryD63D4S02
SweepType       = "VgVd"    # VgVd, HVd, TVd
# The X axis can be automatically extracted from VgVd file names
# For H or T sweeps all we get in the file names is #0001, #0002 etc
SweepStart      = 1         # The H or T value of the first file (#0001)
SweepStep       = -0.05     # The H or T step in Labview units

# Interpolation overhead: how many extra points are we going to sample?
InterpolationOverhead   = 10

# At higher lock-in frequencies we get a shift in phase.
PhaseShift      = 0       # degrees
Rfilters        = 41240/2     # Ohm
Cfilters        = 11.3e-9   # Farad
freq            = 140        # Hertz
Zfilters        = Rfilters/(1+1j*Rfilters*Cfilters*2*np.pi*freq)     # Complex number
alpha           = 0.88
VdOffset        = -0.03

#%% INITIALIZE
# This cell is for initialization: import the relevant libraries/modules, choose what directory we are importing things from

import bisect
import guidata                        # guidata
from guiqwt.plot import ImageDialog   # ImageDialog does ???
from guiqwt.builder import make       # make does ???
from guidata.qt.QtGui import QFont
import numpy as np                    # numpy gives us numerical analysis functions
import os                             # os, "miscellaneous operating system interfaces", is for handling directory/file operations. Useful for getting the files in a directory
import re                             # re for regular expressions (https://docs.python.org/3/howto/regex.html)
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

_app        = guidata.qapplication()

def remove_filters1D( _DataArray ):

    # Create an empty array of the same size as the original data array,
    # - for the adjusted values:
    AdjustedDataArray       = np.zeros( ( len( _DataArray ) ,
                                      len( _DataArray[0] ) ) )
    for VdIterator in range( 0 , len( _DataArray ) ):
        # Values given in the data files
        _Vd             = _DataArray[ VdIterator , 0 ]
        _Id             = _DataArray[ VdIterator , 1 ]
        _GdiffAmplitude = _DataArray[ VdIterator , 2 ] 
        _GdiffPhase     = _DataArray[ VdIterator , 3 ]
        _ILeak          = _DataArray[ VdIterator , 4 ]
        _Vg             = _DataArray[ VdIterator , 5 ]
        _T              = _DataArray[ VdIterator , 6 ]
        
        # Adjusted values
        _Gdiff_complex      = _GdiffAmplitude*np.exp(1j*(_GdiffPhase+PhaseShift)/180*np.pi)+1E-10
        _GdiffDUT_complex   = 1 / ( alpha / _Gdiff_complex - 2*Zfilters ) # Gdiff in S
        _VdDUT              = _Vd -  0.001* 2*Rfilters * _Id  # Factor 0.001: Vd is in mV and Id in µA
        
        # Return the adjusted data array
        AdjustedDataArray[ VdIterator , 0:7 ] = [ _VdDUT , _Vd ,
                            np.abs(_GdiffDUT_complex), np.abs(_Gdiff_complex) , _ILeak , _Vg , _T ]
    return AdjustedDataArray
    
def remove_filters( _DataArray ):
    # Create an empty array of the same size as the original data array,
    # - for the adjusted values:
    AdjustedDataMap         = np.zeros( ( len( _DataArray ) , len ( _DataArray[0] ) ,
                                      len( _DataArray[0,0] ) ) )
    # - and for the interpolated values, with ten times smaller Vd steps: 
    InterpolatedData        = np.zeros( ( len( _DataArray ) , 
                                          len( _DataArray[0] ) * InterpolationOverhead,
                                          len( _DataArray[0,0] ) ) )
    for VgIterator in range( 0 , len( _DataArray ) ):
        # Return the adjusted data array
        AdjustedDataMap[ VgIterator , : , : ] = remove_filters1D( _DataArray[ VgIterator ] )
        
    # The Vd values in AdjustedDataMap are no longer equally spaced.
    # We are going to interpolate between them, and then evaluate that
    # interpolation at evenly spaced Vd values.

    # The Vg values (column 5) of the output are the same as the input
    InterpolatedData[:,:,5] = np.array( [ _DataArray[:,0,5] ] * 
                                len( InterpolatedData[0] ) ).T
    # The output Vd range is given by AdjustedDataMap, and has more points
    InterpolationRangeMin   = min( [ min( AdjustedDataMap[VgIterator,:,0] ) for VgIterator
                                     in range( 0 , len(AdjustedDataMap) ) ] )
    InterpolationRangeMax   = max( [ max( AdjustedDataMap[VgIterator,:,0] ) for VgIterator
                                     in range( 0 , len(AdjustedDataMap) ) ] )
    InterpolatedData[:,:,0] = np.linspace(  InterpolationRangeMin ,
                                            InterpolationRangeMax ,
                                            num=len( InterpolatedData[0] ) ,
                                            endpoint=True )
    # The other columns will be interpolated
    for VgIterator in range( 0 , len( _DataArray ) ):
        AdjustedVd  = AdjustedDataMap[ VgIterator , : , 0 ]
        for Column in [1, 2, 3, 4, 6]:
            AdjustedColumn                  = AdjustedDataMap[ VgIterator , : , Column ]
            InterpolatedFunction            = interp1d( AdjustedVd , AdjustedColumn ,
                                                       kind='nearest',
                                                       bounds_error = False ,
                                                       fill_value=0)
            InterpolatedData[ VgIterator , 
                             : ,Column]     = InterpolatedFunction( 
                                                         InterpolatedData[0,:,0] )
    return InterpolatedData

def groupby_mean(a):
    # Sort array by groupby column
    b = a[a[:,0].argsort()]
    # Get interval indices for the sorted groupby col
    idx = np.flatnonzero(np.r_[True,b[:-1,0]!=b[1:,0],True])
    # Get counts of each group and sum rows based on the groupings & hence averages
    counts = np.diff(idx)
    avg = np.add.reduceat(b[:,1:],idx[:-1],axis=0)/counts.astype(float)[:,None]
    # Finally concatenate for the output in desired format
    return np.c_[b[idx[:-1],0],avg]

def imshow( x, y, data ):
    if SweepType == "VgVd":
        XLabel  = "Gate Voltage (V)"
    elif SweepType == "TVd":
        XLabel  = "T (K)"
    elif SweepType == "HVd":
        XLabel  = "H (mT)"
    else:
        print("Error: SweepType not found")
        exit
    win = ImageDialog(edit=False, toolbar=True,
                      wintitle="Image with custom X/Y axes scales",
                      options=dict(xlabel=XLabel, ylabel="Drain Voltage (mV)",
                                   yreverse=False,show_contrast=True))
    item = make.xyimage(x, y, data)
    plot = win.get_plot()
    plot.set_aspect_ratio(lock=False)
    plot.add_item(item)
    win.show()
    win.exec_()
    return plot

def selectData(  x, y, data, xRange, yRange ):
    xSelection      = x[bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    ySelection      = y[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1])]
    dataSelection   = data[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1]),
                           bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    return np.array([ xSelection , ySelection, dataSelection ])

def imshowRange( x, y, data, xRange, yRange ):
    xSelection      = x[bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    ySelection      = y[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1])]
    dataSelection   = data[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1]),
                           bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    return imshow( xSelection , ySelection, dataSelection )

if User == "Francois":
    # Francois stores the data on his C drive
    RootDirectory   = "C:/Users/FL134692/Documents/ANR - Projets/2019 SUNISIDEUP/Data/C2N_PtSi/"
elif User == "Tom":
    # Tom stores the data on his Z drive
    RootDirectory   = "Z:/Data/Data sorted by lot/W12 Laurie Calvet's wafer/Christophe's fridge/"
    RootDirectoryBigoudene  = "Z:/Data/Data sorted by lot/W12 Laurie Calvet's wafer/Bigoudene/"
else:
    print("Error: User not found")
    exit

#%% IMPORT

# Single file:
DataArray300K   = np.array( np.loadtxt( RootDirectoryBigoudene + DataFileList[0] ) )
DataArray4K     = np.array( np.loadtxt( RootDirectoryBigoudene + DataFileList[1] ) )
DataArray150mK  = np.array( np.loadtxt( RootDirectoryBigoudene + DataFileList[2] ) )

#%% PLOT 1D CURVES

plt.ioff()

# Number 1
fig     = plt.figure(figsize=(4,3.5))
plt.plot( groupby_mean( DataArray300K )[:,0], 1000* DataArray300K[0,1] - 1000* savgol_filter( groupby_mean( DataArray300K )[:,1], 3, 1 ) ,
         label='300 K')
plt.xlabel('$V_\mathrm{g}$ (V)',fontsize=13);
plt.ylabel('$I_d$ (nA)  at  $V_\mathrm{d}=1$ mV',fontsize=13);
plt.gca().set_xlim([-4,0])
plt.legend()
plt.tight_layout()
plt.savefig('img/D84D4_Id_Vg_300K_4K_150mK_1.png',dpi=500,bbox_inches='tight')

# Number 2
fig     = plt.figure(figsize=(4,3.5))
plt.plot( groupby_mean( DataArray300K )[:,0], 1000* DataArray300K[0,1] - 1000* savgol_filter( groupby_mean( DataArray300K )[:,1], 3, 1 ) ,
         label='300 K')
plt.plot( groupby_mean( DataArray4K )[:,0] ,  - 1000* savgol_filter( groupby_mean( DataArray4K )[:,1] , 3 , 1 ) ,
         label='4 K')
plt.xlabel('$V_\mathrm{g}$ (V)',fontsize=13);
plt.ylabel('$I_d$ (nA)  at  $V_\mathrm{d}=1$ mV',fontsize=13);
plt.gca().set_xlim([-4,0])
plt.legend()
plt.tight_layout()
plt.savefig('img/D84D4_Id_Vg_300K_4K_150mK_2.png',dpi=500,bbox_inches='tight')

# Number 3
fig     = plt.figure(figsize=(4,3.5))
plt.plot( groupby_mean( DataArray300K )[:,0], 1000* DataArray300K[0,1] - 1000* savgol_filter( groupby_mean( DataArray300K )[:,1], 3, 1 ) ,
         label='300 K')
plt.plot( groupby_mean( DataArray4K )[:,0] ,  - 1000* savgol_filter( groupby_mean( DataArray4K )[:,1] , 3 , 1 ) ,
         label='4 K')
plt.plot( groupby_mean(DataArray150mK)[:,0] , 1000* groupby_mean(DataArray150mK)[:,1] ,
         label='150 mK')
plt.xlabel('$V_\mathrm{g}$ (V)',fontsize=13);
plt.ylabel('$I_d$ (nA)  at  $V_\mathrm{d}=1$ mV',fontsize=13);
plt.gca().set_xlim([-4,0])
plt.legend()
plt.tight_layout()
plt.savefig('img/D84D4_Id_Vg_300K_4K_150mK_3.png',dpi=500,bbox_inches='tight')