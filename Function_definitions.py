# -*- coding: utf-8 -*-
"""
Created on Wed Apr 07 15:30:42 2021

@author: TV255140
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter
rcParams.update({'figure.autolayout': True})

import bisect
# import guidata                        # guidata
# from guiqwt.plot import ImageDialog   # ImageDialog does ???
# from guiqwt.builder import make       # make does ???
# from guidata.qt.QtGui import QFont
import numpy as np                    # numpy gives us numerical analysis functions
import os                             # os, "miscellaneous operating system interfaces", is for handling directory/file operations. Useful for getting the files in a directory
import re                             # re for regular expressions (https://docs.python.org/3/howto/regex.html)
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

# _app        = guidata.qapplication()

#%% Operations specific to our type of data

def GetDataArray( DataFileList, Directory ):
    return np.array( [ np.loadtxt( Directory + Filename ) for Filename in DataFileList ])

def GetDataFileList( Directory, SweepType ):
     return np.array( [ entry for entry in os.listdir( Directory )[1*bool(SweepType=="TVd" or SweepType=="HVd")::] ] )

def GetRootDirectory( User ):
    if User == "Francois":
        # Francois stores the data on his C drive
        RootDirectory   = "C:/Users/FL134692/Documents/ANR - Projets/2019 SUNISIDEUP/Data/C2N_PtSi/"
    elif User == "Tom":
        # Tom stores the data on his Z drive
        RootDirectory   = "C:/Users/vethaak/OneDrive - Chalmers/Documents/PhD/Data/Data sorted by lot/W12 Laurie Calvet's wafer/"
    else:
        print("Error: User not found")
        exit
    return RootDirectory

def GetXValues( SweepType, SweepStart, SweepStep, DataFileList, DataArray ):
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
    return XValues
    
def remove_filters1D( Global, _DataArray ):
    # Create an empty array of the same size as the original data array,
    # - for the adjusted values:
    AdjustedDataArray       = np.zeros( ( len( _DataArray ) ,
                                      len( _DataArray[0] ) ) )
    for VdIterator in range( 0 , len( _DataArray ) ):
        # Values given in the data files
        _Vd             = _DataArray[ VdIterator , 0 ] # mV
        _Id             = _DataArray[ VdIterator , 1 ] # µA
        _GdiffAmplitude = _DataArray[ VdIterator , 2 ] # S
        _GdiffPhase     = _DataArray[ VdIterator , 3 ] #
        _ILeak          = _DataArray[ VdIterator , 4 ] # A
        _Vg             = _DataArray[ VdIterator , 5 ] # V
        _T              = _DataArray[ VdIterator , 6 ] # K
        
        # Adjusted values
        _Gdiff_complex      = _GdiffAmplitude * np.exp( 1j * ( _GdiffPhase + Global['PhaseShift']) / 180 * np.pi )
        _GdiffDUT_complex   = 1 / ( Global['alpha'] / _Gdiff_complex - 2 * Global['Zfilters'] ) # Gdiff in S
        _VdDUT              = _Vd -  0.001* 2 * Global['Rfilters'] * _Id  # Factor 0.001: Vd is in mV and Id in µA
        
        # Return the adjusted data array
        AdjustedDataArray[ VdIterator , 0:7 ] = [ _VdDUT , _Id ,
                            np.abs(_GdiffDUT_complex), np.abs( _Gdiff_complex ), _ILeak , _Vg , _T ]
    return AdjustedDataArray
    
def remove_filters( Global, _DataArray ):
    # Create an empty array of the same size as the original data array,
    # - for the adjusted values:
    AdjustedDataMap         = np.zeros( ( len( _DataArray ) , len ( _DataArray[0] ) ,
                                      len( _DataArray[0,0] ) ) )
    # - and for the interpolated values, with ten times smaller Vd steps: 
    InterpolatedData        = np.zeros( ( len( _DataArray ) , 
                                          len( _DataArray[0] ) * Global['InterpolationOverhead'],
                                          len( _DataArray[0,0] ) ) )
    for VgIterator in range( 0 , len( _DataArray ) ):
        # Return the adjusted data array
        AdjustedDataMap[ VgIterator , : , : ] = remove_filters1D( Global, _DataArray[ VgIterator ] )
        
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

def selectData( x, y, data, xRange, yRange ):
    xSelection      = x[bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    ySelection      = y[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1])]
    dataSelection   = data[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1]),
                           bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    return np.array([ xSelection , ySelection, dataSelection ], dtype=object )

def Select_Data_XRange( DataArray, XRange ):
    return DataArray[bisect.bisect_left(DataArray[:,0,5], XRange[0]):bisect.bisect_right(DataArray[:,0,5], XRange[1])]

def imshowRange( x, y, data, xRange, yRange ):
    xSelection      = x[bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    ySelection      = y[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1])]
    dataSelection   = data[bisect.bisect_left(y,yRange[0]):bisect.bisect_right(y,yRange[1]),
                           bisect.bisect_left(x,xRange[0]):bisect.bisect_right(x,xRange[1])]
    return imshow( xSelection , ySelection, dataSelection )

#%% Return an adjusted data set

def DataArray_To_2DMap( Global, DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, SetupOrDUT ):
    # Select raw data in one direction: avoid features from averaging
    (DataArraySelectionY,DataArraySelectionX,DataArraySelection)= selectData( DataArray[0,:,0], DataArray[:,0,5],   DataArray[ : , :int(len( DataArray[0,:] ) / 2) , : ],  DataYRange,   DataXRange)
    # Filter the data    
    DataArrayAdjusted       = remove_filters( Global, DataArraySelection )
    VdAdjusted              = DataArrayAdjusted[ 0 , :, 0 ] - YOffset
    # Select filtered data:
    if SetupOrDUT == 'SETUP':
        GdiffAdjusted_SETUP     = np.swapaxes( DataArrayAdjusted[ : , : , 3 ] , 0 , 1 )
        (VgSelection,VdSelection,GdiffSETUPDataSelection)   = selectData(DataArraySelectionX,VdAdjusted,GdiffAdjusted_SETUP ,XRange,YRange)
        return GdiffSETUPDataSelection
    else:
        GdiffAdjusted_DUT       = np.swapaxes( DataArrayAdjusted[ : , : , 2 ] , 0 , 1 )
        (VgSelection,VdSelection,GdiffDUTDataSelection)     = selectData(DataArraySelectionX,VdAdjusted,GdiffAdjusted_DUT   ,XRange,YRange)
        return GdiffDUTDataSelection

#%% Misc functions
    
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
  
#%% List operations

def Select_At_H( DataArray, X ):
    if min( DataArray[:,0,5] ) <= X and max( DataArray[:,0,5] ) >= X:
        return DataArray[ bisect.bisect_left( DataArray[:,0,5], X ), :]
    else:
        return np.array([[0]*len(DataArray[0,0] )]*len(DataArray[0]))

def Select_At_Vg( DataArray, X ):
    if min( DataArray[:,0,5] ) <= X and max( DataArray[:,0,5] ) >= X:
        return DataArray[ bisect.bisect_left( DataArray[:,0,5], X ), :]
    else:
        return np.array([[0]*len(DataArray[0,0] )]*len(DataArray[0]))

def Select_At_X( DataArray, X ):
    SortedDataArray     = groupby_mean( DataArray )
    if min( SortedDataArray[:,0] ) < X and max( SortedDataArray[:,0] ) > X:
        return SortedDataArray[ bisect.bisect_left( SortedDataArray[:,0], X ), :]
    else:
        return [0]*len( SortedDataArray[0] )

#%% BTK

from scipy.optimize import minimize_scalar

def BTK_Normalized_Isn_At_Delta_core( Z ):
    if Z == 0:
        return 2
    else:
        return np.sqrt( 1 + Z**2 ) / ( Z*( 1 + 2*Z**2 ) ) * np.arctanh( 2*Z*np.sqrt( 1 + Z**2 ) / (  1 + 2*Z**2 ) )

def BTK_Normalized_Isn_At_Delta( Z ):
    if isinstance( Z, list ) or isinstance( Z, np.ndarray ):
        OutputArray             = np.array( [0.] * len( Z ) )
        for ii in np.arange(0, len(Z) ):
            OutputArray[ii]     = BTK_Normalized_Isn_At_Delta_core( Z[ii] )
        return OutputArray
    elif isinstance( Z, float ) or isinstance( Z, int ):
        return BTK_Normalized_Isn_At_Delta_core( Z )
    else:
        return "ERROR: data type not recognized"

def BTK_Z_Error( Z, Isn ):
    return np.abs( Isn - BTK_Normalized_Isn_At_Delta( Z ) )

def BTK_Isn_to_Z_core( Isn ):
    Minimize_Scalar_Z   = minimize_scalar( BTK_Z_Error , bounds=(0,100), args = ( Isn ), method='bounded')
    return Minimize_Scalar_Z.x

def BTK_Isn_to_Z( Isn ):
    if isinstance( Isn, list ) or isinstance( Isn, np.ndarray ):
        OutputArray             = np.array( [0.] * len( Isn ) )
        for ii in np.arange(0, len(Isn) ):
            OutputArray[ii]     = BTK_Isn_to_Z_core( Isn[ii] )
        return OutputArray
    elif isinstance( Isn, float ) or isinstance( Isn, int ):
        return BTK_Isn_to_Z_core( Isn )
    else:
        return "ERROR: data type not recognized"

def Z_to_Gamma( Z ):
    return 1 / ( 1 + Z**2 )

def Z_to_BTK_A( Z ):
    return 1 / ( 1 + 2*Z**2 )**2

def Shift_Badness( Shift, DataArray, Sign ):
    ShiftedArray    = DataArray+np.array([[Shift],[0]])
    PositivePart    = ShiftedArray[:,bisect.bisect_left(ShiftedArray[0],0):]
    NegativePart    = ShiftedArray[:,0:bisect.bisect_left(ShiftedArray[0],0)]
    AvgCurve        = np.concatenate((PositivePart,NegativePart*np.array([[-1],[Sign]])),axis=1)
    return np.sum(np.abs(np.diff(AvgCurve)))

def Find_Sign( DataArray ):
    LeftPart        = DataArray[:,bisect.bisect_left(DataArray[0],0):]
    RightPart       = DataArray[:,0:bisect.bisect_right(DataArray[0],0)]
    return np.sign(np.mean(LeftPart[1])/np.mean(RightPart[1]))

def Curve_Shift( DataArray ):
    Sign            = Find_Sign( DataArray )
    Minimize_Scalar_Shift   = minimize_scalar( Shift_Badness , bounds=(-0.1,0.1), args = ( DataArray, Sign ), method='bounded')
    return Minimize_Scalar_Shift.x

def Curve_Shift_Isn_3D( DataArray ):
    ShiftArray  = np.array([0.]*len(DataArray))
    for VgIterator in np.arange(0,len(DataArray)):
        ShiftArray[VgIterator]  = Curve_Shift( DataArray[VgIterator,[0,1]] )
    return np.mean( ShiftArray )

def Curve_Shift_Gdiff_3D( DataArray ):
    ShiftArray  = np.array([0.]*len(DataArray))
    for VgIterator in np.arange(0,len(DataArray)):
        ShiftArray[VgIterator]  = Curve_Shift( DataArray[VgIterator,[0,2]] )
    return np.mean( ShiftArray )

def Center_Curve( DataArray ):
    XShift          = Curve_Shift( DataArray )
    return DataArray+np.array([[XShift],[0]])

def Center_Array( DataArray ):
    XShift          = Curve_Shift( np.transpose(groupby_mean(DataArray)[:,[0,2]]) )
    return DataArray+np.transpose(np.array(np.concatenate([[[XShift]],[[0]]*( len( DataArray[0] ) - 1 )])))

def Mirror_Curve( DataArray ):
    CentredDataArray= Center_Curve( DataArray )
    PositivePart    = CentredDataArray[:,bisect.bisect_left(DataArray[0],0):]
    NegativePart    = CentredDataArray[:,0:bisect.bisect_right(DataArray[0],0)]
    Sign            = Find_Sign( DataArray )
    MirrorCurve     = np.concatenate((PositivePart,NegativePart*np.array([[-1],[Sign]])),axis=1)
    return MirrorCurve

def Extract_RN( DataArray, X ):
    INLeft          = np.abs( Select_At_X( DataArray, -X )[1] )
    INRight         = np.abs( Select_At_X( DataArray, X )[1] )
    if INLeft > 0 and INRight > 0:
        return 1000 * X / np.mean([ INLeft, INRight ]) # Factor 1000: Vd is in mV and Id in µA
    else:
        return 0

def Extract_RNdiff( DataArray2D, Delta, RelativeRange ):
    SortedGdiffCurve    = groupby_mean( DataArray2D[:,[0,2]] )
    SmoothedGdiffCurve  = savgol_filter( ( SortedGdiffCurve[:,0], SortedGdiffCurve[:,1] ), 9, 2 )
    CenteredGdiffCurve  = Center_Curve( SmoothedGdiffCurve )
    LeftSelection       = CenteredGdiffCurve[:,bisect.bisect_left(CenteredGdiffCurve[0],-(1+RelativeRange)*Delta):bisect.bisect_right(CenteredGdiffCurve[0],-Delta)]
    RightSelection      = CenteredGdiffCurve[:,bisect.bisect_left(CenteredGdiffCurve[0],Delta):bisect.bisect_right(CenteredGdiffCurve[0],(1+RelativeRange)*Delta)]
    if len(LeftSelection[1]) > 0 and len(RightSelection[1]) > 0:
        LeftGN          = np.min(LeftSelection[1])
        RightGN         = np.min(RightSelection[1])
        if LeftGN != 0 and RightGN != 0:
            LeftRN      = 1/LeftGN
            RightRN     = 1/RightGN
        else:
            LeftRN      = 0
            RightRN     = 0
    else:
        LeftRN          = 0
        RightRN         = 0
    return np.mean(np.array([LeftRN,RightRN]))

def Extract_Isn( IsnCurve, Delta ):
    SortedIsnCurve      = groupby_mean( IsnCurve )
    SmoothedIsnCurve    = savgol_filter( ( SortedIsnCurve[:,0], SortedIsnCurve[:,1] ), 9, 2 )
    CenteredIsnCurve    = Center_Curve( SmoothedIsnCurve )
    if min(CenteredIsnCurve[0]) < -Delta and max(CenteredIsnCurve[0]) > Delta:
        LeftIsn         = CenteredIsnCurve[1,bisect.bisect_left(CenteredIsnCurve[0],-Delta)]
        RightIsn        = CenteredIsnCurve[1,bisect.bisect_left(CenteredIsnCurve[0],Delta)]
    else:
        LeftIsn         = 0
        RightIsn        = 0
    return np.mean(np.array([np.abs(LeftIsn),np.abs(RightIsn)]))