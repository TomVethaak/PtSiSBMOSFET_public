# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 10:55:18 2021

@author: TV255140
"""

from Function_definitions import *

GlobalDPI   = 600

#%% Initialize colors
rcParams['image.cmap']              ='viridis'  # Color map
rcParams['axes.linewidth']          = 1.2      # frame thickness
plt.rcParams['axes.facecolor']      ='white'    # background color
plt.rcParams['savefig.facecolor']   ='white'    # background color

#%% 2D plot Gdiff

def Plot_2D_Gdiff( Global, IMGfilename, DataArray, VgValues, DataXRange, DataYRange, XRange, YRange, YOffset, SetupOrDUT ):
    DataSet2DPlot   = DataArray_To_2DMap( Global, DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, SetupOrDUT )  # Removes filters
    fig             = plt.figure(figsize=(10,5))
    plt.imshow(1E6*DataSet2DPlot, aspect='auto',
               extent=[XRange[0],   XRange[1],
                       YRange[0],   YRange[1]],
               interpolation = 'none');
    plt.xlabel( '$V_\mathrm{g}$ (V)',   fontsize=11 );
    plt.ylabel( '$V_\mathrm{d}$ (mV)',  fontsize=11 );
    cbar    = plt.colorbar(shrink=0.45,pad=0.02)
    cbar.set_label('$G_\mathrm{diff}\;(\mu\mathrm{S})$', fontsize=11, rotation=90)
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Combo plot of 2D Gdiff and Gdiff selected curves

def ComboPlot_Gdiff( Global, IMGfilename, DataArray, VgValues, DataXRange, DataYRange, XRange0, XRange10, XRange11, YRange0, YRange10, YRange11, YOffset, SetupOrDUT ):
    # 2D map
    DataSet2DPlot   = DataArray_To_2DMap( Global, DataArray, DataXRange, DataYRange, XRange0, YRange0, YOffset, SetupOrDUT )
    #VgValues        = DataArray[:,0,5]
    #VdValues        = np.arange(YRange0[0],YRange0[1],(YRange0[1]-YRange0[0])/len(DataSet2DPlot))

    fig             = plt.figure(figsize=(10,5));
    gs              = GridSpec(nrows=2, ncols=2);
    axes0           = fig.add_subplot(gs[0,:]);
    axes10          = fig.add_subplot(gs[1,0]);
    axes11          = fig.add_subplot(gs[1,1]);
    axes0.imshow(DataSet2DPlot, aspect='auto' ,
                 extent=[XRange0[0]-0.5*(XRange0[1]-XRange0[0])/len(DataArray), XRange0[1]+0.5*(XRange0[1]-XRange0[0])/len(DataArray),
                         YRange0[0], YRange0[1]],
                 interpolation = 'none');
    axes0.set_xlabel( '$V_\mathrm{g}$ (V)',     fontsize=11);
    axes0.set_ylabel( '$V_\mathrm{d}$ (mV)',    fontsize=11);

    axes10.set_xlim(XRange10);
    axes10.set_ylim(YRange10);
    axes10.set_xlabel('$V_\mathrm{d}$ (mV)',        fontsize=11);
    axes10.set_ylabel('$G_\mathrm{diff}$ ($\mu$S)', fontsize=11);
    axes10.axvline(x = 0,       color = 'black',    linestyle = '--');
    axes10.axvline(x = -0.31,   color = 'grey',     linestyle = '--');
    axes10.axvline(x = 0.31,    color = 'grey',     linestyle = '--');

    axes11.set_xlim(XRange11);
    axes11.set_ylim(YRange11);
    axes11.set_xlabel('$V_\mathrm{d}$ (mV)',fontsize=11);
    axes11.axvline(x = 0, color = 'black', linestyle = '--');
    for VgIndex in np.arange(0,len(VgValues)):
        # Draw a vertical line where the cut is made
        axes0.axvline(x=VgValues[VgIndex],color=plt.rcParams['axes.prop_cycle'].by_key()['color'][VgIndex]);#['blue','orange','green'][VgIndex]);
        # Select one Gdiff(vd) curve at Vg=VgValue from the original data array
        DataArraySliceAtVg          = Select_At_Vg(DataArray,VgValues[VgIndex])
        # Process the data: remove the filters
        # Only take the first half of the data points to avoid double ZBCP
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )[int(len(DataArraySliceAtVg)/2):]
        # Apply Savitzky-Golay filter
        FilteredCurve               = savgol_filter( ( groupby_mean(DataArraySliceAtVgAdjusted)[:,0] ,
                                                       groupby_mean(DataArraySliceAtVgAdjusted)[:,2] ) , 3, 1)
        # Plot the Vd values vs the Gdiff values
        axes10.plot( FilteredCurve[0] - YOffset, 1E6*FilteredCurve[1],label='{:3.2f} V'.format(VgValues[VgIndex]));
        axes11.plot( FilteredCurve[0] - YOffset, 1E6*FilteredCurve[1],label='{:3.2f} V'.format(VgValues[VgIndex]));
    handles, labels = plt.gca().get_legend_handles_labels()
    axes10.legend(handles[::-1], labels[::-1],loc='upper left')
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Combo plot of 2D Gdiff and Gamma, BTK A

def ComboPlot_Gdiff_Gamma_A( Global, IMGfilename, DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, SetupOrDUT):
    # 2D map
    DataSet2DPlot   = DataArray_To_2DMap( Global, DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, SetupOrDUT )
    # Plotting
    VgValues        = DataArray[:,0,5]
    VdValues        = np.arange(YRange[0],YRange[1],(YRange[1]-YRange[0])/len(DataSet2DPlot))
    AArray          = np.array([0.]*len(VgValues))
    ZArray          = np.array([0.]*len(VgValues))
    GammaArray      = np.array([0.]*len(VgValues))
    for VgIndex in np.arange(0,len(VgValues)):
        # Select one Gdiff(vd) curve at Vg=VgValue from the original data array
        DataArraySliceAtVg          = Select_At_Vg(DataArray,VgValues[VgIndex])
        # Process the data: remove the filters
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )
        # Extract RN
        RN      = Extract_RNdiff( DataArraySliceAtVgAdjusted, Global['DeltaPtSi'] , 0.1 )
        # Extract Isn
        Isn     = Extract_Isn( DataArraySliceAtVgAdjusted[:,[0,1]], Global['DeltaPtSi'] )
        # Calculate normalized Isn
        NormalizedIsn   = 0.001 * Isn * RN / Global['DeltaPtSi'] # Factor 0.001: Delta is in mV, Isn in muA
        # Solve for Z
        ZArray[VgIndex] = BTK_Isn_to_Z( NormalizedIsn )
        # Calculate Gamma from Z
        GammaArray[VgIndex] = Z_to_Gamma( ZArray[VgIndex] )
        # Calculate A from Z
        AArray[VgIndex]     = Z_to_BTK_A( ZArray[VgIndex] )
    fig, axes   = plt.subplots(2, 1, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(8,5))
    axes[0].set_xticklabels([])
    axes[0].imshow(DataSet2DPlot, aspect='auto' ,
                 extent=[np.min(VgValues),  np.max(VgValues),
                         np.min(VdValues),  np.max(VdValues)],
                 interpolation = 'none');
    axes[0].set_ylabel('$V_\mathrm{d}$ (mV)',fontsize=11);
    axes[1].set_xlim(XRange);
    axes[1].set_ylim([0,max(GammaArray)]);
    axes[1].set_xlabel('$V_\mathrm{g}$ (V)',fontsize=11);
    axes[1].set_ylabel('$A$, $\Gamma$',fontsize=11);
    axes[1].plot( VgValues, GammaArray, label='$\Gamma (V_\mathrm{g})$' );
    axes[1].plot( VgValues, AArray,     label='$A (V_\mathrm{g}, V_\mathrm{d}=0)$' );
    handles, labels = plt.gca().get_legend_handles_labels();
    axes[1].legend(handles, labels,loc='upper right');
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Combo plot of 2D Gdiff and Gamma, BTK A (SMALL)

def ComboPlot_Gdiff_Gamma_A_small( Global, IMGfilename, DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, SetupOrDUT ):
    # 2D map
    DataSet2DPlot   = DataArray_To_2DMap( Global, DataArray, DataXRange, DataYRange, XRange, YRange, YOffset, SetupOrDUT )
    # Plotting
    VgValues        = DataArray[:,0,5]
    VdValues        = np.arange(YRange[0],YRange[1],(YRange[1]-YRange[0])/len(DataSet2DPlot))
    AArray          = np.array([0.]*len(VgValues))
    ZArray          = np.array([0.]*len(VgValues))
    GammaArray      = np.array([0.]*len(VgValues))
    for VgIndex in np.arange(0,len(VgValues)):
        # Select one Gdiff(vd) curve at Vg=VgValue from the original data array
        DataArraySliceAtVg          = Select_At_Vg(DataArray,VgValues[VgIndex])
        # Process the data: remove the filters
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )
        # Extract RN
        RN      = Extract_RNdiff( DataArraySliceAtVgAdjusted, Global['DeltaPtSi'] , 0.1 )
        # Extract Isn
        Isn     = Extract_Isn( DataArraySliceAtVgAdjusted[:,[0,1]], Global['DeltaPtSi'] )
        # Calculate normalized Isn
        NormalizedIsn   = 0.001 * Isn * RN / Global['DeltaPtSi'] # Factor 0.001: Delta is in mV, Isn in muA
        # Solve for Z
        ZArray[VgIndex] = BTK_Isn_to_Z( NormalizedIsn )
        # Calculate Gamma from Z
        GammaArray[VgIndex] = Z_to_Gamma( ZArray[VgIndex] )
        # Calculate A from Z
        AArray[VgIndex]     = Z_to_BTK_A( ZArray[VgIndex] )
    fig, axes   = plt.subplots(2, 1, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(5,5))
    axes[0].set_xticklabels([])
    axes[0].imshow(DataSet2DPlot, aspect='auto' ,
                 extent=[np.min(VgValues),  np.max(VgValues),
                         np.min(VdValues),  np.max(VdValues)],
                 interpolation = 'none');
    axes[0].set_ylabel('$V_\mathrm{d}$ (mV)',fontsize=11);
    axes[1].set_xlim(XRange);
    axes[1].set_ylim([0,max(GammaArray)]);
    axes[1].set_xlabel('$V_\mathrm{g}$ (V)',fontsize=11);
    axes[1].set_ylabel('$A$, $\Gamma$',fontsize=11);
    axes[1].plot( VgValues, GammaArray, label='$\Gamma (V_\mathrm{g})$' );
    axes[1].plot( VgValues, AArray,     label='$A (V_\mathrm{g}, V_\mathrm{d}=0)$' );
    handles, labels = plt.gca().get_legend_handles_labels();
    axes[1].legend(handles, labels,loc='upper right');
    plt.text(0.035, 0.08, '(a)',  transform = axes[0].transAxes,   fontsize=11)
    plt.text(0.035, 0.20, '(b)',  transform = axes[1].transAxes,   fontsize=11)
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% ListPlot: plot a list of curves

def ListPlot_Gdiff_Vg_Manual( Global, IMGfilename, CurveList, LabelList, XRange, YRange, XShifts ):
    fig     = plt.figure(figsize=(5,3.7));
    plt.xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11);
    plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    plt.axvline(x = 0,  color = 'black',    linestyle = '--');
    plt.axvline(x = -Global['DeltaPtSi'], color = 'grey',     linestyle = '--');
    plt.axvline(x = Global['DeltaPtSi'],  color = 'grey',     linestyle = '--');
    plt.xlim( [ -0.39, 0.39 ] )
    plt.ylim( [ 0, 5.6E1 ] )
    for CurveIndex in np.arange( 0, len( CurveList ) ):
        plt.plot( CurveList[CurveIndex][0] + XShifts[CurveIndex], 1E6*CurveList[CurveIndex][1], label=LabelList[ CurveIndex ] )
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend( handles, labels, loc='upper center' )
    plt.savefig( IMGfilename, dpi=GlobalDPI, bbox_inches='tight' )
    plt.close()

#%% Plot Gdiff curves at selected Vg values

def Plot_Gdiff_Vg( Global, IMGfilename, DataArray, VgValues, XRange, YRange ):
    fig     = plt.figure(figsize=(5,3));
    plt.xlim(XRange);
    plt.ylim(YRange);
    plt.xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11);
    plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    plt.axvline(x = 0,      color = 'black',    linestyle = '--');
    plt.axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    plt.axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    for Vg in VgValues:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtVg          = Select_At_Vg(DataArray,Vg)
        # Process the data: remove the filters
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtVgAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = Center_Curve( savgol_filter( ( groupby_mean(
                                        DataArraySliceAtVgAdjusted )[:,0] , SortedCurve ) , 3, 1) )
        # Plot the Vd values vs the Gdiff values
        plt.plot(FilteredCurve[0] , 1E6*FilteredCurve[1], label='{:3.2f} V'.format(Vg));
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1],loc='upper left')
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Plot Gdiff curves at selected Vg values

def Plot_Gdiff_Vg_Manual( Global, IMGfilename, DataArray, VgValues, XRange, YRange, XShifts ):
    fig     = plt.figure(figsize=(5,3.7));
    plt.xlim(XRange);
    plt.ylim(YRange);
    plt.xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11);
    plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    plt.axvline(x = 0,      color = 'black',    linestyle = '--');
    plt.axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    plt.axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    for Vg in VgValues:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtVg          = Select_At_Vg(DataArray,Vg)
        # Process the data: remove the filters
        # Take only the first half of the data to avoid double ZBCP
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )#[int(len(DataArraySliceAtVg)/2):]
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtVgAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = Center_Curve( savgol_filter( ( groupby_mean(
                                        DataArraySliceAtVgAdjusted )[:,0] , SortedCurve ) , 5, 1) )
        # Plot the Vd values vs the Gdiff values
        plt.plot(FilteredCurve[0] + XShifts[VgValues.index(Vg)], 1E6*FilteredCurve[1], label='{:3.2f} V'.format(Vg));
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend( handles[::-1], labels[::-1], loc='upper center' )
    plt.savefig( IMGfilename, dpi=GlobalDPI, bbox_inches='tight' );
    plt.close();

#%% Plot Gdiff curves at selected H values

def Plot_Gdiff_H( Global, IMGfilename, DataArray, HValues, HSelection, XRange, YRange ):
    fig             = plt.figure(figsize=(4,2.5));
    plt.xlim(XRange);
    plt.ylim(YRange);
    plt.xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11);
    plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    plt.axvline(x = 0,      color = 'black',    linestyle = '--');
    plt.axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    plt.axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    for H in HSelection:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtH           = DataArray[ bisect.bisect_left( HValues, H ), :, : ]
        # Process the data: remove the filters
        DataArraySliceAtHAdjusted   = remove_filters1D( Global, DataArraySliceAtH )
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtHAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = Center_Curve( savgol_filter( ( groupby_mean(
                                        DataArraySliceAtHAdjusted )[:,0] , SortedCurve ) , 3, 1) )
        # Plot the Vd values vs the Gdiff values
        plt.plot(FilteredCurve[0] , 1E6*FilteredCurve[1], label='{:2.0f} mT'.format( 1E3 * H ));
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1],loc='upper left')
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Plot Gdiff curves at selected H values with manual X (Vd) shifts

def Plot_Gdiff_H_Manual( Global, IMGfilename, DataArray, HValues, HSelection, XRange, YRange, XShifts ):
    fig             = plt.figure(figsize=(4,2.5));
    plt.xlim(XRange);
    plt.ylim(YRange);
    plt.xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11);
    plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    plt.axvline(x = 0,      color = 'black',    linestyle = '--');
    plt.axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    plt.axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    for H in HSelection:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtH           = DataArray[ bisect.bisect_left( HValues, H ), :, : ]
        # Process the data: remove the filters
        DataArraySliceAtHAdjusted   = remove_filters1D( Global, DataArraySliceAtH )
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtHAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = savgol_filter( ( groupby_mean(
                                        DataArraySliceAtHAdjusted )[:,0] , SortedCurve ) , 3, 1)
        # Plot the Vd values vs the Gdiff values
        plt.plot( FilteredCurve[0] + XShifts[HSelection.index(H)], 1E6*FilteredCurve[1], label='{:2.0f} mT'.format( 1E3 * H ));
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1],loc='upper left')
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Plot Gdiff curves at selected T values

def Plot_Gdiff_T( Global, IMGfilename, DataArray, TValues, TSelection, XRange, YRange ):
    fig             = plt.figure(figsize=(4,2.5));
    plt.xlim(XRange);
    plt.ylim(YRange);
    plt.xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11);
    plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    plt.axvline(x = 0,      color = 'black',    linestyle = '--');
    plt.axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    plt.axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    for T in TSelection:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtT           = DataArray[ bisect.bisect_left( DataArray[:,0,6], T ), :, : ]
        # Process the data: remove the filters
        DataArraySliceAtTAdjusted   = remove_filters1D( Global, DataArraySliceAtT )
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtTAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = Center_Curve( savgol_filter( ( groupby_mean(
                                        DataArraySliceAtTAdjusted )[:,0] , SortedCurve ) , 3, 1) )
        # Plot the Vd values vs the Gdiff values
        plt.plot(FilteredCurve[0] , 1E6*FilteredCurve[1], label='{:3.2f} K'.format( T ));
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1],loc='upper left')
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Plot Gdiff curves at selected T values

def Plot_Gdiff_T_manual( Global, IMGfilename, DataArray, TValues, TSelection, XRange, YRange, XShifts):
    fig             = plt.figure(figsize=(4,2.5));
    plt.xlim(XRange);
    plt.ylim(YRange);
    plt.xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11);
    plt.ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    plt.axvline(x = 0,      color = 'black',    linestyle = '--');
    plt.axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    plt.axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    for T in TSelection:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtT           = DataArray[ bisect.bisect_left( DataArray[:,0,6], T ), :, : ]
        # Process the data: remove the filters
        DataArraySliceAtTAdjusted   = remove_filters1D( Global, DataArraySliceAtT )
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtTAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = savgol_filter( ( groupby_mean(
                                        DataArraySliceAtTAdjusted )[:,0] , SortedCurve ) , 3, 1)
        # Plot the Vd values vs the Gdiff values
        plt.plot(FilteredCurve[0] + XShifts[TSelection.index(T)] , 1E6*FilteredCurve[1], label='{:3.2f} K'.format( T ));
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1],loc='upper left')
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Combine plots of Gdiff curves at selected H and T values

def ComboPlot_Gdiff_H_T_manual( Global, IMGfilename, HDataArray, HValues, HSelection, XRange, HYRange, HXShifts,
                                                     TDataArray, TValues, TSelection,         TYRange, TXShifts):
    fig, axes   = plt.subplots(2, 1, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(5,5))

    # Top plot: scan in T
    axes[0].set_xticklabels([])
    axes[0].set_xlim(XRange);
    axes[0].set_ylim([ TYRange[0], TYRange[1] ]);
    axes[0].set_ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    axes[0].axvline(x = 0,      color = 'black',    linestyle = '--');
    axes[0].axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    axes[0].axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    handlelist=[]
    for T in TSelection:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtT           = TDataArray[ bisect.bisect_left( TDataArray[:,0,6], T ), :, : ]
        # Process the data: remove the filters
        DataArraySliceAtTAdjusted   = remove_filters1D( Global, DataArraySliceAtT )
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtTAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = savgol_filter( ( groupby_mean(
                                        DataArraySliceAtTAdjusted )[:,0] , SortedCurve ) , 3, 1)
        # Plot the Vd values vs the Gdiff values
        handlelist+=axes[0].plot(FilteredCurve[0] + TXShifts[TSelection.index(T)] , 1E6*FilteredCurve[1], label='{:3.2f} K'.format( T ));
    handles, labels = plt.gca().get_legend_handles_labels();
    axes[0].legend(handles=list(reversed(handlelist)), loc='upper left');

    # Bottom plot: scan in H
    axes[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axes[1].set_xlim(XRange);
    axes[1].set_ylim([ HYRange[0], HYRange[1] ]);
    axes[1].set_xlabel( '$V_\mathrm{d}$ (mV)',          fontsize=11 );
    axes[1].set_ylabel( '$G_\mathrm{diff}$ ($\mu$S)',   fontsize=11);
    axes[1].axvline(x = 0,      color = 'black',    linestyle = '--');
    axes[1].axvline(x = -0.31,  color = 'grey',     linestyle = '--');
    axes[1].axvline(x = 0.31,   color = 'grey',     linestyle = '--');
    handlelist=[]
    for H in HSelection:
        # Select one Gdiff(Vd) curve from the original data array
        DataArraySliceAtH           = HDataArray[ bisect.bisect_left( HValues, H ), :, : ]
        # Process the data: remove the filters
        DataArraySliceAtHAdjusted   = remove_filters1D( Global, DataArraySliceAtH )
        # Column 1 is Id, column 2 is GdiffDUT
        SortedCurve                 = groupby_mean( DataArraySliceAtHAdjusted )[:,2]
        # Apply Savitzky-Golay filter
        FilteredCurve               = savgol_filter( ( groupby_mean(
                                        DataArraySliceAtHAdjusted )[:,0] , SortedCurve ) , 3, 1)
        # Plot the Vd values vs the Gdiff values
        handlelist+=axes[1].plot( FilteredCurve[0] + HXShifts[HSelection.index(H)] , 1E6*FilteredCurve[1], label='{:2.0f} mT'.format( 1E3 * H ));
    handles, labels = plt.gca().get_legend_handles_labels();
    axes[1].legend(handles=list(reversed(handlelist)), loc='upper left');
    
    plt.text(0.035, 0.09, '(a)',  transform = axes[0].transAxes,   fontsize=11)
    plt.text(0.035, 0.09, '(b)',  transform = axes[1].transAxes,   fontsize=11)
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Plot Isn curves at selected Vg values

def Plot_Isn( Global, IMGfilename, DataArray, VgValues, XRange, YRange ):
    fig = plt.figure(figsize=(4,2.5));
    plt.xlim(XRange);
    plt.ylim(YRange);
    plt.xlabel( '$eV_\mathrm{d}/\Delta$',               fontsize=11);
    plt.ylabel( '$eI_\mathrm{d}R_\mathrm{N}/\Delta$',   fontsize=11);
    plt.axvline(x = 0, color = 'black', linestyle = '--');
    plt.axvline(x = -1, color = 'grey', linestyle = '--');
    plt.axvline(x = 1, color = 'grey', linestyle = '--');
    for VgIndex in np.arange( 0, len(VgValues) ):
        # Select one Gdiff(vd) curve at Vg=VgValue from the original data array
        # Column 5 is Vg
        DataArraySliceAtVg          = Select_At_Vg(DataArray,VgValues[VgIndex])
        # Process the data: remove the filters
        DataArraySliceAtVgAdjusted  = groupby_mean( remove_filters1D( Global, DataArraySliceAtVg ) )
        # Column 1 is Id, column 2 is GdiffDUT
        RN                          = Extract_RNdiff( DataArraySliceAtVgAdjusted, Global['DeltaPtSi'] , 0.1 )
        # Factor 1000: current is in muA, Delta in mV
        SelectedCurveNormalized     = DataArraySliceAtVgAdjusted[:,1] * RN / ( 1000 * Global['DeltaPtSi'] )
        # Apply Savitzky-Golay filter, and center the curves
        FilteredCurve               = Center_Curve( savgol_filter( ( DataArraySliceAtVgAdjusted[:,0] ,
                                                       SelectedCurveNormalized ) , 3, 1) )
        MirrorCurve                 = Mirror_Curve( FilteredCurve )
        SortedAvgCurve              = np.transpose( groupby_mean( np.transpose( MirrorCurve ) ) )
        # Plot the Vd values vs the Gdiff values
        plt.plot( MirrorCurve[0] / Global['DeltaPtSi'] , MirrorCurve[1], label='{:3.2f} V'.format(VgValues[VgIndex]))
    plt.plot(XRange,XRange)
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1],loc='upper left')
    plt.savefig(IMGfilename,dpi=GlobalDPI,bbox_inches='tight');
    plt.close();

#%% Extract RNdiff and RN, plot both vs Vg

def Plot_RNdiffRN( Global, IMGfilename, DataArray, XRange, YRange ):
    VgValues        = Select_Data_XRange(DataArray , XRange)[:,0,5]
    RdiffArray      = np.array([0.]*len(VgValues))
    RNArray         = np.array([0.]*len(VgValues))
    for VgIndex in np.arange(0,len(VgValues)):
        # Select one Gdiff(vd) curve at Vg=VgValue from the original data array
        # Column 5 is Vg
        DataArraySliceAtVg          = DataArray[VgIndex,:,:]
        # Process the data: remove the filters
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )
        # Extract RN from Gdiff
        RdiffArray[VgIndex]         = Extract_RNdiff( DataArraySliceAtVgAdjusted, DeltaPtSi, 0.1 )
        # Extract RN from I, V
        RNArray[VgIndex]            = Extract_RN( DataArraySliceAtVgAdjusted, DeltaPtSi )
    fig             = plt.figure(figsize=(4,2.5));
    plt.xlabel('$V_\mathrm{g} (V)$',fontsize=11);
    plt.ylabel('$R_\mathrm{diff}$, $R_\mathrm{N}$',fontsize=11);
    plt.ylim(YRange);
    plt.plot( VgValues, RdiffArray, label='$R_\mathrm{diff}(eV=\Delta) (\Omega)$' );
    plt.plot( VgValues, RNArray,    label='$R_\mathrm{N}(eV=\Delta) (\Omega)$' );
    handles, labels = plt.gca().get_legend_handles_labels();
    plt.legend(handles[::-1], labels[::-1],loc='upper left');
    plt.savefig(    IMGfilename,    dpi=GlobalDPI,    bbox_inches='tight');
    plt.close();

def LogPlot_RNdiffRN( Global, IMGfilename, DataArray, XRange, YRange ):
    VgValues        = Select_Data_XRange(DataArray , XRange)[:,0,5]
    RdiffArray      = np.array([0.]*len(VgValues))
    RNArray         = np.array([0.]*len(VgValues))
    for VgIndex in np.arange(0,len(VgValues)):
        # Select one Gdiff(vd) curve at Vg=VgValue from the original data array
        # Column 5 is Vg
        DataArraySliceAtVg          = DataArray[VgIndex,:,:]
        # Process the data: remove the filters
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )
        # Extract RN from Gdiff
        RdiffArray[VgIndex]         = Extract_RNdiff( DataArraySliceAtVgAdjusted, Global['DeltaPtSi'], 0.1 )
        # Extract RN from I, V
        RNArray[VgIndex]            = Extract_RN( DataArraySliceAtVgAdjusted, Global['DeltaPtSi'] )
    fig             = plt.figure(figsize=(4,2.5));
    plt.xlabel('$V_\mathrm{g} (V)$',fontsize=11);
    plt.ylabel('$R$ $(\Omega)$',fontsize=11);
    plt.yscale('log');
    plt.ylim(YRange);
    plt.plot( VgValues, RdiffArray, label='$R_\mathrm{N,diff}(eV=\Delta)$' );
    plt.plot( VgValues, RNArray,    label='$R_\mathrm{N}(eV=\Delta)$' );
    handles, labels = plt.gca().get_legend_handles_labels();
    plt.legend(handles[::-1], labels[::-1],loc='upper left');
    plt.savefig( IMGfilename, dpi=GlobalDPI, bbox_inches='tight');
    plt.close();

def Plot_RNdiff_RN_Ratio( Global, IMGfilename, DataArray, XRange, YRange ):
    VgValues        = Select_Data_XRange(DataArray , XRange)[:,0,5]
    RdiffArray      = np.array([0.]*len(VgValues))
    RNArray         = np.array([0.]*len(VgValues))
    RatioArray      = np.array([0.]*len(VgValues))
    for VgIndex in np.arange(0,len(VgValues)):
        # Select one Gdiff(vd) curve at Vg=VgValue from the original data array
        # Column 5 is Vg
        DataArraySliceAtVg          = DataArray[VgIndex,:,:]
        # Process the data: remove the filters
        DataArraySliceAtVgAdjusted  = remove_filters1D( Global, DataArraySliceAtVg )
        # Extract RN from Gdiff
        RdiffArray[VgIndex]         = Extract_RNdiff( DataArraySliceAtVgAdjusted, Global['DeltaPtSi'], 0.1 )
        # Extract RN from I, V
        RNArray[VgIndex]            = Extract_RN( DataArraySliceAtVgAdjusted, Global['DeltaPtSi'] )
        if RdiffArray[VgIndex]>0:
            RatioArray[VgIndex]     = RNArray[VgIndex] / RdiffArray[VgIndex]
        else:
            RatioArray[VgIndex]     = 0
    fig             = plt.figure(figsize=(4,2.5));
    plt.xlabel('$V_\mathrm{g} (V)$',fontsize=11);
    plt.ylabel('$R_\mathrm{N}$/$R_\mathrm{N,diff}$',fontsize=11);
    plt.ylim(YRange);
    plt.plot( VgValues, RatioArray, label='$R_\mathrm{N}/R_\mathrm{N,diff}$ at $eV=\Delta$' );
    handles, labels = plt.gca().get_legend_handles_labels();
    plt.legend(handles[::-1], labels[::-1],loc='upper left');
    plt.savefig( IMGfilename, dpi=GlobalDPI, bbox_inches='tight');
    plt.close();
