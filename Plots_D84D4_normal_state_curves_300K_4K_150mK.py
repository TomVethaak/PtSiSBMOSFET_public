# -*- coding: utf-8 -*-
from Standardized_plots import *
from DataDirectoryLibrary import *

User            = "Tom"
#User           = "Francois"

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

#%% Initialize 
RootDirectory   = GetRootDirectory( User )

#%% IMPORT

# Single file:
DataArray300K   = np.array( np.loadtxt( RootDirectory + DataFileList[0] ) )
DataArray4K     = np.array( np.loadtxt( RootDirectory + DataFileList[1] ) )
DataArray150mK  = np.array( np.loadtxt( RootDirectory + DataFileList[2] ) )

#%% PLOT 1D CURVES

plt.ioff()

# Number 1
fig     = plt.figure(figsize=(3,2.8))
plt.plot( groupby_mean( DataArray300K )[:,0], 1000* DataArray300K[0,1] - 1000* savgol_filter( groupby_mean( DataArray300K )[:,1], 3, 1 ) ,
         label='300 K')
plt.xlabel('$V_\mathrm{g}$ (V)',fontsize=13);
plt.ylabel('$I_d$ (nA)  at  $V_\mathrm{d}=1$ mV',fontsize=13);
plt.gca().set_xlim([-4,0])
plt.legend()
plt.tight_layout()
plt.savefig('img/D84D4_Id_Vg_300K_4K_150mK_1.png',dpi=500,bbox_inches='tight')
plt.close()

# Number 2
fig     = plt.figure(figsize=(3,2.8))
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
plt.close()

# Number 3
fig     = plt.figure(figsize=(3,2.8))
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
plt.close()