import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate as intg
from lmfit import Model

###  Constants  ###

EnbPoints = 500 # E sampling for the integral
minE = -2e-3 #eV
maxE = 2e-3 #eV
E = np.linspace(minE,maxE,EnbPoints)
Epos = np.linspace(0,maxE,EnbPoints)

Vmin = -1.5e-3 #V  Vbias min
Vmax = 1.5e-3 #V   Vbias max
VnbPoints = 200 # Vbias sampling
Zlist = [1.9] # Barrier height
T = 0.5 #K Temperature
Delta = 130e-6 #eV Superconducting Gap
K = 8.6e-5 #eV/K  Boltzmann
G =68e-6 #eV  Dynes parameter


### BTK Functions  ###

def u(E, Delta, G):
    return (((1+ ( (complex(E,-G)**2-Delta**2)/E**2 )**0.5 )/2)**0.5)
    
def v(E, Delta, G):
    return ((1-u(E,Delta,G)**2)**0.5)
    
def gama(E, Z, Delta, G):
    return (u(E,Delta,G)**2 + (u(E,Delta,G)**2 - v(E,Delta,G)**2)*Z**2)

def A(E, Z, Delta, G):
    if abs(E) > Delta:
        return( (u(E,Delta,G)**2 * v(E,Delta,G)**2) / gama(E,Z,Delta,G)**2) 
    else:
        return( Delta**2 / (complex(E,-G)**2 + (Delta**2-complex(E,-G)**2)*(1+2*Z**2)**2))
        
def B(E, Z, Delta, G):
    if abs(E) > Delta:
        return ( ((u(E,Delta,G)**2 - v(E,Delta,G)**2)**2 * Z**2 * (1+Z**2)) / gama(E,Z,Delta,G)**2)
    else:
        return (1 - A(E,Z,Delta,G))
        
def C(E, Z, Delta, G):
    if abs(E) > Delta:
        return ( (u(E,Delta,G)**2 * (u(abs(E),Delta,G)**2 - v(E,Delta,G)**2) * (1+Z**2)) / gama(E,Z,Delta,G)**2)
    else:
        return (0)

def D(E, Z, Delta, G):
    if abs(E) > Delta:
        return ( (v(E,Delta,G)**2 * (u(E,Delta,G)**2 - v(E,Delta,G)**2) * Z**2) / gama(E,Z,Delta,G)**2)
    else:
        return (0)
        
def Fermi(En, T, V):
    return 1/(1 + np.exp( (En-V)/(K*T) ))
    

def Conductance(V, Z, Delta, T, G, Einterval):
    FermiValues = np.zeros(len(Einterval))
    for i in range(len(Einterval)):
        FermiValues[i] = Fermi(Einterval[i], T, V)
    
    integrandeValues = np.zeros(len(Einterval) - 1)
    for i in range(len(Einterval) - 1):
        integrandeValues[i] = (1 + A(Einterval[i],Z,Delta,G) - B(Einterval[i],Z,Delta,G)).real * (np.diff(FermiValues)/np.diff(Einterval))[i]
        
    return( -(1+Z**2) * intg.trapz(integrandeValues, Einterval[:-1]) ) #Trapeze methode is used here for integration
    
##############################################################################

def ConductanceFit(Vlist, Z, Delta, T, G):
    condList = []
    for i in range(len(Vlist)):
        condList.append(Conductance(Vlist[i], Z, Delta, T, G, E))
    return (condList)

### The ConductanceFit fonction is used to create a lmfit Model object and guess parameters are defined into the object parameters
gmodel = Model(ConductanceFit)
gmodel.set_param_hint('Z', value=1)
gmodel.set_param_hint('T', value=0.35, vary=False)
gmodel.set_param_hint('Delta', value=150e-6)
gmodel.set_param_hint('G', value=70e-6, vary=True)
parameters = gmodel.make_params()

print(f'parameter names: {gmodel.param_names}')
print(f'independent variables: {gmodel.independent_vars}')


### The experimental data are imported
data250 = np.genfromtxt ('BTK_data_vs_Vdg\Gdiff_Vd_g=-2.5V.dat', delimiter=' ')[3200:6200,:]
data275 = np.genfromtxt ('BTK_data_vs_Vdg\Gdiff_Vd_g=-2.75V.dat', delimiter=' ')[3200:6200,:]
data300 = np.genfromtxt ('BTK_data_vs_Vdg\Gdiff_Vd_g=-3.0V.dat', delimiter=' ')[3900:5600,:]
data325 = np.genfromtxt ('BTK_data_vs_Vdg\Gdiff_Vd_g=-3.25V.dat', delimiter=' ')[3800:5700,:]
data350 = np.genfromtxt ('BTK_data_vs_Vdg\Gdiff_Vd_g=-3.5V.dat', delimiter=' ')[4000:5500,:]

### Normalization parameters
mult250 = 66000
mult275 = 38000
mult300 = 28500
mult325 = 25500
mult350 = 24300

multlist = [mult250, mult275, mult300, mult325, mult350]
datalist = [data250, data275, data300, data325, data350]

### Variables to store the fitting parameters (the output of the algo)
Vg_list = [-2.5, -2.75, -3.0, -3.25, -3.5] 
param_result = np.zeros((len(datalist),5))

### Experimental sampling parameters
initial_sampling_lim = 5e-4
sampling_lim = initial_sampling_lim

for x in range(len(datalist)):
    data = datalist[x]
    mult = multlist[x]
    
    ### Normalization ###
    data_norm = [[],[]]
    data_norm[1] = data[:,1]*mult
    data_norm[0] = data[:,0]/1000 +0.063e-3
        
    ### Sampling (higher sampling far from the gap than inside)###
    datater = [[],[]]
    for i in range(len(data_norm[0])):
        if abs(data_norm[0][i]) < sampling_lim: 
            if i % 16 == 0:
                datater[0].append(data_norm[0][i])
                datater[1].append(data_norm[1][i])
        else :
            if i % 34 == 0:
                datater[0].append(data_norm[0][i])
                datater[1].append(data_norm[1][i])  
    
    ### Fit calcultaion ###
    result = gmodel.fit(datater[1],  parameters, Vlist=datater[0])
    
    ### Fit plot ###
    print(result.fit_report())
    result.plot()
    
    plt.plot(datater[0], datater[1],         'bo')
    plt.plot(datater[0], result.init_fit, 'k--')
    plt.plot(datater[0], result.best_fit, 'r-')
    plt.show()
    
    ### Fit parameter save ###
    param_result[x,0] = Vg_list[x]
    param_result[x,1] = result.best_values['Z']
    param_result[x,2] = result.best_values['Delta']
    param_result[x,3] = result.best_values['T']
    param_result[x,4] = result.best_values['G']
    
    ### Initial guess based on the previous fit ###
    #parameters = result.params
    sampling_lim = 2*result.best_values['Delta']
    
    
### Plot the fit parameters ###
fig, axs = plt.subplots(3, sharex=True)
fig.suptitle('Fit parameters')
axs[0].plot(Vg_list, param_result[:,1])
axs[0].set_ylabel('Z')
axs[1].plot(Vg_list, param_result[:,2])
axs[1].set_ylabel('Delta (eV)')
axs[2].plot(Vg_list, param_result[:,4])
axs[2].set_ylabel('G (eV)')
axs[2].set_xlabel('Vg (V)')
plt.show()




"""
new_result = gmodel.fit(datater[1], result.params, Vlist=datater[0])       
print(new_result.fit_report())

plt.plot(datater[0], datater[1],         'bo')
plt.plot(datater[0], new_result.init_fit, 'k--')
plt.plot(datater[0], new_result.best_fit, 'r-')
plt.show()            

"""






