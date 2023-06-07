import math
import numpy as np
import pandas as pd
from AeroPy_submodule.aeropy import xfoil_module
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import functions


#Input data
R=0.125            #radius
R_hub=0.016        #hub radius
c=0.025            #chord
Theta=10           #blade pitch angle
V_inf=1.77         #rate of climb
n=100              #frequency
B=2                #blades
rho=1.225          #density
mu=1.789e-5        #dynamic viscosity
Omega=n*2*math.pi  #tip velocity
D=2*R              #diameter
J=V_inf/ (n*D)     #advance ratio

####################################################################################################################################

M=11  #number of stations along the span
r_ad_ls = np.linspace(R_hub/R, 1,M)  #adimensional local radius
Re_ls=rho*(Omega*r_ad_ls*R)*c / mu
Re_ls=list(map(int, Re_ls))


fig, (ax1,ax2)=plt.subplots(1,2)
legend_labels = []

for i in range(0,len(Re_ls)) :
    Re=Re_ls[i]
    N=151
    alfa_ls = np.linspace(0,15,N)
    alfa_values = []
    cl_values=[] 
    cd_values=[]
    #Create Polar file for each Reynolds number
    #xfoil_module.call('naca0012', alfas=alfa_ls, output='Polar', Reynolds=Re, Mach=0, 
    #                     plots=False, NACA=True, GDES=False,  iteration=100, flap=None, PANE=False, NORM=True) 
    
    
    Reynolds_number=str(Re)

    file_path = 'polar_Cl_alfa_xfoil_n9/Polar_naca0012_{}_0.0_15.0'.format(Reynolds_number) 

    with open(file_path, 'r') as file:
        lines = file.readlines()

    data_start_index = 12  #Row index where reading has to start

    for line in lines[data_start_index:]:
        values = line.split()
        alpha = float(values[0])
        cl = float(values[1])
        cd = float(values[2])
        
        alfa_values.append(alpha)
        cl_values.append(cl)
        cd_values.append(cd)

    #Plot 
    legend_labels.append("Re = {}".format(Re)) 

    ax1.plot(alfa_values,cl_values,label='{}'.format(Re))
    ax1.set_title("Cl-alfa")
    ax1.set_xlabel("alfa") 
    ax1.set_ylabel("Cl")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax2.plot(alfa_values,cd_values)
    ax2.set_title("Cd-alfa")
    ax2.set_xlabel("alfa") 
    ax2.set_ylabel("Cd")

    #Curve 2pi
    alfa_values=np.array(alfa_values)
    ax1.plot(alfa_values,(2*math.pi*alfa_values)/180*math.pi,color="red")

    plt.show()     #INSERT BREAKPOINT



fig, (ax1,ax2)=plt.subplots(1,2)

Re_ls=[134448, 2e5, 5e5, 1e6, 5e6, 1e7]
for i in range(0,len(Re_ls)) :
    Re=Re_ls[i]
    N=151
    alfa_ls = np.linspace(0,15,N)
    alfa_values = []
    cl_values=[] 
    cd_values=[]
    #Create Polar file for each Reynolds number
    #xfoil_module.call('naca0012', alfas=alfa_ls, output='Polar', Reynolds=Re, Mach=0, 
    #                     plots=False, NACA=True, GDES=False,  iteration=100, flap=None, PANE=False, NORM=True) 
    
    
    Reynolds_number=str(Re)

    file_path = 'polar_Cl_alfa_xfoil_n9/Polar_naca0012_{}_0.0_15.0'.format(Reynolds_number) 

    with open(file_path, 'r') as file:
        lines = file.readlines()

    data_start_index = 12  #Row index where reading has to start

    for line in lines[data_start_index:]:
        values = line.split()
        alpha = float(values[0])
        cl = float(values[1])
        cd = float(values[2])

        # #REMEMBER THAT THESE VALUES ARE FOR CHORD = 1 !!! THEN SOLVE THIS !!! PAY ATTENTION TO CD THAT IS NOT PORPORZIONAL TO CHORD
        
        alfa_values.append(alpha)
        cl_values.append(cl)
        cd_values.append(cd)

    #Plot 
    legend_labels.append("Re = {}".format(Re)) 

    ax1.plot(alfa_values,cl_values,label='{}'.format(Re))
    ax1.set_title("Cl-alfa")
    ax1.set_xlabel("alfa") 
    ax1.set_ylabel("Cl")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax2.plot(alfa_values,cd_values)
    ax2.set_title("Cd-alfa")
    ax2.set_xlabel("alfa") 
    ax2.set_ylabel("Cd")

    #Curve 2pi
    alfa_values=np.array(alfa_values)
    ax1.plot(alfa_values,(2*math.pi*alfa_values)/180*math.pi,color="red")

    
    plt.show()     #INSERT BREAKPOINT



