import math
import numpy as np
import pandas as pd
from AeroPy_submodule.aeropy import xfoil_module
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from scipy.integrate import simps
from matplotlib import pyplot as plt
import functions

#########################################################################################################################

#Input data
R=0.125            #radius
R_hub=0.016        #hub radius
c=0.025            #chord
Theta=10           #blade pitch angle
V_inf=1.77         #rate of climb
RPM=6000           #rotations per minute
n=RPM/60           #frequency
B=2                #blades
rho=1.225          #density
mu=1.789e-5        #dynamic viscosity
Omega=n*2*math.pi  #tip velocity
D=2*R              #diameter
J=V_inf/ (n*D)     #advance ratio
#alfa_desired=8

####################################################################################################################################


#Calculate Reynolds number along the span
M=11  #number of stations along the span
r_ad_ls = np.linspace(R_hub/R, 1,M)  #adimensional local radius
Re_ls=rho*(Omega*r_ad_ls*R)*c / mu
Re_ls=list(map(int, Re_ls))

Interp_functions_Cl = {}
Interp_functions_Cd = {}

for Re in Re_ls :
    N=101
    alfa_ls = np.linspace(0,10,N)
    alfa_values = []
    cl_values=[] 
    cd_values=[]
    #Create Polar file for each Reynolds number
    #aeropy.xfoil_module.call('naca0012', alfas=alfa_ls, output='Polar', Reynolds=Re, Mach=0, 
    #                     plots=False, NACA=True, GDES=False,  iteration=100, flap=None, PANE=False, NORM=True) 
    
    
    #Take the data from these files and save them in a dictionary
    Reynolds_number=str(Re)

    file_path = 'Polar_naca0012_{}_0.0_15.0'.format(Reynolds_number) 

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
    f =  interp1d(alfa_values, cl_values, kind='linear',fill_value="extrapolate")
    Interp_functions_Cl[Re] = f
    g = interp1d(alfa_values, cd_values, kind='linear',fill_value="extrapolate")
    Interp_functions_Cd[Re] = g
#Now I have two functions that give me Cl and Cd for every alfa I want, also not present in xfoil tables
print(Interp_functions_Cd[Re_ls[5]](2.12))


#######################################################################################################################################################


df = pd.DataFrame(columns=['r_ad','Reynolds', 'alfa', 'J','a','a_prime','dCT_dr_ad','dCP_dr_ad'])

for i in range(M) : #for each r_ad there is a Reynolds number
    r_ad=r_ad_ls[i]
    Re = int(Re_ls[i])
    Solidity = B*c/(2*math.pi*r_ad*R)
    a_ls=[]
    a_prime_ls=[]
    dCT_dr_ad_ls=[]
    dCP_dr_ad_ls=[]
    J_ls=[]

    for j in range(N):
        alfa=alfa_ls[j]
        Cl = Interp_functions_Cl[Re](alfa)
        Cd = Interp_functions_Cd[Re](alfa)
        #Theta=alfa_desired+(math.atan(V_inf/(Omega*r_ad*R)))*180/math.pi     #iperbolic blade pitch angle distribution in order to have alfa=alfa_desired (without considering inductions)
        inflow_angle = math.radians(Theta - alfa)
        lambda1 = Cl*math.cos(inflow_angle)-Cd*math.sin(inflow_angle)
        lambda2 = Cl*math.sin(inflow_angle)+Cd*math.cos(inflow_angle)
        a = functions.solve_a(Solidity,lambda1,inflow_angle)
        a_ls.append(a)
        a_prime = functions.solve_a_prime(Solidity,lambda2,inflow_angle)
        a_prime_ls.append(a_prime)
        dCT_dr_ad = functions.solve_dCT_dr_ad(Solidity,lambda1,r_ad,a_prime,inflow_angle)
        dCT_dr_ad_ls.append(dCT_dr_ad)
        dCP_dr_ad = functions.solve_dCP_dr_ad(Solidity,lambda2,r_ad,a_prime,inflow_angle)
        dCP_dr_ad_ls.append(dCP_dr_ad)
        J_calculated = functions.solve_J_calculated(r_ad,a_prime,a,inflow_angle)
        J_ls.append(J_calculated)
    near_J,index = functions.nearest_value(J,J_ls)
    alfa=alfa_ls[index]    
    a=a_ls[index]
    a_prime=a_prime_ls[index]
    dCT_dr_ad=dCT_dr_ad_ls[index]
    dCP_dr_ad=dCP_dr_ad_ls[index]

    df.loc[i] = [r_ad, Re, alfa, near_J, a, a_prime, dCT_dr_ad, dCP_dr_ad]

    print(df)
dCT_dr_ad=df["dCT_dr_ad"].values
dCP_dr_ad=df["dCP_dr_ad"].values
r_ad=df["r_ad"].values
CT=simps(dCT_dr_ad,r_ad)
CP=simps(dCP_dr_ad,r_ad)
print(CT)
print(CP)



