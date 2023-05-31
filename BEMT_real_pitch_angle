import math
import numpy as np
import pandas as pd
from AeroPy_submodule.aeropy import xfoil_module
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import functions

#########################################################################################################################

#Input data
R=0.125            #radius
R_hub=0.016        #hub radius
c=0.025            #chord
V_inf=30           #rate of climb
alfa_desired=8     #angle of attach at each station    
n=100              #frequency
B=2                #blades
rho=1.225          #density
mu=1.789e-5        #dynamic viscosity
Omega=n*2*math.pi  #tip velocity
D=2*R              #diameter
J=V_inf/ (n*D)     #advance ratio

####################################################################################################################################


#Calculate Reynolds number along the span
M=11  #number of stations along the span
r_ad_ls = np.linspace(R_hub/R, 1,M)  #adimensional local radius
Re_ls=rho*(Omega*r_ad_ls*R)*c / mu
Re_ls=list(map(int, Re_ls))

Interp_functions_Cl = {}
Interp_functions_Cd = {}

for Re in Re_ls :
    N=151
    alfa_ls = np.linspace(0,15,N)
    alfa_values = []
    cl_values=[] 
    cd_values=[]
    #Create Polar file for each Reynolds number
    #xfoil_module.call('naca0012', alfas=alfa_ls, output='Polar', Reynolds=Re, Mach=0, 
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



#######################################################################################################################################################


df = pd.DataFrame(columns=['N_it','r_ad','Reynolds', 'alfa', 'inflow', 'Theta' ,'a', 'a_prime','dCT_dr_ad', 'dCP_dr_ad' , 'J_calculated', 'real J'])

iterations=20
a=np.zeros(M)
a_prime=np.zeros(M)
for k in range(iterations):
    
    for i in range(M) : #for each r_ad there is a Reynolds number
        r_ad=r_ad_ls[i]
        Re = Re_ls[i]
        Theta=alfa_desired+(math.atan(V_inf/(Omega*r_ad*R)))*180/math.pi     #iperbolic blade pitch angle distribution in order to have alfa=alfa_desired (without considering inductions)
        Solidity = B*c/(2*math.pi*r_ad*R)
        inflow_angle= math.atan((V_inf*(1+a[i]))/(Omega*r_ad*R*(1-a_prime[i])))*180/math.pi
        alfa=Theta-inflow_angle
        Cl = Interp_functions_Cl[Re](alfa)
        Cd = Interp_functions_Cd[Re](alfa)
        inflow_angle = math.radians(inflow_angle)
        lambda1 = Cl*math.cos(inflow_angle)-Cd*math.sin(inflow_angle)
        lambda2 = Cl*math.sin(inflow_angle)+Cd*math.cos(inflow_angle)
        a[i] = functions.solve_a(Solidity,lambda1,inflow_angle)
        a_prime[i] = functions.solve_a_prime(Solidity,lambda2,inflow_angle)
        dCT_dr_ad = functions.solve_dCT_dr_ad(Solidity,lambda1,r_ad,a_prime[i],inflow_angle)
        dCP_dr_ad = functions.solve_dCP_dr_ad(Solidity,lambda2,r_ad,a_prime[i],inflow_angle)
        J_calculated = functions.solve_J_calculated(r_ad,a_prime[i],a[i],inflow_angle)
        alfa="{:.2f}".format(alfa)
        Theta="{:.2f}".format(Theta)
        inflow_angle="{:.2f}".format(math.degrees(inflow_angle))  
        df.loc[i] = [k+1,r_ad,Re, alfa, inflow_angle,Theta, a[i], a_prime[i],dCT_dr_ad, dCP_dr_ad ,J_calculated, J]
    print(df)