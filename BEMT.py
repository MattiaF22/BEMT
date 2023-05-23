import math
import numpy as np
import pandas as pd
import aeropy.xfoil_module
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import functions

#########################################################################################################################

#Input data
R=0.125            #radius
R_hub=0.016        #hub radius
c=0.025            #chord
Theta=25           #blade pitch angle
V_inf=5           #rate of climb
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
    N=21
    alfa_ls = np.linspace(0,15,N)
    alfa_values = []
    cl_values=[] 
    cd_values=[]
    #Create Polar file for each Reynolds number
    aeropy.xfoil_module.call('naca0012', alfas=alfa_ls, output='Polar', Reynolds=Re, Mach=0, 
                         plots=False, NACA=True, GDES=False,  iteration=50, flap=None, PANE=False, NORM=True) 
    
    
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


df = pd.DataFrame(columns=['r_ad','Reynolds', 'alfa', 'a', 'a_prime','dCT_dr_ad', 'dCP_dr_ad' , 'J_calculated' ])

iterations=20
a_new=np.zeros(M)
a_prime_new=np.zeros(M)
for k in range(iterations):
    a=a_new
    a_prime=a_prime_new
    alfa_seq=[]
    a_seq=[]
    
    for i in range(M) : #for each r_ad there is a Reynolds number
        r_ad=r_ad_ls[i]
        Re = int(Re_ls[i])
        Solidity = B*c/(2*math.pi*r_ad*R)
        inflow_angle= math.atan((V_inf*(1+a[i]))/(Omega*r_ad*R*(1-a_prime[i])))*180/math.pi
        alfa=Theta-inflow_angle
        alfa_seq.append(alfa)
        Cl = Interp_functions_Cl[Re](alfa)
        Cd = Interp_functions_Cd[Re](alfa)
        inflow_angle = math.radians(Theta - alfa)
        lambda1 = Cl*math.cos(inflow_angle)-Cd*math.sin(inflow_angle)
        lambda2 = Cl*math.sin(inflow_angle)+Cd*math.cos(inflow_angle)
        z=Solidity*lambda1/(4*(math.sin(inflow_angle))**2)
        a_c = z/(1-z)
        a_seq.append(a_c)
        a_new[i]=a_c
        q=Solidity*lambda2/(2*math.sin(2*inflow_angle))
        a_prime_c = q/(1+q)
        a_prime_new[i]=a_prime_c
        #dCT_dr_ad = float(functions.solve_dCT_dr_ad(Solidity,lambda1,r_ad,a_prime_c,inflow_angle))
        #dCP_dr_ad = float(functions.solve_dCP_dr_ad(Solidity,lambda2,r_ad,a_prime_c,inflow_angle))
        #J_calculated = float(functions.solve_J_calculated(r_ad,a_prime_c,a_c,inflow_angle))
    print(alfa_seq)
    print(a_seq)
  
        #df.loc[i] = [r_ad,Re, alfa, a, a_prime,dCT_dr_ad, dCP_dr_ad ,J_calculated ]
    #print(df)


