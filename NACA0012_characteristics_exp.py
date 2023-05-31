import math
import numpy as np
import pandas as pd
from AeroPy_submodule.aeropy import xfoil_module
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import functions


Re_xfoil=[53000,140000]

fig,(ax1,ax2,ax3)=plt.subplots(3,1)
cl_values={}


for i in range(len(Re_xfoil)) :
    Re=Re_xfoil[i]
    N=151
    alfa_ls = np.linspace(0,15,N)
    alfa_values = []
    cl_values[Re] = []


    #Create Polar file for each Reynolds number
    #xfoil_module.call('naca0012', alfas=alfa_ls, output='Polar', Reynolds=Re, Mach=0, 
    #                    plots=False, NACA=True, GDES=False,  iteration=100, flap=None, PANE=False, NORM=True) 
    
    
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
        
        alfa_values.append(alpha)
        cl_values[Re].append(cl)

    #Plot 
    color=['blue','red']
    ax1.plot(alfa_values,cl_values[Re],label='Re = {} XFOIL'.format(Re),color=color[i])
    ax1.set_title("Cl-alfa")
    ax1.set_ylabel("Cl")


#Curve 2pi
alfa_values=np.array(alfa_values)
ax1.plot(alfa_values,(2*math.pi*alfa_values)/180*math.pi,label='2pi',color='black')

#Experimental curves
df1 = pd.read_csv(r"C:\Users\matti\Desktop\BEMT\Cl_Re53000_Valerie.csv")
df2 = pd.read_csv(r"C:\Users\matti\Desktop\BEMT\Cl_Re140000_Valerie.csv")

ax1.plot(df1.iloc[:,0],df1.iloc[:,1],label='Re=53000 Experimental Data',color="blue",linestyle='dashed')
ax1.plot(df2.iloc[:,0],df2.iloc[:,1],label='Re=140000 Experimental Data',color="red",linestyle='dashed')
ax1.legend(loc='center left', bbox_to_anchor=(0.95, 0.5))

#Difference EXP - XFOIL
f1=interp1d(df1.iloc[:,0],df1.iloc[:,1], kind='linear',fill_value="extrapolate")
f2=interp1d(df2.iloc[:,0],df2.iloc[:,1], kind='linear',fill_value="extrapolate")

ax2.plot(alfa_values, f1(alfa_values)-cl_values[53000],label="(EXP - XFOIL) Re=53000")
ax2.axhline(0, color='black')
ax2.legend()
ax2.set_ylabel("delta Cl")


ax3.plot(alfa_values, f2(alfa_values)-cl_values[140000],label="(EXP - XFOIL) Re=140000")
ax3.axhline(0, color='black')
ax3.legend()
ax3.set_ylabel("delta Cl")
ax3.set_xlabel("alfa °")



#Show plot
plt.show()