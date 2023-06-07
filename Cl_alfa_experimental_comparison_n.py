import math
import numpy as np
import pandas as pd
from AeroPy_submodule.aeropy import xfoil_module
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import functions


fig,(ax1)=plt.subplots(1,1)
cl_values={}

Re=53000
n_ls=[7,8,9,10,11]
cl_values[Re]={}

for n in n_ls :
    N=151
    alfa_ls = np.linspace(0,15,N)
    alfa_values = []
    cl_values[Re][n] = []
    
    file_path = 'polar_Cl_alfa_xfoil_n7-8-9-10-11_exp/Polar_naca0012_0_15_05_{}_N{}.txt'.format(Re,n) 
    with open(file_path, 'r') as file:
        lines = file.readlines()
    data_start_index = 12  #Row index where reading has to start

    for line in lines[data_start_index:]:
        values = line.split()
        alpha = float(values[0])
        cl = float(values[1])
        
        alfa_values.append(alpha)
        cl_values[Re][n].append(cl)

    #Plot 
    ax1.plot(alfa_values,cl_values[Re][n],label='n={}'.format(n))
    ax1.legend()
    ax1.set_title("Cl-alfa Re=53000")
    ax1.set_ylabel("Cl")
    ax1.set_xlabel("alpha (°)")


# #Curve 2pi
# alfa_values=np.array(alfa_values)
# ax1.plot(alfa_values,(2*math.pi*alfa_values)/180*math.pi,label='2pi',color='black')
# ax1.legend()

#Experimental curve
df1 = pd.read_csv(r"C:\Users\matti\Desktop\BEMT\polar_Cl_alfa_xfoil_n7-8-9-10-11_exp\Cl_Re53000_Valerie.csv")
ax1.plot(df1.iloc[:,0],df1.iloc[:,1],label='Experimental Data',linestyle='dashed')
ax1.legend()

plt.show()

################################################################################################################################################
#REYNOLDS NUMBER = 140000
fig,(ax1)=plt.subplots(1,1)
cl_values={}

Re=140000
n_ls=[7,8,9,10,11]
cl_values[Re]={}

for n in n_ls :
    N=151
    alfa_ls = np.linspace(0,15,N)
    alfa_values = []
    cl_values[Re][n] = []
    
    file_path = 'polar_Cl_alfa_xfoil_n7-8-9-10-11_exp/Polar_naca0012_0_15_05_{}_N{}.txt'.format(Re,n) 
    with open(file_path, 'r') as file:
        lines = file.readlines()
    data_start_index = 12  #Row index where reading has to start

    for line in lines[data_start_index:]:
        values = line.split()
        alpha = float(values[0])
        cl = float(values[1])
        
        alfa_values.append(alpha)
        cl_values[Re][n].append(cl)

    #Plot 
    ax1.plot(alfa_values,cl_values[Re][n],label='n={}'.format(n))
    ax1.legend()
    ax1.set_title("Cl-alfa Re=140000")
    ax1.set_ylabel("Cl")
    ax1.set_xlabel("alpha (°)")

# #Curve 2pi
# alfa_values=np.array(alfa_values)
# ax1.plot(alfa_values,(2*math.pi*alfa_values)/180*math.pi,label='2pi',color='black')
# ax1.legend()

#Experimental curve
df2 = pd.read_csv(r"C:\Users\matti\Desktop\BEMT\polar_Cl_alfa_xfoil_n7-8-9-10-11_exp\Cl_Re140000_Valerie.csv")
ax1.plot(df2.iloc[:,0],df2.iloc[:,1],label='Experimental Data',linestyle='dashed')
ax1.legend()

plt.show()



# 
# ax1.plot(df2.iloc[:,0],df2.iloc[:,1],label='Re=140000 Experimental Data',color="red",linestyle='dashed')
# ax1.legend(loc='center left', bbox_to_anchor=(0, 0.65))

# #Difference EXP - XFOIL
# f1=interp1d(df1.iloc[:,0],df1.iloc[:,1], kind='linear',fill_value="extrapolate")
# f2=interp1d(df2.iloc[:,0],df2.iloc[:,1], kind='linear',fill_value="extrapolate")

# ax2.plot(alfa_values, f1(alfa_values)-cl_values[53000],label="(EXP - XFOIL) Re=53000",color='blue')
# ax2.axhline(0, color='black',linewidth=0.8)
# ax2.legend()
# ax2.set_ylabel("delta Cl")
# ax2.set_ylim(-0.3,0.3)


# ax3.plot(alfa_values, f2(alfa_values)-cl_values[140000],label="(EXP - XFOIL) Re=140000",color='red')
# ax3.axhline(0, color='black',linewidth=0.8)
# ax3.legend()
# ax3.set_ylabel("delta Cl")
# ax3.set_xlabel("alfa (°)")
# ax3.set_ylim(-0.3,0.3)




#Show plot
plt.show()