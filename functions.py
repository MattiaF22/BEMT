import math

#FIND a
def solve_a(Solidity,lambda1,inflow_angle):
    k=Solidity*lambda1/(4*(math.sin(inflow_angle))**2)
    a=k/(1-k)
    return a

#FIND a_prime
def solve_a_prime(Solidity,lambda2,inflow_angle):
    k=Solidity*lambda2/(2*math.sin(2*inflow_angle))
    a_prime=k/(1+k)
    return a_prime

#FIND dCT_dr_ad
def solve_dCT_dr_ad(Solidity,lambda1,r_ad,a_prime,inflow_angle):
    dCT_dr_ad = (math.pi**3)*Solidity*lambda1*(r_ad**3)*((1-a_prime)**2)/(4*((math.cos(inflow_angle))**2))
    return dCT_dr_ad
    
#FIND dCP_dr_ad
def solve_dCP_dr_ad(Solidity,lambda2,r_ad,a_prime,inflow_angle):
    dCP_dr_ad = (math.pi**4)*Solidity*lambda2*(r_ad**4)*((1-a_prime)**2)/(4*((math.cos(inflow_angle))**2))
    return dCP_dr_ad

#FIND J_calculated
def solve_J_calculated(r_ad,a_prime,a,inflow_angle): 
    J_calculated = math.pi*r_ad*(1-a_prime)*(math.tan(inflow_angle))/(1+a)
    return J_calculated

def nearest_value(reference_value, list_of_values):
    near_value = min(list_of_values, key=lambda x: abs(x - reference_value))
    index = list_of_values.index(near_value)
    return near_value,index



