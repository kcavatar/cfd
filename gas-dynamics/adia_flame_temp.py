## Calculation of Adiabatic Flame Temperature for LPG-air mixture
#         Using General expression for specific heat, Cp
#           G R Krishna Chand Avatar, SR : 16558

import numpy as np

##***********************************************************************************##

def enthalpy(a, T, To):

    # Module to evaluate enthalpy 
    # Inputs : a = coefficients, T = temperature
    #          To = reference temp
    Ru = 8.314  # J/K/mol
    part1 = T*(a[1] + a[2]*T/2 + a[3]*T**2/3 + a[4]*T**3/4 + a[5]*T**4/5)
    part2 = To*(a[1] + a[2]*To/2 + a[3]*To**2/3 + a[4]*To**3/4 + a[5]*To**4/5)
    h = Ru * (part1 - part2)
    return h

def func_evaluate(T, phi):
    # Cp = Ru(a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4)
    # hof = heat of formation
    # h = enthalpy

    # C_p curve fitting coefficients (Obtained from Steven Turns book)
    a1 = np.zeros(6)
    a2 = np.zeros(6)
    
    # Reference temperature
    To = 298  # K

    # CO_2
    # Temp range: 298 - 1000 K
    a1[1] = 0.02275724 * 10**2; a1[2] = 0.09922072 * 10**-1; a1[3] = -0.10409113 * 10**-4;
    a1[4] = 0.06866686 * 10**-7; a1[5] = -0.02117280 * 10**-10
    
    # Temp range: 1000-5000K
    a2[1] = 0.04453623 * 10**2; a2[2] = 0.03140168 * 10**-1; a2[3] = -0.12784105 * 10**-5;
    a2[4] = 0.02393996 * 10**-8; a2[5] = -0.16690333 * 10**-13
    
    hof_co2 = -393546  # J/mol
    
    if (T<=1000):
       h_co2 = enthalpy(a1, T, To)
    else:
       h_co2 = enthalpy(a1, 1000, To) + enthalpy(a2, T, 1000) 
        
    # H20
    # Temp range: 298 - 1000 K
    a1[1] = 0.03386842 * 10**2; a1[2] = 0.03474982 * 10**-1; a1[3] = -0.06354696 * 10**-4  ;
    a1[4] = 0.06968581 * 10**-7; a1[5] = -0.02506588 * 10**-10
    
    # Temp range: 1000-5000K
    a2[1] = 0.02672145 * 10**2; a2[2] = 0.03056293 * 10**-1; a2[3] = -0.08730260 * 10**-5;
    a2[4] = 0.12009964 * 10**-9; a2[5] = -0.06391618 * 10**-13

    hof_h2o = -241845
    if (T<=1000):
       h_h2o = enthalpy(a1, T, To)
    else:
       h_h2o = enthalpy(a1, 1000, To) + enthalpy(a2, T, 1000) 

    # O_2
    # Temp range: 298 - 1000 K
    a1[1] = 0.03212936 * 10**2; a1[2] = 0.11274864 * 10**-2; a1[3] = -0.05756150 * 10**-5;
    a1[4] = 0.13138773 * 10**-8; a1[5] = -0.08768554 * 10**-11
    
    # Temp range: 1000-5000K
    a2[1] = 0.03697578 * 10**2; a2[2] = 0.06135197 * 10**-2; a2[3] = -0.12588420 * 10**-6;
    a2[4] = 0.01775281 * 10**-9; a2[5] = -0.11364354 * 10**-14
    hof_o2 = 0
    if (T<=1000):
       h_o2 = enthalpy(a1, T, To)
    else:
       h_o2 = enthalpy(a1, 1000, To) + enthalpy(a2, T, 1000) 

    # N_2
    # Temp range: 298 - 1000 K
    a1[1] = 0.03298677 * 10**2; a1[2] = 0.14082404 * 10**-2; a1[3] = -0.03963222 * 10**-4;
    a1[4] = 0.05641515 * 10**-7; a1[5] = -0.02444854 * 10**-10
    
    # Temp range: 1000-5000K
    a2[1] = 0.02926640 * 10**2; a2[2] = 0.14879768 * 10**-2; a2[3] = -0.05684760 * 10**-5; 
    a2[4] = 0.10097038 * 10**-9; a2[5] = -0.06753351 * 10**-13
    hof_n2 = 0
    if (T<=1000):
       h_n2 = enthalpy(a1, T, To)
    else:
       h_n2 = enthalpy(a1, 1000, To) + enthalpy(a2, T, 1000) 

    # CO
    # Temp range: 298 - 1000 K
    a1[1] = 0.03262451 * 10**2; a1[2] = 0.15119409 * 10**-2; a1[3] = -0.03881755 * 10**-4;
    a1[4] = 0.05581944 * 10**-7; a1[5] = -0.02474951 * 10**-10
    
    # Temp range: 1000-5000K
    a2[1] = 0.03025078 * 10**2; a2[2] = 0.14426885 * 10**-1; a2[3] = -0.05630827 * 10**-5;
    a2[4] = 0.10185813 * 10**-9; a2[5] = -0.06910951 * 10**-13

    hof_co = -110541
    if(phi > 1.0):
       if (T<=1000):
          h_co = enthalpy(a1, T, To)
       else:
          h_co = enthalpy(a1, 1000, To) + enthalpy(a2, T, 1000) 

    # C4H10 (Butane)
    hof_but = -124733 # J/mol

    # C3H8 (Propane)
    hof_prop = -103847 # J/mol

    ## Reaction computations

    # Moles of the reacting substances
    n_but = phi * 0.6; n_prop = phi * 0.4; n_o2 = 5.9; n_n2 = 3.76 * 5.9

    # Heat of reaction
    H_reac = n_but * hof_but + n_prop * hof_prop + n_o2 * hof_o2 + n_n2 * hof_n2

    if (phi <= 1.0): # Fuel lean

           n_co2 = 3.6 * phi; n_h2o=4.6 * phi; n_o2 = 5.9*(1-phi) # O2 after reaction
           hof_prod = n_co2 * hof_co2 + n_h2o * hof_h2o + n_n2 * hof_n2
           h_prod = n_co2 * h_co2 + n_h2o * h_h2o + n_o2*h_o2 +  n_n2 * h_n2
           H_prod = hof_prod + h_prod

    else:
        if (phi > 1.0): # Fuel rich but phi must be less than 1.44
          if (phi <= 1.44):
            n_co2=(11.8 - 8.2 * phi); n_co=(11.8 * phi - 11.8); n_h2o=4.6 * phi;
            hof_prod=n_co2 * hof_co2 + n_co * hof_co + n_h2o * hof_h2o + n_n2 * hof_n2
            h_prod=n_co2 * h_co2 + n_co * h_co + n_h2o * h_h2o + n_n2 * h_n2
            H_prod=hof_prod + h_prod
          else:
            print("Equivalence ratio must be less than 1.44")
            exit() # Stop execution

    # Value = H_prod - H_reac
    value = H_prod - H_reac
    return value

##***********************************************************************************##

## Main program

print(" ***  Adiabatic Flame Temperature computation for LPG mixture ***")
phi = float(input("Enter equivalence ratio :    "))
 
# phi = equivalence ratio

## Computation using Bisection method 

# Guesses for Bisection method
rlim=1000; # in Kelvins
llim=3000; # K

err=10.0 # Error counter

while(abs(err) > 10**-1):
    f_rlim=func_evaluate(rlim, phi)
    f_llim=func_evaluate(llim, phi)

    if(f_rlim * f_llim < 0):
        mid=(rlim + llim) * 0.5
    f_mid=func_evaluate(mid, phi)
    
    err = abs(rlim - mid)
    print("Error = "+ str(err) +"\n")
    if(f_llim * f_mid < 0):
        rlim=mid
    else:
        llim=mid

## Print result on screen
print("The value of Adiabatic flame temp is : " + str(mid))
