# -*- coding: utf-8 -*-
#----------------------------------------
#Discretisation
#----------------------------------------


import matplotlib.pyplot as plt
import numpy as np
from BEM_Code import solveStreamtube as SS
import pandas as pd

#Setup
N_ann_reg = 900              #Number of annuli
N_ann_list = [1, 5, 10, 50, 250, 500, 750, 800, 950]
Dis_method = 'constant'   #Discretisation method, ('cosine' or 'constant')
case = 'turbine'        #'turbine' or 'propeller'
convergence_test = True


#Import correct geometry file
if case == 'turbine':
    import BEM_Geometry_turbine as V
    print('Turbine')
elif case == 'propeller':
    import BEM_Geometry_propeller as V
    print('Propellor')
else:
    print('You are DOOMED')


        

if Dis_method == 'constant':
    r_R = np.linspace(0.2,1,N_ann_reg+1)


#Calculation
if case == 'turbine':
    Omega =V.U0*V.TSR/V.R
    airfoil = 'polar_DU95W180.xlsx' # Where are we using this ?
    sheet   = "Blade1"
    #chord_distribution = 3*(1-r_R)+1 # !!need optimised geometry
    #twist_distribution = -14*(1-r_R)+V.Pblade # need optimised geometry
    
elif case == 'propeller':
    Omega = V.rpm/60*2*np.pi
    airfoil = 'ARA_airfoil_data.xlsx'
    sheet   = "Sheet1"
    #chord_distribution = 0.18-0.06*r_R # need optimised geometry
    #twist_distribution = -50*r_R+V.Pblade # need optimised geometry


data1=pd.read_excel( airfoil, sheet)#airfoil, header=0,
                    #names = ["Alfa", "Cl", "Cd", "Cm"], sep='\s+')
polar_alpha = data1['Alfa']#[:]
polar_cl = data1['Cl']#[:]
polar_cd = data1['Cd']#[:]

results = np.zeros((N_ann_reg,13))                                     #***CHANGE***
r_R_centroid = np.zeros((N_ann_reg,6))                                #***CHANGE***
for annulus in range(0,N_ann_reg):                                    #***CHANGE***
    r1_R = r_R[annulus]                                           #***CHANGE***
    r2_R = r_R[annulus+1]                                         #***CHANGE***
    r_R_centroid[annulus] = (r1_R + r2_R)/2                       #***CHANGE***              
    #chord = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R.transpose(), chord_distribution.transpose())#***CHANGE***
    #twist = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R, twist_distribution)                        #***CHANGE***
    results[annulus,:] = SS(case, V.rho, V.U0, r1_R, r2_R, V.rootradius_R, V.tipradius_R, Omega, V.R, V.Nblades, chord, twist, polar_alpha, polar_cl, polar_cd)

