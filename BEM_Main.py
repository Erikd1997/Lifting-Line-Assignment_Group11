# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 19:45:23 2019

@author: ErikD
"""

import matplotlib.pyplot as plt
import numpy as np
from BEM_Code import solveStreamtube as SS
from BEM_Code import Optimise as Opt
import pandas as pd

#--------------------------------------------------------------
#----------Start: user should only change this block-----------
#--------------------------------------------------------------
#Setup
N_ann_reg = 50         #Number of annuli
N_ann_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 12]
Dis_method = 'constant'   #Discretisation method, ('cosine' or 'constant')
case = 'turbine'        #'turbine' or 'propeller'
convergence_test = False
Optimise = True
#--------------------------------------------------------------
#----------End: user should only change this block-------------
#--------------------------------------------------------------

#--------------------------------------------------------------
#----------User should not ALTER code hereafter----------------
#--------------------------------------------------------------
 
#Import correct geometry file
if case == 'turbine':
    import BEM_Geometry_turbine as V
    print('Turbine')
elif case == 'propeller':
    import BEM_Geometry_propeller as V
    print('Propeller')

#----------------------------------------
#Convergence Block
#----------------------------------------
if convergence_test:
    CP_conv = np.zeros(np.shape(N_ann_list))
    CT_conv = np.zeros(np.shape(N_ann_list))
    #Discretisation
    
    index = 0
    for N_ann in N_ann_list:
        N_ann = int(N_ann)
        print(N_ann)
        if Dis_method == 'cosine':
            r_R = np.linspace(0,1,N_ann+1)
            for n in range(1,N_ann+2):
                theta = np.pi/N_ann * (n-1)
                r_R[n-1] = V.rootradius_R+((0.5-V.rootradius_R/2)-(0.5-V.rootradius_R/2)*np.cos(theta))
        elif Dis_method == 'constant':
            r_R = np.linspace(V.rootradius_R,1,N_ann+1)
        
        
        #Calculation
        if case == 'turbine':
            Omega =V.U0*V.TSR/V.R
            airfoil = 'polar_DU95W180.xlsx' # Where are we using this ?
            sheet   = "Blade1"
            chord_distribution = 3*(1-r_R)+1 # meters
            twist_distribution = -14*(1-r_R)+V.Pblade # degrees !!SEND QUESTION ABOUT -14 or +14
            
        elif case == 'propeller':
            Omega = V.rpm/60*2*np.pi
            airfoil = 'ARA_airfoil_data.xlsx'
            sheet   = "Sheet1"
            chord_distribution = 0.18-0.06*r_R # meters
            twist_distribution = -50*(r_R)+V.Pblade # degrees
        
        
        data1=pd.read_excel( airfoil, sheet)#airfoil, header=0,
                            #names = ["Alfa", "Cl", "Cd", "Cm"], sep='\s+')
        polar_alpha = data1['Alfa']#[:]
        polar_cl = data1['Cl']#[:]
        polar_cd = data1['Cd']#[:]
        
        results = np.zeros((N_ann,14))                                    #***CHANGE***
        r_R_centroid = np.zeros((N_ann,6))                                #***CHANGE***
        for annulus in range(0,N_ann):                                    #***CHANGE***
            r1_R = r_R[annulus]                                           #***CHANGE***
            r2_R = r_R[annulus+1]                                         #***CHANGE***
            r_R_centroid[annulus] = (r1_R + r2_R)/2                       #***CHANGE***              
            chord = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R.transpose(), chord_distribution.transpose())#***CHANGE***
            twist = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R.transpose(), twist_distribution.transpose())                        #***CHANGE***
            results[annulus,:] = SS(case, V.rho, V.U0, r1_R, r2_R, V.rootradius_R, V.tipradius_R, Omega, V.R, V.Nblades, chord, twist, polar_alpha, polar_cl, polar_cd, V.TSR)
        
        areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*V.R**2
        dr = (r_R[1:]-r_R[:-1])*V.R
        CT_conv[index] = np.sum(dr*results[:,3]*V.Nblades/(0.5*V.U0**2*V.rho**np.pi*V.R**2))
        CP_conv[index] = np.sum(dr*results[:,4]*results[:,2]*V.Nblades*V.R*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))
        index = index + 1   
        #***CHANGE***
        
#----------------------------------------
#Main Block
#----------------------------------------
if Dis_method == 'cosine':
    r_R = np.linspace(0,1,N_ann_reg+1)
    for n in range(1,N_ann_reg+2):
        theta = np.pi/N_ann_reg * (n-1)
        r_R[n-1] = V.rootradius_R+((0.5-V.rootradius_R/2)-(0.5-V.rootradius_R/2)*np.cos(theta))
elif Dis_method == 'constant':
    r_R = np.linspace(V.rootradius_R,1,N_ann_reg+1)


#Calculation
if case == 'turbine':
    Omega =V.U0*V.TSR/V.R
    airfoil = 'polar_DU95W180.xlsx' # Where are we using this ?
    sheet   = "Blade1"
    chord_distribution = 3*(1-r_R)+1 # meters
    twist_distribution = -14*(1-r_R)+V.Pblade # degrees !!SEND QUESTION ABOUT -14 or +14
    
elif case == 'propeller':
    Omega = V.rpm/60*2*np.pi
    airfoil = 'ARA_airfoil_data.xlsx'
    sheet   = "Sheet1"
    chord_distribution = 0.18-0.06*r_R # meters
    twist_distribution = 50*(r_R)+V.Pblade # degrees


data1=pd.read_excel( airfoil, sheet)#airfoil, header=0,
                    #names = ["Alfa", "Cl", "Cd", "Cm"], sep='\s+')
polar_alpha = data1['Alfa']#[:]
polar_cl = data1['Cl']#[:]
polar_cd = data1['Cd']#[:]

results = np.zeros((N_ann_reg,16))  
if Optimise:
    results_opt = np.zeros((N_ann_reg,9))
r_R_centroid = np.zeros((N_ann_reg,6))                                #***CHANGE***
for annulus in range(0,N_ann_reg):                                    #***CHANGE***
    r1_R = r_R[annulus]                                           #***CHANGE***
    r2_R = r_R[annulus+1]                                         #***CHANGE***
    r_R_centroid[annulus] = (r1_R + r2_R)/2                       #***CHANGE***    
    chord = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R.transpose(), chord_distribution.transpose())#***CHANGE***
    twist = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R, twist_distribution)                        #***CHANGE***
    results[annulus,:] = SS(case, V.rho, V.U0, r1_R, r2_R, V.rootradius_R, V.tipradius_R, Omega, V.R, V.Nblades, chord, twist, polar_alpha, polar_cl, polar_cd, V.TSR)
    if Optimise:
        Area = np.pi*((r2_R*V.R)**2-(r1_R*V.R)**2) #  area streamtube
        fnorm_opti = V.CT_opt*0.5*V.rho*V.U0**2*Area/V.Nblades/V.R/(r2_R-r1_R)
        results_opt[annulus,:] = Opt(case, fnorm_opti, V.U0, r1_R, r2_R, V.R, V.rootradius_R, V.tipradius_R, Omega, V.Nblades, V.rho, polar_alpha, polar_cl, polar_cd, V.TSR)

#Location Upwind
P_inf = 101325 #Pascal
P0 = P_inf + (0.5*V.rho*V.U0)
P_Total_Uinf = P0*np.ones(N_ann_reg) 

#Location Upwind Rotor
One = np.ones(N_ann_reg)
P_RotorU = P_Total_Uinf - (0.5*V.rho*(results[:,14]*2+results[:,15]*2))
P_Total_RU = P_Total_Uinf

#Location Downwind Rotor
if case == 'turbine':
    P_RotorD = P_RotorU - (results[:,3]/areas)
elif case == 'propeller':
    P_RotorD = P_RotorU + (results[:,3]/areas)
        
P_Total_RD = P_RotorD + (0.5*V.rho*(results[:,14]*2+results[:,15]*2))
        
#Location DownWind
P_Total_Dinf = P_Total_RD

#Plot results
areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*V.R**2
dr = (r_R[1:]-r_R[:-1])*V.R
CT = np.sum(dr*results[:,3]*V.Nblades/(0.5*V.U0**2*V.rho*np.pi*V.R**2))
CP = np.sum(dr*results[:,4]*results[:,2]*V.Nblades*V.R*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))
if Optimise:
    CT_opt = np.sum(dr*results_opt[:,3]*V.Nblades/(0.5*V.U0**2*V.rho*np.pi*V.R**2))
    CP_opt = np.sum(dr*results_opt[:,4]*results_opt[:,2]*V.Nblades*V.R*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))

print('CP is ', CP)
print('CT is ', CT)
if Optimise:
    print('Optimal CT is ', CT_opt)
    print('Optimal CP is ', CP_opt)

#Axial and tangential induction factor over radius
fig_a = plt.figure(figsize=(12,6))
plt.title('Axial and tangential induction factor')
plt.plot(results[:,2], results[:,0], 'r-', label=r'$a$')          #***CHANGE***
plt.plot(results[:,2], results[:,1], 'g--', label=r'$a^,$')       #***CHANGE***
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Normal and tangential forces
fig_force = plt.figure(figsize=(12,6))
plt.title(r'Normal and tangential force, non-dimensionalised by $\frac{1}{2} \rho U_\infty^2 R$')
plt.plot(results[:,2], results[:,3]/(0.5*V.U0**2*V.rho*V.R), 'r-', label=r'Fnorm')#***CHANGE***
plt.plot(results[:,2], results[:,4]/(0.5*V.U0**2*V.rho*V.R), 'g--', label=r'Ftan')#***CHANGE***
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()
    
#Circulation distribution
fig_circ = plt.figure(figsize=(12,6))
plt.title(r'Circulation distribution, non-dimensionalised by $\frac{\pi U_\infty^2}{\Omega * NBlades}$')
plt.plot(results[:,2], results[:,5]/(np.pi*V.U0**2/(V.Nblades*Omega)), 'r-', label=r'$\Gamma$')#***CHANGE***
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#alpha distribution
fig_alpha = plt.figure(figsize=(12,6))
plt.title('Alpha distribution')
plt.plot(results[:,2], results[:,6], 'r-', label=r'$\alpha$')#***CHANGE***
plt.grid()
plt.ylabel(r'$\alpha$ (deg)')
plt.xlabel('r/R')
plt.legend()
plt.show()

#Inflow distribution
fig_inflow = plt.figure(figsize=(12,6))
plt.title('Angle distribution')
plt.plot(results[:,2], results[:,7]*180/np.pi, 'r-', label='Inflowangle')#***CHANGE***
plt.plot(r_R, twist_distribution, 'g--', label='Twist')
plt.plot(results[:,2], results[:,6], 'k-.', label=r'$\alpha$')
plt.grid()
plt.ylabel('(deg)')
plt.xlabel('r/R')
plt.legend()
plt.show()

#Ct and Cn distribution
fig_ct_cn = plt.figure(figsize=(12,6))
plt.title(r'Distribution of $C_{n}$ and $C_{t}$')
plt.plot(results[:,2], results[:,9], 'r-' , label=r'$C_{n}$')#***CHANGE***
plt.plot(results[:,2], results[:,8], 'g--', label=r'$C_{t}$')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Cq distribution
fig_cq = plt.figure(figsize=(12,6))
plt.title(r'Distribution of $C_{Q}$')
plt.plot(results[:,2], results[:,13], 'r-' , label=r'$C_{n}$')#***CHANGE***
#plt.plot(results[:,2], results[:,8], 'g--', label=r'$C_{t}$')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Prandtl correction
fig_prandtl = plt.figure(figsize=(12,6))
plt.title(r'Prandtl correction distribution')
plt.grid()
plt.plot(results[:,2], results[:,10], 'r-', label='Prandtl')
plt.plot(results[:,2], results[:,11], 'g.', label='Prandtl tip')
plt.plot(results[:,2], results[:,12], 'b.', label='Prandtl root')
plt.xlabel('r/R')
plt.legend()
plt.show()

#Convergence plot
if convergence_test:
    fig_prandtl = plt.figure(figsize=(12,6))
    plt.title(r'Convergence error of $C_{t}$ and $C_{p}$')
    plt.grid()
    plt.yscale('log')
    plt.plot(N_ann_list[1:-1], abs(CT_conv[-1] - CT_conv[1:-1]), 'r-', label=r'$C_{t}$')
    plt.plot(N_ann_list[1:-1], abs(CP_conv[-1] - CP_conv[1:-1]), 'g--', label=r'$C_{p}$')
    plt.xlabel('# of annuli')
    plt.legend()
    plt.show()
    
    
    
#Axial and tangential induction factor over radius (Optimal case)
if Optimise:
    fig_a_opt = plt.figure(figsize=(12,6))
    plt.title('Axial and tangential induction factor for optimised turbine')
    plt.plot(results_opt[:,2], results_opt[:,0], 'r-', label=r'$a$')          #***CHANGE***
    plt.plot(results_opt[:,2], results_opt[:,1], 'g--', label=r'$a^,$')       #***CHANGE***
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()

#Twist distribution (Optimal case)
    fig_a_opt = plt.figure(figsize=(12,6))
    plt.title('Twist distribution for optimised turbine')
    plt.plot(results_opt[:,2], results_opt[:,6]*180/np.pi, 'r-.', label='inflowangle')          #***CHANGE***
    plt.plot(results_opt[:,2], results_opt[:,5], 'k-', label=r'$\alpha_{opt}$')
    plt.plot(results_opt[:,2], results_opt[:,7], 'g--', label='twist')       #***CHANGE***
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()

    fig_force_opt = plt.figure(figsize=(12,6))
    plt.title(r'Normal and tangential force for optimised turbine, non-dimensionalised by $\frac{1}{2} \rho U_\infty^2 R$')
    plt.plot(results_opt[:,2], results_opt[:,3]/(0.5*V.U0**2*V.rho*V.R), 'k-', label=r'$\alpha_{opt}$')
    plt.plot(results_opt[:,2], results_opt[:,4]/(0.5*V.U0**2*V.rho*V.R), 'g--', label='optimal')       #***CHANGE***
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()

    fig_chord_opt = plt.figure(figsize=(12,6))
    plt.title(r'Chord distribution for optimised turbine')
    plt.plot(results_opt[:,2], results_opt[:,8], 'r-', label='chord')       #***CHANGE***
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()   
    
    fig_circ = plt.figure(figsize=(12,6))
    plt.title(r'Circulation distribution for optimised turbine, non-dimensionalised by $\frac{\pi U_\infty^2}{\Omega * NBlades}$')
    plt.plot(results_opt[:,2], results_opt[:,5]/(np.pi*V.U0**2/(V.Nblades*Omega)), 'r-', label=r'$\Gamma$')#***CHANGE***
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()
    
#-----------------------------------------------------------------------------
#**PLOT**#

if convergence_test==False:
    # Stagnation Enthalphy
    fig_a = plt.figure(figsize=(12,6))
    plt.title('Stagnation Enthalphy at Four Locations')
    plt.plot(results[:,2], P_Total_Uinf , 'ro-', label=r'Infinitly Upwind')          #**CHANGE**
    plt.plot(results[:,2], P_Total_RU , 'b.-', label=r'Rotor Upwind') 
    plt.plot(results[:,2], P_Total_RD , 'g-', label=r'Rotor Downwind') 
    plt.plot(results[:,2], P_Total_Dinf , 'y*', label=r'Infinitly Downwind') 
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()