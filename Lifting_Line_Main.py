# -*- coding: utf-8 -*-
"""
Created on Thu May 16 20:20:15 2019

@author: ErikD
"""

#Still to do
#-Different CT and CP calculation for double rotor configuration

#IF IT ALL WORKS:
#-Add propellers

#TO DO: GET THE RIGHT RESULTS!
from scipy.interpolate import InterpolatedUnivariateSpline
import time
import numpy as np
from BEM_Main import BEM_Model as BEM
from GammaFunction import newGamma as GF
from Lifting_Line_Loads import loadBladeOverAllElements as LBE
from Lifting_Line_Rings import WakeGeometry as WG
from Vortex_Ring_System import induced_velocity_vortex_system as IV
import matplotlib.pyplot as plt

import pandas as pd

t_start = time.time()
#-------------------------------------------------------
#------(Start) Setup Block: Can be altered by user------
#-------------------------------------------------------

#Blade discretisation parameters
Ncp = 15
dis_mode = 'constant'
reload = False
order = 3

#Vortex wake parameters
n_r = 25
n_pr = 40
a_w_init = 0.1
case = 'turbine'           #Choice between 'turbine' and 'propeller'

#Iteration parameters
n_iterations = 20
Error_margin = 0.0001

BEM_compare = True         #Always false if Double_rotor = True
Plot = False
show_results_BEM = False

#Model setup for double rotor configuration
Double_rotor = True         #Set to True for double rotor configuration
L_sep = 1000000000          #Separation distance between two rotors expressed 
                            #in rotor diameter
phase_dif = 0               #Phase difference in degrees                            

#-----------------------------------------------
#-------------(End) Setup Block-----------------
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#--User should NOT alter any code hereonafter---
#-----------------------------------------------
#------(Start) Initialise problem setup---------
#-----------------------------------------------

#load in inflow conditions
if case == 'turbine':
    import Geometry_turbine as V
    Omega = V.U0*V.TSR/V.R
    airfoil = 'polar_DU95W180.xlsx'
    sheet   = "Blade1"
elif case == 'propeller':
    Omega = V.rpm/60*2*np.pi
    import Geometry_propeller as V
    airfoil = 'ARA_airfoil_data.xlsx'
    sheet   = "Sheet1"

#Pre-calculations in case of double rotor configuration
if Double_rotor:
    #BEM_compare = False
    L_sep = L_sep*(V.R*2)
    n_rotors = 2
else:
    n_rotors = 1

#Load in blade geometry
Geo_filename = 'GeoOptimal'+case+'.dat'
Geo = np.loadtxt(Geo_filename , dtype='float',skiprows = 1)
N_ann = len(Geo[:,0])

#Run BEM model program if number of controlpoints of lifting line program does
#not coincide with number of annuli in BEM program OR
#if we want to compare to BEM model calculations
BEM_t0 = time.time()
if N_ann != Ncp or reload:
    print('BEM model will be run first')
    BEM(Ncp, dis_mode, case, Optimise=True, show_results=False)
    print('')
    
    #Reload geo file with correct number of controlpoints
    Geo = np.loadtxt(Geo_filename , dtype='float',skiprows = 1)
   
if BEM_compare:
    results_BEM = BEM(Ncp, dis_mode, case, Optimise=False, show_results=show_results_BEM)  
BEM_tend = time.time()

print('BEM method took ', BEM_tend-BEM_t0,' seconds to complete')

#Continue with setting up lifting line model
controlpoints = Geo[:,0].reshape((len(Geo[:,0]),1))
Twist_cp = Geo[:,1].reshape(len(Geo[:,1]),1)
Chord_cp = Geo[:,2].reshape(len(Geo[:,2]),1)

#Determine chord, twist and radii at discretised locations
#The matrix Geo only yields the geometry at the controlpoints, which lie inside
#the discretised locations
R0 = V.rootradius_R
Rend = V.tipradius_R
if dis_mode == 'constant':
    R_disL = np.linspace(R0, Rend, Ncp+1)
elif dis_mode == 'cosine':
    R_mid = (Rend-R0)/2 
    R_disL = R0 + R_mid*(1-np.cos(np.pi/Ncp*np.arange(Ncp+1))) 

s_chord = InterpolatedUnivariateSpline(controlpoints, Chord_cp, k=order)
Chord_disL = s_chord(R_disL)
s_twist = InterpolatedUnivariateSpline(controlpoints, Twist_cp, k=order)
Twist_disL = s_twist(R_disL)


#Create twist, chord and controlpoint vectors that contain the controlpoints on
#all the blades
controlpoints_all = controlpoints
Twist_all_cp = Twist_cp
chord_all_cp = Chord_cp
for i in range(1,V.Nblades*n_rotors):
    controlpoints_all     = np.vstack((controlpoints_all, controlpoints))
    Twist_all_cp = np.vstack((Twist_all_cp,Twist_cp))
    chord_all_cp = np.vstack((chord_all_cp,Chord_cp))

#Load in polar data
data1=pd.read_excel(airfoil, sheet)
polar_alpha = data1['Alfa']
polar_cl = data1['Cl']
polar_cd = data1['Cd']

#Define total number of points in on lifting line wake
n_t = n_r*n_pr

print('Lifting line method for '+case+' has started')

#-----------------------------------------------
#-------(End) Initialise problem setup----------
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#-(Start) Main block 1: Initialise circulation--
#-----------------------------------------------

#Setup initial vortex wake structure
t_VW_0 = time.time()
if Double_rotor:
    Vortex_Wake = WG(V.U0, Omega, n_t, n_r, a_w_init, V.Nblades, R_disL*V.R, Chord_disL, Twist_disL, double=True, S_sep=L_sep, phase_dif=phase_dif, plot=True)
else:
    Vortex_Wake = WG(V.U0, Omega, n_t, n_r, a_w_init, V.Nblades, R_disL*V.R, Chord_disL, Twist_disL, plot=True)


t_VW_end = time.time()
print('Vortex wake geometry is determined in ', t_VW_end-t_VW_0,' seconds')

#Setup Biot-Savart induction matrix for Gamma = 1
if Double_rotor:
    [Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints_all*V.R, 1, double=True, phase_dif=phase_dif, d_sep=L_sep)
else:
    [Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints_all*V.R, 1)
t_ind_end = time.time()
print('Induced velocity matrices are calculated in ',t_ind_end-t_VW_end,' seconds')

#Initial estimate for circulation (with zero induced velocity)
t_Gamma_0 = time.time()
Uind = np.mat(np.zeros((Ncp*V.Nblades*n_rotors,1)))
Vind = np.mat(np.zeros((Ncp*V.Nblades*n_rotors,1)))
Wind = np.mat(np.zeros((Ncp*V.Nblades*n_rotors,1)))
if Double_rotor:
    Gamma = GF(V.rho, V.U0, Uind, Vind, Wind, Omega, controlpoints*V.R, Vortex_Wake, Twist_all_cp, polar_alpha, polar_cl, chord_all_cp, double=True, phase_dif=phase_dif)
else:
    Gamma = GF(V.rho, V.U0, Uind, Vind, Wind, Omega, controlpoints*V.R, Vortex_Wake, Twist_all_cp, polar_alpha, polar_cl, chord_all_cp)
t_Gamma_end = time.time()
print('Gamma is calculated in ',t_Gamma_end-t_Gamma_0,' seconds')

#-----------------------------------------------
#--(End) Main block 1: Initialise circulation---
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#---(Start) Main block 2: Iterate circulation---
#-----------------------------------------------

#Calculate induced velocity
i_iter = 1
Error_old = 1
while i_iter < n_iterations:
    Uind = Ind_Vel_Mat_u*Gamma
    Vind = Ind_Vel_Mat_v*Gamma
    Wind = Ind_Vel_Mat_w*Gamma
    
    #Update vortex rings with new a_w
    if Double_rotor:
        a_w = float(np.average(Uind)/V.U0)
        Vortex_Wake = WG(V.U0, Omega, n_t, n_r, a_w, V.Nblades, R_disL*V.R, Chord_disL, Twist_disL, double=True, S_sep=L_sep, phase_dif=phase_dif)
    else:
        a_w_weights = np.sin(np.pi/Ncp*np.arange(Ncp)).reshape((Ncp,1))
        a_w = float(np.average(Uind[:Ncp],axis=0,weights=a_w_weights)/V.U0)
        Vortex_Wake = WG(V.U0, Omega, n_t, n_r, a_w, V.Nblades, R_disL*V.R, Chord_disL, Twist_disL)
    
    if Double_rotor:
        [Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints_all*V.R, 1, double=True, phase_dif=phase_dif, d_sep=L_sep)
    else:
        [Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints_all*V.R, 1)
        
    #Determine new circulation
    if Double_rotor:
        GammaNew = GF(V.rho, V.U0, Uind, Vind, Wind, Omega, controlpoints*V.R, Vortex_Wake, Twist_all_cp, polar_alpha, polar_cl, chord_all_cp, double=True, phase_dif=phase_dif)
    else:
        GammaNew = GF(V.rho, V.U0, Uind, Vind, Wind, Omega, controlpoints*V.R, Vortex_Wake, Twist_all_cp, polar_alpha, polar_cl, chord_all_cp)

    #Calculate error
#    Error = np.sqrt(sum(np.multiply((GammaNew - Gamma),(GammaNew - Gamma))))/len(Gamma)

    #Different solution convergence check
    referror = max(abs(GammaNew))
    referror = max(referror,0.001)
    Error = max(abs(GammaNew-Gamma))
    Error = Error/referror   
    
    #Update the circulation
    Gamma = GammaNew
    
    if Error < Error_margin:
        print('Error is ', np.round(float(Error),5))
        break


    #Determine if gamma is converging
    if Error < Error_old:
        print('Error is ', np.round(float(Error),5), 'and Gamma is converging with: ',np.round(float(Error-Error_old),5))
    else:
        print('Error is ', np.round(float(Error),5), 'and Gamma is diverging with: +',np.round(float(Error-Error_old),5))
        
    Error_old = Error
    print(i_iter)
    i_iter += 1

#-----------------------------------------------
#----(End) Main block 2: Iterate circulation----
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#-----(Start) Main block 3: Calculate Loads-----
#-----------------------------------------------

#One simple function to obtain [fnorm, ftan, AoA, AngleInflow]
if Double_rotor:
    results = LBE(V.rho, V.U0, Uind, Vind, Wind, Omega, controlpoints*V.R, Twist_all_cp, polar_alpha, polar_cl, polar_cd, chord_all_cp, Vortex_Wake, double=False, phase_dif=phase_dif)
else:
    results = LBE(V.rho, V.U0, Uind, Vind, Wind, Omega, controlpoints*V.R, Twist_all_cp, polar_alpha, polar_cl, polar_cd, chord_all_cp, Vortex_Wake)

#-----------------------------------------------
#------(End) Main block 3: Calculate Loads------
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#---(Start) Post: Process results and display---
#-----------------------------------------------

#Average results if double rotor configuration is used
if Double_rotor:
    Fnorm = np.zeros((V.Nblades*n_rotors, Ncp))
    Ftan = np.zeros((V.Nblades*n_rotors, Ncp))
    alpha = np.zeros((V.Nblades*n_rotors, Ncp))
    inflow = np.zeros((V.Nblades*n_rotors, Ncp))
    Vax = np.zeros((V.Nblades*n_rotors, Ncp))
    Vtan = np.zeros((V.Nblades*n_rotors, Ncp))
    Gamma_blade = np.zeros((V.Nblades*n_rotors, Ncp))
    for i in range(V.Nblades*n_rotors):
        Fnorm[i,:] = results[0][i*Ncp:(i+1)*Ncp].reshape((1,Ncp))
        Ftan[i,:] = results[1][i*Ncp:(i+1)*Ncp].reshape((1,Ncp))
        alpha[i,:] = results[2][i*Ncp:(i+1)*Ncp].reshape((1,Ncp))
        inflow[i,:] = results[3][i*Ncp:(i+1)*Ncp].reshape((1,Ncp))
        Vax[i,:] = results[4][i*Ncp:(i+1)*Ncp].reshape((1,Ncp))
        Vtan[i,:] = results[5][i*Ncp:(i+1)*Ncp].reshape((1,Ncp))
        Gamma_blade[i,:] = Gamma[i*Ncp:(i+1)*Ncp].reshape((1,Ncp))
    Fnorm = np.average(Fnorm, axis=0)
    Ftan = np.average(Ftan, axis=0)
    alpha = np.average(alpha, axis=0)
    inflow = np.average(inflow, axis=0)
    Vax = np.average(Vax, axis=0)
    Vtan = np.average(Vtan, axis=0)
    Gamma_blade = np.average(Gamma_blade, axis=0)
        
else:
    Fnorm = results[0][:Ncp]
    Ftan = results[1][:Ncp]
    alpha = results[2][:Ncp] 
    inflow = results[3][:Ncp]
    Vax = results[4][:Ncp]
    Vtan = results[5][:Ncp]   
    Gamma_blade = Gamma[:Ncp]

#Calculate thrust and power coefficient
dr = (R_disL[1:]-R_disL[:-1])*V.R
CT = np.sum(dr*np.array(Fnorm).reshape((Ncp,))*V.Nblades/(0.5*V.U0**2*V.rho*np.pi*V.R**2))
CP = np.sum(dr*np.array(Ftan).reshape((Ncp,))*np.array(controlpoints).reshape((Ncp,))*V.R*V.Nblades*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))
print('LL CT = ', np.round(CT,4))
print('LL CP = ', np.round(CP,4))
if Plot:
    #Axial induction factor versus radius  
    a_i = 1 - Vax/V.U0
    aline = Vtan.reshape((Ncp,1))/(Omega*controlpoints*V.R) - 1
    
    fig_a = plt.figure(figsize=(12,6))
    ax_a = plt.gca()
    plt.title('Axial induction factor')
    plt.plot(controlpoints, a_i, 'r-x', label=r'$a$ - LL')
    plt.plot(controlpoints, aline, 'g--x', label=r'$a^,$')
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()
    
    #Normal and tangential forces
    fig_force = plt.figure(figsize=(12,6))
    ax_force = plt.gca()
    plt.title(r'Normal and tangential force, non-dimensionalised by $\frac{1}{2} \rho U_\infty^2 R$')
    plt.plot(controlpoints, Fnorm/(0.5*V.U0**2*V.rho*V.R), 'r-x', label=r'Fnorm - LL')
    plt.plot(controlpoints, Ftan/(0.5*V.U0**2*V.rho*V.R), 'g--x', label=r'Ftan - LL')
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()
        
    #Circulation distribution
    fig_circ = plt.figure(figsize=(12,6))
    ax_circ = plt.gca()
    plt.title(r'Circulation distribution, non-dimensionalised by $\frac{\pi U_\infty^2}{\Omega * NBlades}$')
    plt.plot(controlpoints, Gamma_blade/(np.pi*V.U0**2/(V.Nblades*Omega)), 'r-x', label=r'$\Gamma$ - LL')
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.show()
    
    #Inflow distribution
    fig_inflow = plt.figure(figsize=(12,6))
    ax_inflow = plt.gca()
    plt.title('Angle distribution')
    plt.plot(controlpoints, inflow*180/np.pi, 'r-x', label='Inflowangle - LL')
    plt.plot(controlpoints, Twist_cp, 'g-s', label='Twist - LL')
    plt.plot(controlpoints, alpha, 'k-x', label=r'$\alpha$ - LL')
    plt.grid()
    plt.ylabel('(deg)')
    plt.xlabel('r/R')
    plt.legend()
    plt.show()

if BEM_compare:
    CT_BEM = np.sum(dr*results_BEM[:,3]*V.Nblades/(0.5*V.U0**2*V.rho*np.pi*V.R**2))
    CP_BEM = np.sum(dr*results_BEM[:,4]*np.array(controlpoints).reshape((Ncp,))*V.R*V.Nblades*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))

    print('BEM CT = ', np.round(CT_BEM,4))
    print('BEM CP = ', np.round(CP_BEM,4))
    if Plot:
        ax_a.plot(controlpoints, results_BEM[:,0], 'k-o', label=r'$a$ - BEM')
        ax_a.plot(controlpoints, results_BEM[:,1], 'k--o', label=r'$a^,$ - BEM')
        ax_a.legend()
    
        ax_force.plot(controlpoints, results_BEM[:,3]/(0.5*V.U0**2*V.rho*V.R), 'k-o', label=r'Fnorm - BEM')
        ax_force.plot(controlpoints, results_BEM[:,4]/(0.5*V.U0**2*V.rho*V.R), 'k--o', label=r'Ftan - BEM')
        ax_force.legend()    
        
        ax_circ.plot(controlpoints, results_BEM[:,5]/(np.pi*V.U0**2/(V.Nblades*Omega)), 'k--o', label=r'$\Gamma$ - BEM')
        ax_circ.legend()
        
        ax_inflow.plot(controlpoints, results_BEM[:,7]*180/np.pi, 'r--o', label='Inflowangle - BEM')
        ax_inflow.plot(controlpoints, results_BEM[:,6], 'k--o', label=r'$\alpha$ - BEM')
        ax_inflow.legend()
#-----------------------------------------------
#----(End) Post: Process results and display----
#-----------------------------------------------
t_end = time.time()
print('Total elapsed time is ',t_end-t_start,' seconds')