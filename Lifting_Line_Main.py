# -*- coding: utf-8 -*-
"""
Created on Thu May 16 20:20:15 2019

@author: ErikD
"""

#Still to do
#-Obtain controlpoints and points at discretised locations DONE
#-Block calculate results                                  DONE
#-Block display results                                    DONE
#-Calculation of induction factor to implement in the loop DONE

#TO DO: DEBUG AND GET IT WORKING!
import numpy as np
from GammaFunction import newGamma as GF
from Lifting_Line_Loads import loadBladeOverAllElements as LBE
from Lifting_Line_Rings import WakeGeometry as WG
from Vortex_Ring_System import induced_velocity_vortex_system as IV
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

#-------------------------------------------------------
#------(Start) Setup Block: Can be altered by user------
#-------------------------------------------------------

#Vortex wake parameters
n_r = 1
n_t = 50
a_w_init = 0.1
case = 'Turbine'

#Iteration parameters
n_iterations = 1000
Error_margin = 0.0001

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
if case == 'Turbine':
    import Geometry_turbine as V
    Omega = V.U0*V.TSR/V.R
    airfoil = 'polar_DU95W180.xlsx'
    sheet   = "Blade1"
elif case == 'Propeller':
    Omega = V.rpm/60*2*np.pi
    import Geometry_propeller as V
    airfoil = 'ARA_airfoil_data.xlsx'
    sheet   = "Sheet1"

#Load in blade geometry
Geo_filename = 'GeoOptimal'+case+'.dat'
Geo = np.loadtxt(Geo_filename , dtype='float',skiprows = 1)
controlpoints = Geo[:,0]
Twist_cp = Geo[:,1]
Chord_cp = Geo[:,2]
Ncp = len(controlpoints)

#Determine chord, twist and radii at discretised locations
#The matrix Geo only yields the geometry at the controlpoints, which lie inside
#the discretised locations
R0 = controlpoints[0]-(controlpoints[1]-controlpoints[0])/2
Rend = controlpoints[-1]+(controlpoints[-1]-controlpoints[-2])/2
R_disL = np.linspace(R0, Rend, Ncp+1)     #This is possible since we will only
                                           #treat linearly discretised blades

Twist0 = Geo[0,1]-(Geo[1,1]-Geo[0,1])/2
Twistend = Geo[-1,1]+(Geo[-1,1]-Geo[-2,1])/2
Twist_disL = np.zeros((Ncp+1,1))

Chord0 = Geo[0,2]-(Geo[1,2]-Geo[0,2])/2
Chordend = Geo[-1,2]+(Geo[-1,2]-Geo[-2,2])/2
Chord_disL = np.zeros((Ncp+1,1))
for i in range(1,Ncp):
    Twist_disL[i] = (Geo[i-1,1]+Geo[i,1])/2
    Chord_disL[i] = (Geo[i-1,2]+Geo[i,1])/2
Twist_disL[-1] = Twistend
Chord_disL[-1] = Chordend

#Create twist, chord and controlpoint vectors that contain the controlpoints on
#all the blades
controlpoints_all = controlpoints
Twist_all_cp = Twist_cp
chord_all_cp = Chord_cp
for i in range(1,V.Nblades):
    controlpoints_all     = np.vstack((controlpoints_all, controlpoints))
    Twist_all_cp = np.vstack((Twist_all_cp,Twist_cp))
    chord_all_cp = np.vstack((chord_all_cp,Chord_cp))

#Load in polar data
data1=pd.read_excel(airfoil, sheet)
polar_alpha = data1['Alfa']
polar_cl = data1['Cl']
polar_cd = data1['Cd']

print('Lifting line method for '+case+' has started')

#-----------------------------------------------
#-------(End) Initialise problem setup----------
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#-(Start) Main block 1: Initialise circulation--
#-----------------------------------------------

#Setup initial vortex wake structure
Vortex_Wake = WG(V.U0, Omega, n_t, n_r, a_w_init, V.Nblades, R_disL*V.R, Chord_disL, Twist_disL)

# =============================================================================
# ###########################################STUCK HERE
# =============================================================================

#Setup Biot-Savart induction matrix for Gamma = 1
[Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints*V.R, 1)

#Initial estimate for circulation (with zero induced velocity)
uind = np.mat(np.zeros((Ncp*V.Nblades,1)))
vind = np.mat(np.zeros((Ncp*V.Nblades,1)))
wind = np.mat(np.zeros((Ncp*V.Nblades,1)))
Gamma = GF(V.rho, V.U0, uind, vind, wind, Omega, controlpoints_all*V.R, Vortex_Wake, Twist_all_cp, polar_alpha, polar_cl, chord_all_cp)#Create a function for determining gamma (Output should be a numpy matrix!)

#-----------------------------------------------
#--(End) Main block 1: Initialise circulation---
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#---(Start) Main block 2: Iterate circulation---
#-----------------------------------------------

#Calculate induced velocity
i_iter = 1
while i_iter < n_iterations:
    Uind = Ind_Vel_Mat_u*Gamma
    Vind = Ind_Vel_Mat_v*Gamma
    Wind = Ind_Vel_Mat_w*Gamma
    
    #Update vortex rings with new a_w
    a_w = 1 - (np.mean(Uind)+V.U0)/V.U0
    Vortex_Wake = WG(V.U0, Omega, n_t, n_r, a_w, V.Nblades, R_disL*V.R, Chord_disL, Twist_disL)
    [Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints*V.R, 1)
    
    #Determine new circulation
    GammaNew = GF(V.rho, V.U0, Uind, Vind, Wind, Omega, controlpoints_all*V.R, Vortex_Wake, Twist_all_cp, polar_alpha, polar_cl, chord_all_cp)
    
    #Calculate error
    Error = np.sqrt(sum(np.multiply((GammaNew - Gamma),(GammaNew - Gamma)))/len(Gamma))

    Gamma = GammaNew
    if Error < Error_margin:
        break
    
    print(i)
    i += 1


#-----------------------------------------------
#----(End) Main block 2: Iterate circulation----
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#-----(Start) Main block 3: Calculate Loads-----
#-----------------------------------------------

#One simple function to obtain [fnorm, ftan, AoA, AngleInflow]
results = LBE(V.rho, V.U0, Uind[:Ncp], Vind[:Ncp], Omega, controlpoints*V.R, Twist_cp, polar_alpha, polar_cl, polar_cd, Chord_cp)
    
#-----------------------------------------------
#------(End) Main block 3: Calculate Loads------
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#---(Start) Post: Process results and display---
#-----------------------------------------------

#Calculate thrust and power coefficient
dr = (R_disL[1:]-R_disL[:-1])*V.R
CT = np.sum(dr*results[:,0]*V.Nblades/(0.5*V.U0**2*V.rho*np.pi*V.R**2))
CP = np.sum(dr*results[:,1]*controlpoints*V.Nblades*V.R*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))

#Axial induction factor over radius
a_i = 1 - (Uind+V.U0)/V.U0
fig_a = plt.figure(figsize=(12,6))
plt.title('Axial induction factor')
plt.plot(results[:,2], a_i, 'r-', label=r'$a$')
#plt.plot(results[:,2], results[:,1], 'g--', label=r'$a^,$')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Normal and tangential forces
fig_force = plt.figure(figsize=(12,6))
plt.title(r'Normal and tangential force, non-dimensionalised by $\frac{1}{2} \rho U_\infty^2 R$')
plt.plot(controlpoints, results[:,0]/(0.5*V.U0**2*V.rho*V.R), 'r-', label=r'Fnorm')
plt.plot(controlpoints, results[:,1]/(0.5*V.U0**2*V.rho*V.R), 'g--', label=r'Ftan')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()
    
#Circulation distribution
fig_circ = plt.figure(figsize=(12,6))
plt.title(r'Circulation distribution, non-dimensionalised by $\frac{\pi U_\infty^2}{\Omega * NBlades}$')
plt.plot(controlpoints, Gamma[:Ncp]/(np.pi*V.U0**2/(V.Nblades*Omega)), 'r-', label=r'$\Gamma$')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Inflow distribution
fig_inflow = plt.figure(figsize=(12,6))
plt.title('Angle distribution')
plt.plot(controlpoints, results[:,3]*180/np.pi, 'r-', label='Inflowangle')
plt.plot(controlpoints, Twist_cp, 'g--', label='Twist')
plt.plot(controlpoints, results[:,2], 'k-.', label=r'$\alpha$')
plt.grid()
plt.ylabel('(deg)')
plt.xlabel('r/R')
plt.legend()
plt.show()
    
#-----------------------------------------------
#----(End) Post: Process results and display----
#-----------------------------------------------