# -*- coding: utf-8 -*-
"""
Created on Thu May 16 20:20:15 2019

@author: ErikD
"""
import numpy as np
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
n_t = 400
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
    Omega = V.rpm/60*2*np.pi
    airfoil = 'polar_DU95W180.xlsx'
    sheet   = "Blade1"
elif case == 'Propeller':
    import Geometry_propeller as V
    airfoil = 'ARA_airfoil_data.xlsx'
    sheet   = "Sheet1"

#Load in blade geometry
Geo_filename = 'GeoOptimal'+case+'.dat'
Geo = np.loadtxt(Geo_filename , dtype='float',skiprows = 1)
print('Lifting line method for '+case+' has started')

#Load in polar data
data1=pd.read_excel(airfoil, sheet)
polar_alpha = data1['Alfa']
polar_cl = data1['Cl']
polar_cd = data1['Cd']

#-----------------------------------------------
#-------(End) Initialise problem setup----------
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#-(Start) Main block 1: Initialise circulation--
#-----------------------------------------------

#Setup initial vortex wake structure
Vortex_Wake = WG(V.U_infty, Omega, n_t, n_r, a_w, V.Nblades, V.R, chord, Twist)

#Setup Biot-Savart induction matrix for Gamma = 1
[Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints, 1)

#Initial estimate for circulation (Only contributions from)
Vaxial = V.U_infty
Vtan = Omega*controlpoints
Vp = np.sqrt(Vaxial^2 + Vtan^2)
Gamma = np.matrix(0.5*Chord.reshape((len(controlpoints*V.Nblades),1))*Vp*Cl) #Create a function for determining gamma (Output should be a numpy matrix!)

#-----------------------------------------------
#--(End) Main block 1: Initialise circulation---
#-----------------------------------------------
###############################################################################
#-----------------------------------------------
#---(Start) Main block 2: Iterate circulation---
#-----------------------------------------------

#Calculate induced velocity   (Will possibly become a loop)
i_iter = 1
while i_iter < n_iterations:
    Uind = Ind_Vel_Mat_u*Gamma
    Vind = Ind_Vel_Mat_v*Gamma
    Wind = Ind_Vel_Mat_w*Gamma
    
    #Update vortex rings with new a_w
    Vortex_Wake = WG(V.U_infty, Omega, n_t, n_r, a_w, V.Nblades, V.R, chord, Twist)
    [Ind_Vel_Mat_u, Ind_Vel_Mat_v, Ind_Vel_Mat_w] = IV(Vortex_Wake, controlpoints, 1)
    
    GammaNew = 1#Create new gamma (function)
    Error = GammaNew - Gamma

    Gamma = Gammanew
    if Error < Error_margin:
        break



#-----------------------------------------------
#----(End) Main block 2: Iterate circulation----
#-----------------------------------------------
###############################################################################