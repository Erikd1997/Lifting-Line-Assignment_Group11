# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 19:45:23 2019

@author: ErikD
"""

from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import numpy as np
from BEM_Code import solveStreamtube as SS
from BEM_Code import Optimise as Opt
import pandas as pd

def BEM_Model(N_ann_reg, Dis_method, case, Optimise=False, show_results = False):
    #--------------------------------------------------------------
    #----------Start: user should only change this block-----------
    #--------------------------------------------------------------
    #Setup
    order = 3
    plot = show_results
    SaveGeo = True
    #--------------------------------------------------------------
    #----------End: user should only change this block-------------
    #--------------------------------------------------------------
    
    #--------------------------------------------------------------
    #----------User should not ALTER code hereafter----------------
    #--------------------------------------------------------------
     
    #Import correct geometry file
    if case == 'turbine':
        import Geometry_turbine as V
        print('BEM method for the turbine has commenced with N =', N_ann_reg)
    elif case == 'propeller':
        import Geometry_propeller as V
        print('BEM method for the propeller has commenced with N =', N_ann_reg)
        
    #----------------------------------------
    #Main Block
    #----------------------------------------
    if Dis_method == 'cosine':
        r_R = np.linspace(0,1,N_ann_reg+1)
        for n in range(1,N_ann_reg+2):
            theta = np.pi/N_ann_reg * (n-1)
            r_R[n-1] = V.rootradius_R+((0.5-V.rootradius_R/2)-(0.5-V.rootradius_R/2)*np.cos(theta))
    elif Dis_method == 'constant':
        r_R = np.linspace(V.rootradius_R,V.tipradius_R,N_ann_reg+1)
        
    #Calculation
    if case == 'turbine':
        Omega =V.U0*V.TSR/V.R
        airfoil = 'polar_DU95W180.xlsx' # Where are we using this ?
        sheet   = "Blade1"

#        chord_distribution = 3*(1-r_R)+1 # meters
#        twist_distribution = -14*(1-r_R)+V.Pblade # degrees !!SEND QUESTION ABOUT -14 or +14
#        
    elif case == 'propeller':
        Omega = V.rpm/60*2*np.pi
        airfoil = 'ARA_airfoil_data.xlsx'
        sheet   = "Sheet1"
#        chord_distribution = 0.18-0.06*r_R # meters
#        twist_distribution = 50*(r_R)+V.Pblade # degrees
#    
    if not Optimise:
        Geo_filename = 'GeoOptimal'+case+'.dat'
        Geo = np.loadtxt(Geo_filename , dtype='float',skiprows = 1)
        R_cp = Geo[:,0]
        twist_cp = Geo[:,1]
        chord_cp = Geo[:,2]

        #Interpolate to find chord and twist distributions
        s_chord = InterpolatedUnivariateSpline(R_cp, chord_cp, k=order)
        chord_distribution = np.array(s_chord(r_R)).reshape((N_ann_reg+1,))
        s_twist = InterpolatedUnivariateSpline(R_cp, twist_cp, k=order)
        twist_distribution = np.array(s_twist(r_R)).reshape((N_ann_reg+1,)) 
        
    data1=pd.read_excel(airfoil, sheet)
    polar_alpha = data1['Alfa']
    polar_cl = data1['Cl']
    polar_cd = data1['Cd']
    
    if Optimise:
        results = np.zeros((N_ann_reg,9))
    else:
        results = np.zeros((N_ann_reg,16))  
    r_R_centroid = np.zeros((N_ann_reg,6))                                #***CHANGE***
    for annulus in range(0,N_ann_reg):                                    #***CHANGE***
        r1_R = r_R[annulus]                                           #***CHANGE***
        r2_R = r_R[annulus+1]                                         #***CHANGE***
        r_R_centroid[annulus] = (r1_R + r2_R)/2                       #***CHANGE***    
        if Optimise and case == 'turbine':
            Area = np.pi*((r2_R*V.R)**2-(r1_R*V.R)**2) #  area streamtube
            fnorm_opti = V.CT_opt*0.5*V.rho*V.U0**2*Area/V.Nblades/V.R/(r2_R-r1_R)
            results[annulus,:] = Opt(case, fnorm_opti, V.U0, r1_R, r2_R, V.R, V.rootradius_R, V.tipradius_R, Omega, V.Nblades, V.rho, polar_alpha, polar_cl, polar_cd, V.TSR)
        elif Optimise and case == 'propeller':
            results[annulus,2] = (r1_R+r2_R)/2
            results[annulus,7] = 50*(r1_R+r2_R)/2+V.Pblade
            results[annulus,8] = 0.18-0.06*(r1_R+r2_R)/2
        else:
            chord = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R, chord_distribution)
            twist = np.interp((r_R[annulus]+r_R[annulus+1])/2, r_R, twist_distribution)
            results[annulus,:] = SS(case, V.rho, V.U0, r1_R, r2_R, V.rootradius_R, V.tipradius_R, Omega, V.R, V.Nblades, chord, twist, polar_alpha, polar_cl, polar_cd, V.TSR)
        
    #Plot results
    if not Optimise or case == 'turbine':
        dr = (r_R[1:]-r_R[:-1])*V.R
        CT = np.sum(dr*results[:,3]*V.Nblades/(0.5*V.U0**2*V.rho*np.pi*V.R**2))
        CP = np.sum(dr*results[:,4]*results[:,2]*V.Nblades*V.R*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))

    if plot:
        if Optimise and case == 'turbine':
            print('Optimal CT is ', CT)
            print('Optimal CP is ', CP)
            
            #Axial and tangential induction factor over radius (Optimal case)
            fig_a_opt = plt.figure(figsize=(12,6))
            plt.title('Axial and tangential induction factor for optimised turbine')
            plt.plot(results[:,2], results[:,0], 'r-', label=r'$a$')          #***CHANGE***
            plt.plot(results[:,2], results[:,1], 'g--', label=r'$a^,$')       #***CHANGE***
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()
        
            #Twist distribution (Optimal case)
            fig_geo_opt = plt.figure(figsize=(12,6))
            plt.title('Twist distribution for optimised turbine')
            plt.plot(results[:,2], results[:,6]*180/np.pi, 'r-.', label='inflowangle')          #***CHANGE***
            plt.plot(results[:,2], results[:,5], 'k-', label=r'$\alpha_{opt}$')
            plt.plot(results[:,2], results[:,7], 'g--', label='twist')       #***CHANGE***
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()
        
            fig_force_opt = plt.figure(figsize=(12,6))
            plt.title(r'Normal and tangential force for optimised turbine, non-dimensionalised by $\frac{1}{2} \rho U_\infty^2 R$')
            plt.plot(results[:,2], results[:,3]/(0.5*V.U0**2*V.rho*V.R), 'k-', label=r'$\alpha_{opt}$')
            plt.plot(results[:,2], results[:,4]/(0.5*V.U0**2*V.rho*V.R), 'g--', label='optimal')       #***CHANGE***
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()
        
            fig_chord_opt = plt.figure(figsize=(12,6))
            plt.title(r'Chord distribution for optimised turbine')
            plt.plot(results[:,2], results[:,8], 'r-', label='chord')       #***CHANGE***
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()   
            
            fig_circ = plt.figure(figsize=(12,6))
            plt.title(r'Circulation distribution for optimised turbine, non-dimensionalised by $\frac{\pi U_\infty^2}{\Omega * NBlades}$')
            plt.plot(results[:,2], results[:,5]/(np.pi*V.U0**2/(V.Nblades*Omega)), 'r-', label=r'$\Gamma$')#***CHANGE***
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()
            
        elif not Optimise: 
            CT = np.sum(dr*results[:,3]*V.Nblades/(0.5*V.U0**2*V.rho*np.pi*V.R**2))
            CP = np.sum(dr*results[:,4]*results[:,2]*V.Nblades*V.R*Omega/(0.5*V.U0**3*V.rho*np.pi*V.R**2))
            print('CT is ', CT)
            print('CP is ', CP)
            
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
            
        #-----------------------------------------------------------------------------

    #--------------------------------------------------------------
    #----------Save chord and twist distribution-------------------
    #--------------------------------------------------------------
    if Optimise:
        if SaveGeo:
            Head = 'X    Twist    Chord'
            file_path = 'GeoOptimal'+case+'.dat'
            SaveMat = results[:,[2,7,8]]
            np.savetxt(file_path, SaveMat, fmt='%f', header=Head)            
    
    return results
#BEM_Model(20, 'constant', 'propeller', show_results=True)