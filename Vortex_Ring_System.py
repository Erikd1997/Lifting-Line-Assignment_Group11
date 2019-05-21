#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 12:00:17 2019

@author: TannerRusk
"""

import numpy as np

def induced_velocity_vortex_system(BigMatrix, controlpoints, gamma, Omega):
    NBlades = len(BigMatrix[0,0,:,0])
    Ncp_Blade = len(BigMatrix[:,0,0,0])
    Ncp = Ncp_Blade * NBlades
    Umatrix = np.zeros((Ncp,Ncp))
    Vmatrix = np.zeros((Ncp,Ncp))
    Wmatrix = np.zeros((Ncp,Ncp))
    
    for Blade_cp in range(NBlades):
        theta0 = Blade_cp* 2*np.pi/NBlades
        for icp_Blade in range(Ncp_Blade):
            for Blade_ring in range(NBlades):
                for jring_Blade in range(Ncp_Blade):
                    icp = Blade_cp * Ncp_Blade + icp_Blade
                    jring = Blade_ring * Ncp_Blade + jring_Blade
                    X1 = BigMatrix[jring_Blade,:-1,Blade_ring,:].reshape((len(BigMatrix[0,1:,0,0]),3))
                    X2 = BigMatrix[jring_Blade,1:,Blade_ring,:].reshape((len(BigMatrix[0,1:,0,0]),3)) 
                    Xp = [0, controlpoints[icp_Blade]*np.sin(theta0), controlpoints[icp_Blade]*np.cos(theta0)]
                    [U, V, W] = velocity_matrix_from_vortex_filament(X1, X2, Xp, gamma, Omega, theta0)
                    Umatrix[icp,jring] = U
                    Vmatrix[icp,jring] = V
                    Wmatrix[icp,jring] = W
    return [np.matrix(Umatrix), np.matrix(Vmatrix), np.matrix(Wmatrix)]




def velocity_matrix_from_vortex_filament(X1, X2, Xp, gamma, omega, theta0):
    R1 = np.sqrt((Xp[0]-X1[:,0])**2+(Xp[1]-X1[:,1])**2+(Xp[2]-X1[:,2])**2)
    R2 = np.sqrt((Xp[0]-X2[:,0])**2+(Xp[1]-X2[:,1])**2+(Xp[2]-X2[:,2])**2)
    R_12_X = (Xp[1]-X1[:,1])*(Xp[2]-X2[:,2])-(Xp[2]-X1[:,2])*(Xp[1]-X2[:,1])
    R_12_Y = -(Xp[0]-X1[:,0])*(Xp[2]-X2[:,2]) + (Xp[2]-X1[:,2])*(Xp[0]-X2[:,0])
    R_12_Z = (Xp[0]-X1[:,0])*(Xp[1]-X2[:,1]) - (Xp[1]-X1[:,1])*(Xp[0]-X2[:,0])
    R_12_sqr = R_12_X**2 + R_12_Y**2 + R_12_Z**2
    R0_1 = (X2[:,0] - X1[:,0])*(Xp[0] - X1[:,0]) + (X2[:,1] - X1[:,1])*(Xp[1] - X1[:,1]) + (X2[:,2] - X1[:,2])*(Xp[2] - X1[:,2])
    R0_2 = (X2[:,0] - X1[:,0])*(Xp[0] - X2[:,0]) + (X2[:,1] - X1[:,1])*(Xp[1] - X2[:,1]) + (X2[:,2] - X1[:,2])*(Xp[2] - X2[:,2])
    
    #Use solid body rotation for velocity if points lies on the vortex
    #(Could only occur for one vortex filament)
    j = 0       #Index to keep track of deleted elements
    U_add = 0
    V_add = 0
    W_add = 0
    for k in range(len(R_12_sqr)):
        i = k-j
        if np.round(R_12_sqr[i],6) == 0:
            j += 1
            R_12_sqr = np.delete(R_12_sqr, i)
            R0_1 = np.delete(R0_1, i)
            R1 =np.delete(R1, i)
            R0_2 =np.delete(R0_2, i)
            R2 =np.delete(R2, i)
            R_12_X =np.delete(R_12_X, i)
            R_12_Y =np.delete(R_12_Y, i)
            R_12_Z =np.delete(R_12_Z, i)
            if not isinstance(gamma, int):
                del gamma[i]
            U_add = 0
            V_add = 0*omega*np.sqrt(Xp[0]**2+Xp[1]**2+Xp[2]**2)*np.cos(theta0)
            W_add = 0*omega*np.sqrt(Xp[0]**2+Xp[1]**2+Xp[2]**2)*np.sin(theta0)
    
    K = gamma / (4*np.pi*R_12_sqr)*(R0_1/R1 - R0_2/R2)
    U = K*R_12_X + U_add
    V = K*R_12_Y + V_add
    W = K*R_12_Z + W_add
    return [sum(U), sum(V), sum(W)]