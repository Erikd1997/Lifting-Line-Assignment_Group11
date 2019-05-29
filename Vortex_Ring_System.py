#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 12:00:17 2019

@author: TannerRusk
"""

import numpy as np

def induced_velocity_vortex_system(BigMatrix, controlpoints, gamma):
    NBlades = len(BigMatrix[0,0,:,0])
    Nwp = len(BigMatrix[0,:,0,0])
    Ncp_Blade = len(BigMatrix[:,0,0,0])
    Ncp = Ncp_Blade * NBlades
    Umatrix = np.zeros((Ncp,Ncp))
    Vmatrix = np.zeros((Ncp,Ncp))
    Wmatrix = np.zeros((Ncp,Ncp))
    
    #Define Gamma in the case of Gamma = 1
    gamma = np.ones([1,Ncp,1])*gamma

    #Create X1 and X2 $X1[controlpoint,jring,wakepoint]$
    X1_x = BigMatrix[:, :-1, 0, 0].reshape([1,Ncp_Blade,Nwp-1])
    X1_y = BigMatrix[:, :-1, 0, 1].reshape([1,Ncp_Blade,Nwp-1])
    X1_z = BigMatrix[:, :-1, 0, 2].reshape([1,Ncp_Blade,Nwp-1])
    
    X2_x = BigMatrix[:, 1:, 0, 0].reshape([1,Ncp_Blade,Nwp-1])
    X2_y = BigMatrix[:, 1:, 0, 1].reshape([1,Ncp_Blade,Nwp-1])
    X2_z = BigMatrix[:, 1:, 0, 2].reshape([1,Ncp_Blade,Nwp-1])
    
    #Create theta0 (for blade distribution over 2pi radians)
    theta0 = np.ones([Ncp_Blade,1,1])*(0) * 2*np.pi/NBlades
    for i in range(1,NBlades):
        X1_x = np.concatenate((X1_x, BigMatrix[:, :-1, i, 0].reshape([1,Ncp_Blade,Nwp-1])),axis=1)
        X1_y = np.concatenate((X1_y, BigMatrix[:, :-1, i, 1].reshape([1,Ncp_Blade,Nwp-1])),axis=1)
        X1_z = np.concatenate((X1_z, BigMatrix[:, :-1, i, 2].reshape([1,Ncp_Blade,Nwp-1])),axis=1)
        
        X2_x = np.concatenate((X2_x, BigMatrix[:, 1:, i, 0].reshape([1,Ncp_Blade,Nwp-1])),axis=1)    
        X2_y = np.concatenate((X2_y, BigMatrix[:, 1:, i, 1].reshape([1,Ncp_Blade,Nwp-1])),axis=1)    
        X2_z = np.concatenate((X2_z, BigMatrix[:, 1:, i, 2].reshape([1,Ncp_Blade,Nwp-1])),axis=1)    
        
        theta0 = np.concatenate((theta0, np.ones([Ncp_Blade,1,1])*(i)*2*np.pi/NBlades),axis=0)

    #Create Xp $Xp[controlpoint,jring,wakepoint]$
    Xp_x = np.zeros([Ncp,1,1])
    Xp_y = -controlpoints.reshape([Ncp,1,1])*np.sin(theta0)
    Xp_z = controlpoints.reshape([Ncp,1,1])*np.cos(theta0)
    
    [Umatrix, Vmatrix, Wmatrix] = velocity_matrix_from_vortex_filament(X1_x, X1_y, X1_z, X2_x, X2_y, X2_z, Xp_x, Xp_y, Xp_z, gamma)
    return [np.matrix(Umatrix), np.matrix(Vmatrix), np.matrix(Wmatrix)]




def velocity_matrix_from_vortex_filament(X1_x, X1_y, X1_z, X2_x, X2_y, X2_z, Xp_x, Xp_y, Xp_z, gamma):
    R1 = (np.multiply((Xp_x-X1_x),(Xp_x-X1_x))+np.multiply((Xp_y-X1_y),(Xp_y-X1_y))+np.multiply((Xp_z-X1_z),(Xp_z-X1_z)))**0.5 #5seconds
    R2 = (np.multiply((Xp_x-X2_x),(Xp_x-X2_x))+np.multiply((Xp_y-X2_y),(Xp_y-X2_y))+np.multiply((Xp_z-X2_z),(Xp_z-X2_z)))**0.5
    R_12_X = np.multiply((Xp_y-X1_y),(Xp_z-X2_z)) - np.multiply((Xp_z-X1_z),(Xp_y-X2_y))
    R_12_Y = -np.multiply((Xp_x-X1_x),(Xp_z-X2_z)) + np.multiply((Xp_z-X1_z),(Xp_x-X2_x))
    R_12_Z = np.multiply((Xp_x-X1_x),(Xp_y-X2_y)) - np.multiply((Xp_y-X1_y),(Xp_x-X2_x))
    R_12_sqr = np.multiply(R_12_X,R_12_X) + np.multiply(R_12_Y,R_12_Y) + np.multiply(R_12_Z,R_12_Z)
    R0_1 = (X2_x - X1_x)*(Xp_x - X1_x) + (X2_y - X1_y)*(Xp_y - X1_y) + (X2_z - X1_z)*(Xp_z - X1_z)
    R0_2 = (X2_x - X1_x)*(Xp_x - X2_x) + (X2_y - X1_y)*(Xp_y - X2_y) + (X2_z - X1_z)*(Xp_z - X2_z)
    
    #If filament intersects with the controlpoint, it has no effect on the
    #induced velocity. This is represented by setting its distance to infinity
    n = np.argwhere(np.round(R_12_sqr,8) == 0)
    R_12_sqr[n[:,0],n[:,1],n[:,2]] = np.inf #2 seconds
    K = gamma / (4*np.pi*R_12_sqr)*(R0_1/R1 - R0_2/R2)
    U = np.multiply(K,R_12_X)
    V = np.multiply(K,R_12_Y)
    W = np.multiply(K,R_12_Z)
    return [np.sum(U,axis=2), np.sum(V,axis=2), np.sum(W,axis=2)]