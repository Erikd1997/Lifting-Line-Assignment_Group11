#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 12:00:17 2019

@author: TannerRusk
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def induced_velocity_vortex_system(BigMatrix, controlpoints, gamma):
    NBlades = len(BigMatrix[0,0,:,0])
    Ncp_Blade = len(BigMatrix[:,0,0,0])
    Ncp = Ncp_Blade * NBlades
    Umatrix = np.zeros((Ncp,Ncp))
    Vmatrix = np.zeros((Ncp,Ncp))
    Wmatrix = np.zeros((Ncp,Ncp))
    
    for Blade in range(NBlades):
        for icp_Blade in range(Ncp_Blade):
            for jring in range(Ncp):
                icp = Blade * Ncp_Blade + icp_Blade
                X1 = BigMatrix[icp_Blade,:-1,Blade,:].reshape((len(BigMatrix[0,:-1,0,0]),3))
                X2 = BigMatrix[icp_Blade,1:,Blade,:].reshape((len(BigMatrix[0,:-1,0,0]),3)) 
                Xp = controlpoints[icp_Blade,:]
                [U, V, W] = velocity_matrix_from_vortex_filament(X1, X2, Xp, gamma)
                Umatrix[icp,jring] = U
                Vmatrix[icp,jring] = V
                Wmatrix[icp,jring] = W
    return [Umatrix, Vmatrix, Wmatrix]




def velocity_matrix_from_vortex_filament(X1, X2, Xp, gamma):
    R1 = np.sqrt((Xp[0]-X1[:,0])^2+(Xp[1]-X1[:,1])^2+(Xp[2]-X1[:,2])^2)
    R2 = np.sqrt((Xp[0]-X2[:,0])^2+(Xp[1]-X2[:,1])^2+(Xp[2]-X2[:,2])^2)
    R_12_X = (Xp[1]-X1[:,1])*(Xp[2]-X2[:,2])-(Xp[2]-X1[:,2])*(Xp[1]-X2[:,1])
    R_12_Y = -(Xp[0]-X1[:,0])*(Xp[2]-X2[:,2]) + (Xp[2]-X1[:,2])*(Xp[0]-X2[:,0])
    R_12_Z = (Xp[0]-X1[:,0])*(Xp[1]-X2[:,1]) - (Xp[1]-X1[:,1])*(Xp[0]-X2[:,0])
    R_12_sqr = R_12_X^2 + R_12_Y^2 + R_12_Z^2
    R0_1 = (X2[:,0] - X1[:,0])*(Xp[0] - X1[:,0]) + (X2[:,1] - X1[:,1])*(Xp[1] - X1[:,1]) + (X2[:,2] - X1[:,2])*(Xp[2] - X1[:,2])
    R0_2 = (X2[:,0] - X1[:,0])*(Xp[0] - X2[:,0]) + (X2[:,1] - X1[:,1])*(Xp[1] - X2[:,1]) + (X2[:,2] - X1[:,2])*(Xp[2] - X2[:,2])
    K = gamma / (4*np.pi*R_12_sqr)*(R0_1/R1 - R0_2/R2)
    U = K*R_12_X
    V = K*R_12_Y
    W = K*R_12_Z
    return [sum(U), sum(V), sum(W)]