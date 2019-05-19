# -*- coding: utf-8 -*-
"""
Created on Fri May 17 13:31:34 2019

@author: ErikD
"""
import numpy as np

def newGamma(rho, U_infty, u, v, w, Omega, controlpoints, BigMatrix, twist, polar_alpha, polar_cl, chord):
    """
    calculates the new circulation distribution along the blade
    """
    Vaxial = U_infty + u
    
    Blades = BigMatrix[0,0,:,0].reshape((1,len(BigMatrix[0,0,:,0])))
    theta_0 = (Blades) * 2*np.pi/(Blades[-1]+1)
    Vtan =  np.mat(np.zeros([len(controlpoints)*len(Blades), 1]))
    for i in range(len(controlpoints)):
        for j in range(len(Blades)):
            n_times_vt = np.cos(theta_0[j])*v[i] + np.sin(theta_0[j])*w[i]
            Vtan[j*len(controlpoints)+i] = Omega*controlpoints[i] + n_times_vt
    
    Vp = np.sqrt(np.multiply(Vaxial, Vaxial) + np.multiply(Vtan, Vtan))
    inflowangle = np.arctan2(Vaxial,Vtan)
    alpha = inflowangle*180/np.pi + twist
    cl = np.interp(alpha, polar_alpha, polar_cl)
    gamma = 0.5*np.multiply(np.multiply(Vp,cl),chord)
    return gamma