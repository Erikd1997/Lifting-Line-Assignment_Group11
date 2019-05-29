# -*- coding: utf-8 -*-
"""
Created on Sun May 19 23:24:01 2019

@author: ErikD
"""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

def loadBladeOverAllElements(rho, U_infty, u, v, Omega, controlpoints, twist, polar_alpha, polar_cl, polar_cd, chord):
    """
    calculates the loads on all blade elements in a single blade
    """
    #First determine velocities
    Vaxial = U_infty - u
    Vtan = Omega*controlpoints + v              #Only take first blade wich
                                                #stands vertical
    Vp = np.sqrt(np.multiply(Vaxial, Vaxial) + np.multiply(Vtan, Vtan))
    
    #Next determine force coefficients on the blade elements
    inflowangle = np.arctan2(Vaxial,Vtan)
    alpha = inflowangle*180/np.pi+twist
    s_cl = InterpolatedUnivariateSpline(polar_alpha, polar_cl, k=1)
    s_cd = InterpolatedUnivariateSpline(polar_alpha, polar_cd, k=1)
    cl = s_cl(alpha)
    cd = s_cd(alpha)
    
    #From the force coefficients determine the lift and drag
    lift = 0.5*rho*np.multiply(np.multiply(np.multiply(Vp,Vp),cl),chord)
    drag = 0.5*rho*np.multiply(np.multiply(np.multiply(Vp,Vp),cd),chord)
    fnorm = np.multiply(lift,np.cos(inflowangle))+np.multiply(drag,np.sin(inflowangle))
    ftan = np.multiply(lift,np.sin(inflowangle))-np.multiply(drag,np.cos(inflowangle))
    
    return [fnorm, ftan, alpha, inflowangle, Vaxial, Vtan]