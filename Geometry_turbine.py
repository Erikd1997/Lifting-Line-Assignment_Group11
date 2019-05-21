# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 21:22:14 2019

@author: ErikD
"""
# geometry wind turbine
# airfoil is DU 95-W-180 (polar in attachment on brightspace)
R = 50  # in meters
rho = 1.225 # Density in kg/m^3
Nblades = 4  # number of blades
Pblade = 2  # blade pitch in degrees
RotorYaw = 0  # rotor yaw angle in degrees
rootradius_R = 0.2  # root radius
tipradius_R = 1     # Tip radius
U0 = 10  # in m/s
CT_opt = 0.75
TSR = 10  # lambda (6, 8 or 10)??