# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def WakeGeometry(U_infty, omega, n_t, n_r, a_w, NBlades, R, chord, plot=False):
    #Calculate intermediate parameters
    L = 2*n_t+2
    Rings = np.zeros([len(R)-1, L, NBlades, 3])
    U_w = U_infty*(1-a_w)
    
    #Blade geometry distribution
    r_U = R[1:].reshape([len(R)-1, 1, 1])
    r_B = R[:-1].reshape([len(R)-1, 1, 1])
    chord_U = chord[1:].reshape([len(R)-1, 1, 1])
    chord_B = chord[:-1].reshape([len(R)-1, 1, 1])
    
    #Blade number
    Blades = np.arange(NBlades)
    Blades = Blades.reshape([1, 1, NBlades])
    theta_0 = (Blades) * 2*np.pi/NBlades
    
    #Setup time array with extra zeros in middle
    t_end = n_r*2*np.pi/omega
    t_int = np.linspace(0,t_end,n_t)
    t_int = list(t_int)
    t_int.insert(0, 0)
    t = np.array(t_int)
    t = t.reshape([1, n_t+1, 1])
    
    #Calculate 4-dimensional matrix
    zero_offset = np.zeros([len(chord_U), 1, 1])
    
    offset_U = np.ones([1, n_t, 1])*chord_U.reshape([len(chord_U), 1, 1])
    offset_U = np.concatenate((zero_offset,offset_U),axis=1)
    
    offset_B = np.ones([1, n_t, 1])*np.fliplr(chord_B).reshape([len(chord_B), 1, 1])
    offset_B = np.concatenate((offset_B, zero_offset),axis=1)
    
    xw_U = t*U_w + offset_U
    yw_U = -r_U*np.sin(omega*t+theta_0)
    zw_U = r_U*np.cos(omega*t+theta_0)
    
    xw_B = np.fliplr(t)*U_w + offset_B
    yw_B = -r_B*np.sin(omega*np.fliplr(t)+theta_0)
    zw_B = r_B*np.cos(omega*np.fliplr(t)+theta_0)
    
    xw = np.concatenate((xw_B, xw_U), axis=1)
    yw = np.concatenate((yw_B, yw_U), axis=1)
    zw = np.concatenate((zw_B, zw_U), axis=1)
            
    Rings[:,:,:,0] = xw
    Rings[:,:,:,1] = yw
    Rings[:,:,:,2] = zw
    if plot:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
#        ax.plot(Rings[25,:,0,0], Rings[25,:,0,1], Rings[25,:,0,2], 'b-')

        for Blade in range(NBlades):
            for Radial_Point in range(0,len(r_U),10):
                ax.plot(Rings[Radial_Point,:,Blade,0], Rings[Radial_Point,:,Blade,1], Rings[Radial_Point,:,Blade,2], 'b-')
                
        plt.show()
    return Rings

#Defined elsewhere
omega = 10
U_infty = 10

#Define input parameters
n_t = 200
n_r = 1
a_w = 0

NBlades = 3
filename = 'GeoOptimalturbine.dat'
Mat = np.loadtxt(filename, dtype='float',skiprows = 1)

R = Mat[:,0]
R_int = list(R)
R_int.insert(0, R[0]-(R[5]-R[4])/2)
R_int.insert(len(R_int),R[-1]+(R[5]-R[4])/2)
R = np.array(R_int)

Twist = Mat[:,1]

Chord = Mat[:,2]
Chord_int = list(Chord)
Chord_int.insert(0, Chord[0]-(Chord[1]-Chord[0])/2)
Chord_int.insert(len(Chord_int),Chord[-1]+(Chord[-1]-Chord[-2])/2)
Chord = np.array(Chord_int)

WakeGeometry(U_infty, omega, n_t, n_r, a_w, NBlades, R, Chord, plot=True)