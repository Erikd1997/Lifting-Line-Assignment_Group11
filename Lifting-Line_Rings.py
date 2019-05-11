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
    r_U = R[1:].reshape([len(R)-1, 1, 1, 1])
    r_B = R[:-1].reshape([len(R)-1, 1, 1, 1])
    chord_U = chord[1:].reshape([len(R)-1, 1, 1, 1])
    chord_B = chord[:-1].reshape([len(R)-1, 1, 1, 1])
    
    #Blade number
    Blades = np.arange(NBlades)
    Blades = Blades.reshape([1, 1, NBlades, 1])
    theta_0 = (Blades) * 2*np.pi/NBlades
    
    #Setup time array with extra zeros in middle
    t_end = n_r*2*np.pi/omega
    t_int = np.linspace(0,t_end,n_t)
    t_int = list(t_int)
    t_int.insert(0, 0)
    t = np.array(t_int)
    t = t.reshape([1, n_t, 1, 1])
    
    #Calculate 4-dimensional matrix
    zero_offset = np.zeros([len(chord_U), 1, 1, 1])
    
    offset_U_int = n_t*[chord_U]
    offset_U = np.array(offset_U_int)
    offset_U = offset_U.reshape([len(chord_U), n_t, 1, 1])
    offset_U = np.concatenate((zero_offset,offset_U))
    
    offset_B_int = n_t*[chord_B]
    offset_B = np.array(offset_B_int)
    offset_B = offset_B.reshape([len(chord_B), n_t, 1, 1])
    offset_B = np.concatenate((zero_offset,offset_B))
    
    xw_U = t*U_w + offset_U
    yw_U = -r_U*np.sin(omega*t+theta_0)
    zw_U = r_U*np.cos(omega*t+theta_0)
    
    xw_B = np.fliplr(t*U_w + offset_B)
    yw_B = np.fliplr(-r_B*np.sin(omega*t+theta_0))
    zw_B = np.fliplr(r_B*np.cos(omega*t+theta_0))
    
    xw = np.transpose(np.hstack((xw_B, xw_U)))
    yw = np.transpose(np.hstack((yw_B, yw_U)))
    zw = np.transpose(np.hstack((zw_B, zw_U)))
            
    Rings[:,:,:,0] = xw
    Rings[:,:,:,1] = yw
    Rings[:,:,:,2] = zw
    if plot:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for Blade in range(NBlades):
            for Radial_Point in range(len(r_U)):
                ax.plot(Rings[Radial_Point,:,Blade,0], Rings[Radial_Point,:,Blade,1], Rings[Radial_Point,:,Blade,2], 'b-')
                
        plt.show()
    return Rings

#Defined elsewhere
omega = 10
U_infty = 10

#Define input parameters
n_t = 50
n_r = 10
a_w = 0.2

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