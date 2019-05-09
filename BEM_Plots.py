# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:28:06 2019

@author: ErikD
"""
import matplotlib.pyplot as plt
import numpy as np

#Results for thrust and power factor
areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*Radius**2
dr = (r_R[1:]-r_R[:-1])*Radius
CT = np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2))
CP = np.sum(dr*results[:,4]*results[:,2]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))

print(r'C_{t} is ', CT)
print(r'C_{p} is ', CP)
#Axial and tangential induction factor over radius
fig_a = plt.figure(figsize=(12,6))
plt.title('Axial and tangential induction factor')
plt.plot(results[:,2], results[:,0], 'r-', label=r'$a$')
plt.plot(results[:,2], results[:,1], 'g--', label=r'$a^,$')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Normal and tangential forces
fig_force = plt.figure(figsize=(12,6))
plt.title(r'Normal and tangential force, non-dimensionalised by $\frac{1}{2} \rho U_\infty^2 R$')
plt.plot(results[:,2], results[:,3]/(0.5*Uinf**2*Radius), 'r-', label=r'Fnorm')
plt.plot(results[:,2], results[:,4]/(0.5*Uinf**2*Radius), 'g--', label=r'Ftan')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Circulation distribution
fig_circ = plt.figure(figsize=(12,6))
plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega * NBlades } $')
plt.plot(results[:,2], results[:,5]/(np.pi*Uinf**2/(NBlades*Omega)), 'r-', label=r'$\Gamma$')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Tip correction influence
fig_tip = plt.figure(figsize=(12,6))
plt.title('Influence of tip correction')
plt.plot(XXX, YYY, 'r-', label=r'C_{t}')
plt.plot(XXX, ZZZ, 'g--', label=r'C_{p}')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Number of annuli influence
fig_annuli = plt.figure(figsize=(12,6))
plt.title('Influence of number of annuli')
plt.plot(XXX, YYY, 'r-', label=r'C_{t}')
plt.plot(XXX, ZZZ, 'g--', label=r'C_{p}')
plt.grid()
plt.xlabel('# of annuli')
plt.legend()
plt.show()

#Convergence history
fig_convergence = plt.figure(figsize=(12,6))
plt.title('Convergence history of total thrust')
plt.plot(XXX, YYY, 'r-', label=r'T_{total}')
plt.grid()
plt.xlabel('# of in-loop iterations')
plt.legend()
plt.show()

#Stagnation enthalpy at four locations
fig_enthalpy = plt.figure(figsize=(12,6))
plt.title('Stagnation enthalpy as function of radius')
plt.plot(XXX, YYY, label='infinity upwind')
plt.plot(XXX, ZZZ, label='at rotor (upwind side)')
plt.plot(XXX, AAA, label='at rotor (downwind side)')
plt.plot(XXX, BBB, label='infinity downwind')
plt.grid()
plt.xlabel('r/R')
plt.legend()
plt.show()

#Maximising Cp or efficiency
fig_optimisation = plt.figure(figsize=(12,6))
plt.title(r'Effect of optimisation on C_{p} or efficiency')
plt.plot(XXX, YYY, '', label=r'C_{p}')
plt.plot(XXX, ZZZ, '', label='efficiency')
plt.grid()
plt.xlabel('# of optimisation iterations')
plt.legend()
plt.show()