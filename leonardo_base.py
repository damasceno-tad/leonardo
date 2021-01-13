# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

print('Bora men')

# My first goal is to create a 1D simulation to evaluate the final velocity
# at 100 km altitude

import math
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

pi = math.pi

# Environment
gravity = 9.81      # m/s2
air_density = 1.0   # kg/m3 #TO DO: function of H, using ISA

# Launcher          # TODO: make this an objec
launcher_empty_mass = 10000.0   # kg
launcher_propl_mass = 25000.0   # kg
launcher_diameter = 2.0         # m
launcher_surface = 0.25 * pi * pow(launcher_diameter, 2.0) # m2
launcher_dragcoeff = 0.4 # nodim # TODO: function of Mach


# Propulsion        # TODO: make this an obj
motor_rate = 250.0   # kg/s
motor_Isp = 350.0    # s
motor_thrust = gravity * motor_Isp * motor_rate     # N


# Satellite
satellite_mass = 10.0 # kg

## TODO: create simulation with odeint

def dmass(y, t, mass_rate):
    if y > 0:
        dmdt = -mass_rate
    else:
        dmdt = 0
        
    return dmdt

t = np.arange(0, 105+1, 1)

mass_init = launcher_propl_mass
mass = odeint(dmass, mass_init, t, args = (motor_rate,))
# plt.plot(t, mass)
# plt.axis([90, 105, -5, 5])
# plt.grid(True)
# print(mass.min())


# Simplified simulation

def drag(x_vel, QnoVxS, d_coeff):
    D = QnoVxS * pow(x_vel, 2.0) * d_coeff
    return D
    

def flight(y, t, data):
    alt =  y[0]
    vel =  y[1]
    mass = y[2]
    
    thrust    = data[0]
    mass_rate = data[1]
    qs        = data[2]    
    cd        = data[3]
    g         = data[4] 
    
    drag_force = drag(vel, qs, cd)
    accel = (thrust - drag_force)/mass - g
    # TODO verify ascending or descending movement to correct drag orientation
    
    if mass > 0:
        dmdt = -mass_rate
    else:
        dmdt = 0
        
    dHdt = vel
    dvdt = accel
        
    return [dHdt, dvdt, dmdt]


t_f = np.arange(0, 40+0.1, 0.1)

qnvs = 0.5 * air_density * launcher_surface
fdata = odeint(flight, [0,0,mass_init], t_f,
               args = ([motor_thrust, motor_rate,
                        qnvs, launcher_dragcoeff, gravity],))

plt.plot(t_f, fdata[:,0], 'b')
plt.grid(True)
# plt.axis([90, 105, -5, 5])

plt.figure()
plt.plot(t_f, fdata[:,1], 'r')
# plt.axis([90, 105, -5, 5])
plt.grid(True)

  





