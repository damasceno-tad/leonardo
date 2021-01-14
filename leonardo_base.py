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
import ISA_leo as isa 

pi = math.pi

# Environment
gravity = 9.81      # m/s2
air_density = 1.0   # kg/m3 #TO DO: function of H, using ISA
isa.atmosphere_setup()

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

# Simplified simulation

def drag(x_vel, rho, S, d_coeff):
    D = 0.5 * rho * pow(x_vel, 2.0) * S * d_coeff
    return D
# Changed the drag function to use rho and surface separately     


def flight(y, t, data):
    alt =  y[0]
    vel =  y[1]
    mass = y[2]
    
    theor_thrust = data[0]
    mass_rate    = data[1]
    S            = data[2]    
    cd           = data[3]
    g            = data[4]
    dry_mass     = data[5]
    
    # Propulsion verification: stop when propellant mass is finished
    if mass > dry_mass:
        dmdt = -mass_rate
        thrust = theor_thrust
    else:
        dmdt = 0
        thrust = 0
   
    T = isa.get_temperature(alt)
    p = isa.get_pressure(alt)
    rho = isa.density(T, p, isa.R)
    drag_force = drag(vel, rho, S, cd)
    accel = (thrust - drag_force)/mass - g
    # TODO verify ascending or descending movement to correct drag orientation
    
       
    dHdt = vel
    dvdt = accel
        
    return [dHdt, dvdt, dmdt]


t_f = np.arange(0, 120+0.1, 0.1)


# TODO use total mass and differentiate propellant vs total mass in simul
mass_init = launcher_propl_mass + launcher_empty_mass

fdata = odeint(flight, [0,0,mass_init], t_f,
               args = ([motor_thrust, motor_rate,
                        launcher_surface, launcher_dragcoeff, gravity,
                        launcher_empty_mass],))

plt.plot(t_f, fdata[:,0], 'b')
plt.grid(True)
# plt.axis([90, 105, -5, 5])

plt.figure()
plt.plot(t_f, fdata[:,1], 'r')
# plt.axis([90, 105, -5, 5])
plt.grid(True)

plt.figure()
plt.plot(t_f, fdata[:,2], 'k')
# plt.axis([90, 105, -5, 5])
plt.grid(True)
  





