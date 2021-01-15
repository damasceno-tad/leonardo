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
import matplotlib.pyplot as plt
import ISA_leo as isa 
import leo_sim as ls

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


t_f = np.arange(0, 120+0.1, 0.1)
mass_init = launcher_propl_mass + launcher_empty_mass
inp_data = [mass_init, motor_thrust, motor_rate, launcher_surface,
            launcher_dragcoeff, launcher_empty_mass ]


fdata2 = ls.sim(t_f, inp_data)


plt.plot(t_f, fdata2[:,0], 'b')
plt.grid(True)
# plt.axis([90, 105, -5, 5])

plt.figure()
plt.plot(t_f, fdata2[:,1], 'r')
# plt.axis([90, 105, -5, 5])
plt.grid(True)

plt.figure()
plt.plot(t_f, fdata2[:,2], 'k')
# plt.axis([90, 105, -5, 5])
plt.grid(True)
  





