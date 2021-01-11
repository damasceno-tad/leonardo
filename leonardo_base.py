# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

print('Bora men')

drag_coeff = 0.4

print(drag_coeff)


# My first goal is to create a 1D simulation to evaluate the final velocity
# at 100 km altitude

import math
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



