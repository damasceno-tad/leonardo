# My first goal is to create a 1D simulation to evaluate the final velocity
# at 100 km altitude

from numpy import pi, arange
import matplotlib.pyplot as plt
import ISA_leo as isa 
import leo_sim as ls

# Apogee
apogee = 100e3

# Environment
isa.atmosphere_setup()
latitude = 5 * pi / 180
gravity = isa.gravity_sea_level(latitude)
print(gravity)
# gravity = 9.81      # m/s2
# air_density = 1.0   # kg/m3 #TO DO: function of H, using ISA

# Launcher
launcher = {}

## Propulsion
motor = {}
motor['rate'] = 250.0   # kg/s
motor['isp'] = 350.0    # s
motor['thrust'] = gravity * motor['isp'] * motor['rate']     # N
motor['propl_mass'] = 25000.0   # kg

## Geometry
geometry = {}
geometry['empty_mass'] = 10000.0   # kg
geometry['diameter'] = 2.0         # m
geometry['surface'] = 0.25 * pi * geometry['diameter'] ** 2.0 # m2
geometry['dragcoeff'] = 0.4 # nodim # TODO: function of Mach

launcher['motor'] = motor
launcher['geometry'] = geometry

# Satellite
satellite = {}
satellite['mass'] = 10.0 # kg


t_f = arange(0, 120+0.1, 0.1)
mass_init = launcher['motor']['propl_mass'] + launcher['geometry']['empty_mass']
inp_data = [mass_init, launcher['motor']['thrust'], launcher['motor']['rate'], 
            launcher['geometry']['surface'], launcher['geometry']['dragcoeff'],
            launcher['geometry']['empty_mass']]


fdata2 = ls.sim(t_f, inp_data)


# Plots
plt.figure(1)
plt.plot(t_f, fdata2[:, 0]*1e-3, 'b')
plt.xlabel('Time [s]')
plt.ylabel('Altitude [km]')
plt.grid(True)
plt.show()
# plt.axis([90, 105, -5, 5])

plt.figure(2)
plt.plot(t_f, fdata2[:,1]*1e-3, 'r')
plt.xlabel('Time [s]')
plt.ylabel('Velocity [km/s]')
plt.grid(True)
plt.show()
# plt.axis([90, 105, -5, 5])


plt.figure(3)
plt.plot(t_f, fdata2[:,2]*1e-3, 'k')
plt.xlabel('Time [s]')
plt.ylabel('Mass [1000 kg]')
plt.grid(True)
plt.show()
# plt.axis([90, 105, -5, 5])
