# Basic source: http://www.braeunig.us/space/atmmodel.htm

from numpy import sqrt, exp, sin, pi, arange, array, savetxt
import matplotlib.pyplot as plt
# import export_data

outfile_path = '../outfile.csv'

# Earth atmosphere layers
layers = ['Troposphere', 'Tropopause', 'Stratosphere1',
          'Stratosphere2', 'Stratospause', 'Mesosphere1', 'Mesosphere2']
L = len(layers)
min_altitude, max_altitude = 0, 86e3
layer_h0 = [_ * 1e3 for _ in [0, 11, 20, 32, 47, 51, 71]]
layer_temp_rate = [_ * 1e-3 for _ in [-6.5, 0, 1, 2.8, 0, -2.8, -2]]
layer_T0 = []
layer_p0 = []

# Earth properties
r_earth = 6356766
g0 = 9.80665

# Air (dry ideal gas) properties
gamma = 1.4
molar_mass = 28.96442
R_universal = 8314.32
R = R_universal / molar_mass


def atmosphere_setup():
    T_sl, p_sl = sea_level()
    T0 = T_sl
    p0 = p_sl
    layer_T0.append(T0)
    layer_p0.append(p0)
    for i in range(L - 1):
        h0 = layer_h0[i]
        temp_rate = layer_temp_rate[i]
        h = layer_h0[i + 1]
        p0 = pressure(h, p0, T0, h0, temp_rate)
        T0 = temperature(h, T0, h0, temp_rate)
        layer_T0.append(T0)
        layer_p0.append(p0)


def k2degc(K):
    """ Kelvin to degree Celsius."""
    C = K - 273.15
    return C


def sea_level():
    """ Mean sea level properties."""
    T0 = 288.15
    p0 = 1.01325e5
    return T0, p0


def get_layer(h):
    index = [i for i, h0 in enumerate(layer_h0) if h0 <= h][-1]
    return index


def get_layer_params(layer):
    T0 = layer_T0[layer]
    p0 = layer_p0[layer]
    h0 = layer_h0[layer]
    temp_rate = layer_temp_rate[layer]
    return T0, p0, h0, temp_rate


def get_temperature(h):
    layer = get_layer(h)
    T0, p0, h0, temp_rate = get_layer_params(layer)
    T = temperature(h, T0, h0, temp_rate)
    return T


def get_pressure(h):
    layer = get_layer(h)
    T0, p0, h0, temp_rate = get_layer_params(layer)
    p = pressure(h, p0, T0, h0, temp_rate)
    return p


def temperature(h, T0, h0, temp_rate):
    T = T0 + temp_rate * (h - h0)
    return T


def pressure(h, p0, T0, h0, temp_rate):
    T = temperature(h, T0, h0, temp_rate)
    if temp_rate == 0:  # Isothermal
        p = p0 * exp((-g0 / (R * T)) * (h - h0))
    else:               # Gradient
        p = p0 * (T / T0) ** (-g0 / (R * temp_rate))
    return p


def density(T, p, R):
    rho = p / (R * T)
    return rho


def sound_speed(gamma, R, T):
    a = sqrt(gamma * R * T)
    return a


def gravity_sea_level(latitude):
    g0 = 9.780356 * (1 + 0.0052885 * (sin(latitude) ** 2 -
                                      0.0000059 * (sin(2 * latitude)) ** 2))
    return g0


def geometric_altitude(h):
    z = r_earth * h / (r_earth - h)
    return z


def gravity(z, latitude = 45 * pi / 180):
    g = gravity_sea_level(latitude) * (r_earth / (r_earth + z)) ** 2
    return g 


def main():
    T_sl, p_sl = sea_level()
    rho_sl = density(T_sl, p_sl, R)
    a_sl = sound_speed(gamma, R, T_sl)

    altitudes = arange(max_altitude, step=1000)
    sound_speeds = []
    densities = []
    pressures = []
    temperatures = []
    geometric_altitudes = []
    gravities = []

    for h in altitudes:    
        T = get_temperature(h)
        p = get_pressure(h)
        rho = density(T, p, R)
        a = sound_speed(gamma, R, T)
        z = geometric_altitude(h)
        g = gravity(h)
        temperatures.append(T)
        pressures.append(p)
        densities.append(rho)
        sound_speeds.append(a)
        geometric_altitudes.append(z)
        gravities.append(g)

    params = ['Geopotential Height', 'Geometric Height', 'Temperature', 'Pressure', 'Density', 'Sound Speed', 'Gravity']
    units = ['m\'', 'm', 'K', 'Pa', 'kg/m^3', 'm/s', 'm/s^2']
    labels = [param + ' [' + unit + ']' for param, unit in zip(params, units)]
    data = [altitudes, geometric_altitudes, temperatures, pressures, densities, sound_speeds, gravities]
    data_dict = {}
    for label, datum in zip(labels, data):
        data_dict[label] = datum
    # export_data.as_dict(data_dict, outfile_path)

    fig = plt.figure('sound_speeds', dpi=80)
    plt.plot(sound_speeds, altitudes * 1e-3, 'r-', linewidth=2, label='ISA sound speed')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('Sound speed [m/s]')
    plt.ylabel('Altitude [km]')
    fig.tight_layout()
    plt.show()

    fig = plt.figure('density', dpi=80)
    plt.plot(densities, altitudes * 1e-3, 'r-', linewidth=2, label='ISA Density')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('Density [kg/m^3]')
    plt.ylabel('Altitude [km]')
    fig.tight_layout()
    plt.show()

    fig = plt.figure('Pressure', dpi=80)
    plt.plot(array(pressures) * 1e-3, altitudes * 1e-3, 'r-', linewidth=2, label='ISA Pressure')
    plt.plot(array(layer_p0) * 1e-3, array(layer_h0) * 1e-3, 'bo')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('Pressure [kPa]')
    plt.ylabel('Altitude [km]')
    fig.tight_layout()
    plt.show()

    fig = plt.figure('Temperature', dpi=80)
    plt.plot(temperatures, altitudes * 1e-3, 'r-', linewidth=2, label='ISA Temperature')
    plt.plot(layer_T0, array(layer_h0) * 1e-3, 'bo')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('Temperature [K]')
    plt.ylabel('Altitude [km]')
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    atmosphere_setup()

    print('*** ISA CALCULATOR ***')
    h = float(input('Enter altitude [m]: '))

    layer = get_layer(h)
    print(layers[layer])
    
    T_sl, p_sl = sea_level()
    rho_sl = density(T_sl, p_sl, R)
    a_sl = sound_speed(gamma, R, T_sl)
    
    T = get_temperature(h)
    p = get_pressure(h)
    rho = density(T, p, R)
    a = sound_speed(gamma, R, T)
    z = geometric_altitude(h)

    props = [T, p * 1e-3, rho, a, z]
    names = ['Temperature', 'Pressure', 'Density', 'Sound speed', 'Geometric height']
    units = ['K', 'kPa', 'kg/m^3', 'm/s', 'm']
    sub_props = [k2degc(T), 100 * p / p_sl, 100 * rho / rho_sl, 100 * a / a_sl, z * 1e-3]
    sub_units = ['degC', '% SL', '% SL', '% SL', 'km']

    for prop, name, unit, sub_prop, sub_unit in zip(props, names, units, sub_props, sub_units):
        print((name + ':').rjust(18), '%10.3f' % prop, unit.ljust(8), '[', '%6.2f' % sub_prop, sub_unit.ljust(4), ']')

    main()
