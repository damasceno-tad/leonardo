# leo_sim is the library of a simplified 1D simulation to evaluate altitude
# and velocity for a satellite launcher 

from scipy.integrate import odeint
import ISA_leo as isa 


# The drag function calculates the drag force     
def drag(x_vel, rho, S, d_coeff):
    D = 0.5 * rho * x_vel ** 2.0 * S * d_coeff
    return D


def flight(y, t, data):
    alt =  y[0]
    vel =  y[1]
    mass = y[2]
    
    theor_thrust = data[0]
    mass_rate    = data[1]
    S            = data[2]    
    cd           = data[3]
    dry_mass     = data[4]

    g = 9.81
    
    # Propulsion verification: stop when propellant mass is finished
    if mass > dry_mass:
        dmdt = -mass_rate
        thrust = theor_thrust
    else:
        dmdt = 0
        thrust = 0
        
    isa.atmosphere_setup()
    T = isa.get_temperature(alt)
    p = isa.get_pressure(alt)
    rho = isa.density(T, p, isa.R)
    drag_force = drag(vel, rho, S, cd)
    accel = (thrust - drag_force) / mass - g
    # TODO verify ascending or descending movement to correct drag orientation
       
    dHdt = vel
    dvdt = accel
        
    return [dHdt, dvdt, dmdt]


def sim(t, inp_data):
    mass_init           = inp_data[0]
    motor_thrust        = inp_data[1]
    motor_rate          = inp_data[2]
    launcher_surface    = inp_data[3]
    launcher_dragcoeff  = inp_data[4]
    launcher_empty_mass = inp_data[5]
    
    fdata = odeint (  flight, [0, 0, mass_init], t,
            args = ( [motor_thrust, motor_rate,
            launcher_surface, launcher_dragcoeff, launcher_empty_mass],)  )
    
    return fdata
