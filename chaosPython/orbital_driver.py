#python code to store the orbital driver function

import numpy as np 
# import scipy.integrate as integrate

from .quaternions import body_to_inertial
# from integration_wrapper import ivp_solver
import os
# from transformations import kep_to_eq, equi_to_r_cart, equi_to_v_cart


def time_out(thrust_nominal, firing_nominal, I):
    total = thrust_nominal * firing_nominal
    left = total - I
    if left > 0:
        return 1
    else:
        return 0
   
def modified_equinoctial(t, satellite_state, thrust, orbital_perturb, satellite, grid, environment):
    #define constants:
    equinoctial = satellite_state[:6]
    quaternion = satellite_state[6:10]
    mu = environment.mu
    #update equinoctial element in satellite class 
    satellite.modEquinoctialElements = equinoctial
    p, f, g, h, k, L = equinoctial

    # print('ODE, force:', orbital_perturb)
    # define w, s, alpha
    w = 1 + f * np.cos(L) + g * np.sin(L)
    s = 1 + (h**2) + (k**2)

    
    #get ECI state vector
    r_ECI = satellite.eci_position()

    #compute velocity vector
    v_ECI = satellite.eci_velocity()

    #rotate thrust in inertial direction

    #temporary: constant force thrust in velocity direction:
    # thrust_dir = v_ECI / np.linalg.norm(v_ECI)
    # const_thrust_vec = -((thrust  / satellite.mass) / 1000) * thrust_dir

    #actual acceleration vector
    #flip thrust vector when grid changes
    thrust = np.linalg.norm(thrust) * grid.thrustAxis
    # print('thrust', thrust)
    acc = (thrust / satellite.mass)/1000 # ac, in km/s^2
    acc_in = body_to_inertial(acc, quaternion)#
    # coss = np.dot(acc_in, v_ECI)/ (np.linalg.norm(acc_in) * np.linalg.norm(v_ECI))
    # dir = np.arccos(coss) 

    # file = open('run0/test' + '.txt', 'a')
    # file.write('{}\t{}\t{}\n'.format(t, environment.percentShadow(t), 0))
    # file.close()
    # orbital_perturb, torque_perturb = satellite.envPerturbations(t)#
    #compute vec_a_rsw
    #perturbation in cartesian, km/s^2
    a_perturb_ECI = orbital_perturb +acc_in #+ const_thrust_vec * time_out(t)# #egm08(t, v_ECI, r_ECI)##drag(r_ECI, 1100, v_ECI, 2.2, 1, 100)## 
    # print('ODE, acc', a_perturb_ECI)
    #compute rsw axis 
    ur = r_ECI / np.linalg.norm(r_ECI) # radial direction

    cross = np.cross(r_ECI, v_ECI)
    un = ( cross/ np.linalg.norm(cross)) # normal direction
    
    cross2 = np.cross(cross, r_ECI)
    ut = (cross2 / (np.linalg.norm(cross) * np.linalg.norm(r_ECI))) # tangential direction

    RSW_matrix = np.array([[ur[0], ut[0], un[0]], [ur[1], ut[1], un[1]], [ur[2], ut[2], un[2]]])

    RSW_inverse = np.linalg.inv(RSW_matrix)

    vec_a_rsw = np.matmul(RSW_inverse, a_perturb_ECI)


    dp = np.sqrt(p/mu) *((2*p)/(w)) *vec_a_rsw[1]

    df = np.sqrt(p/mu) *(vec_a_rsw[0] * np.sin(L) + ( (w+1) * np.cos(L) + f ) *( vec_a_rsw[1]/w ) - (h * np.sin(L) - k * np.cos(L)) * (g/w) * (vec_a_rsw[2]))

    dg = np.sqrt(p/mu) *(-vec_a_rsw[0] * np.cos(L) + ( (w+1) * np.sin(L) + g ) * (vec_a_rsw[1]/w) + (h*np.sin(L) - k*np.cos(L))*(f/w) * vec_a_rsw[2])

    dh = np.sqrt(p/mu) *((s) * np.cos(L)/(2*w)) * vec_a_rsw[2]

    dk = np.sqrt(p/mu) *((s) * np.sin(L)/(2*w)) * vec_a_rsw[2]

    dL = np.sqrt(mu * p) * ((w/p)**2) + (1/w) * np.sqrt(p/mu) * (h * np.sin(L) - k *np.cos(L)) * vec_a_rsw[2]
    
    osculating = np.array([dp, df, dg, dh, dk, dL])
    # print('ODE, oscualting', osculating)
    # print('simulated days', t/86400)
    
    
    return osculating


# def stdalone_modified_equinoctial(t, satellite_state, thrust_nominal, satellite, grid, environment, alpha_fun, firing_nominal):
#     #define constants:
#     equinoctial = satellite_state[:8]
#     # quaternion = satellite_state[6:10]
#     mu = environment.mu
#     #update equinoctial element in satellite class 
#     satellite.modEquinoctialElements = equinoctial[:-2]
#     p, f, g, h, k, L, I, I_v= equinoctial

    
#     #define w, s, alpha
#     w = 1 + f * np.cos(L) + g * np.sin(L)
#     s = 1 + (h**2) + (k**2)

    
#     #get ECI state vector
#     r_ECI = satellite.eci_position()

#     #compute velocity vector
#     v_ECI = satellite.eci_velocity()


#     #get the altitude:
#     altitude = np.linalg.norm(r_ECI) - 6378

#     #get angle at altitude, based on function
#     alpha = alpha_fun(altitude) * np.pi/180

#     # get effective thrust:
#     tau_eff = (1 - np.cos(alpha)**2) /4
#     tau_t = (1 - np.cos(alpha)) /2
#     thrust = thrust_nominal * tau_eff

#     #rotate thrust in inertial direction

#     #temporary: constant force thrust in velocity direction:
#     thrust_dir = v_ECI / np.linalg.norm(v_ECI)
#     const_thrust_vec = -((thrust  / satellite.mass) / 1000) * thrust_dir



#     # file = open('run0/test' + '.txt', 'a')
#     # file.write('{}\t{}\t{}\n'.format(t, environment.percentShadow(t), 0))
#     # file.close()
#     orbital_perturb, torque_perturb = satellite.envPerturbations(t)#
    
#     #compute vec_a_rsw
#     #perturbation in cartesian, km/s^2
#     a_perturb_ECI = orbital_perturb + const_thrust_vec * time_out(thrust_nominal, firing_nominal, I)# +acc_in#egm08(t, v_ECI, r_ECI)##drag(r_ECI, 1100, v_ECI, 2.2, 1, 100)## 

#     #compute rsw axis 
#     ur = r_ECI / np.linalg.norm(r_ECI) # radial direction

#     cross = np.cross(r_ECI, v_ECI)
#     un = ( cross/ np.linalg.norm(cross)) # normal direction
    
#     cross2 = np.cross(cross, r_ECI)
#     ut = (cross2 / (np.linalg.norm(cross) * np.linalg.norm(r_ECI))) # tangential direction

#     RSW_matrix = np.array([[ur[0], ut[0], un[0]], [ur[1], ut[1], un[1]], [ur[2], ut[2], un[2]]])

#     RSW_inverse = np.linalg.inv(RSW_matrix)

#     vec_a_rsw = np.matmul(RSW_inverse, a_perturb_ECI)


#     dp = np.sqrt(p/mu) *((2*p)/(w)) *vec_a_rsw[1]

#     df = np.sqrt(p/mu) *(vec_a_rsw[0] * np.sin(L) + ( (w+1) * np.cos(L) + f ) *( vec_a_rsw[1]/w ) - (h * np.sin(L) - k * np.cos(L)) * (g/w) * (vec_a_rsw[2]))

#     dg = np.sqrt(p/mu) *(-vec_a_rsw[0] * np.cos(L) + ( (w+1) * np.sin(L) + g ) * (vec_a_rsw[1]/w) + (h*np.sin(L) - k*np.cos(L))*(f/w) * vec_a_rsw[2])

#     dh = np.sqrt(p/mu) *((s) * np.cos(L)/(2*w)) * vec_a_rsw[2]

#     dk = np.sqrt(p/mu) *((s) * np.sin(L)/(2*w)) * vec_a_rsw[2]

#     dL = np.sqrt(mu * p) * ((w/p)**2) + (1/w) * np.sqrt(p/mu) * (h * np.sin(L) - k *np.cos(L)) * vec_a_rsw[2]

#     #impulse consumed:
#     dI = thrust_nominal * tau_t * time_out(thrust_nominal, firing_nominal, I)
#     #impulse against velocity
#     dI_v = thrust_nominal * tau_eff *time_out(thrust_nominal, firing_nominal, I)


#     osculating = np.array([dp, df, dg, dh, dk, dL, dI, dI_v])

    

#     # print('simulated days', t/86400)
#     return osculating










