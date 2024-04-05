"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


orbital_driver.py

This Python code defines the function used to solve the orbital mechanics 
equations of motion for a spacecraft. The function is an ODE  
 that considers the effects of thrust and orbital perturbations.

The function takes the following arguments:

- t (float): Current simulation time (seconds).
- satellite_state (numpy.ndarray): Array containing the current orbital state 
  of the spacecraft (equinoctial elements and quaternion).
- thrust (numpy.ndarray): Thrust vector acting on the spacecraft (in body frame).
- orbital_perturb (numpy.ndarray): Perturbation forces acting on the spacecraft 
  in the ECI frame (e.g., drag, atmospheric effects).
- satellite (class object): Instance of the `Satellite` class representing the 
  spacecraft properties.
- grid (class object): Instance of the `Grid` class representing the environment 
  interaction grid.
- environment (class object): Instance of the `Environment` class representing 
  the environmental conditions.

The function returns the time derivatives of the equinoctial elements 
(rates of change).
"""
import numpy as np 

from .quaternions import body_to_inertial



def modified_equinoctial(t, satellite_state, thrust, orbital_perturb, satellite, grid, environment):
    
    #split state
    equinoctial = satellite_state[:6]
    quaternion = satellite_state[6:10]
    mu = environment.mu

    #update equinoctial element in satellite class 
    satellite.modEquinoctialElements = equinoctial
    p, f, g, h, k, L = equinoctial


    # define w, s
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
    acc = (thrust / satellite.mass)/1000 # ac, in km/s^2
    acc_in = body_to_inertial(acc, quaternion)#



    #compute acc in rsw frame

    #perturbation in cartesian, km/s^2
    a_perturb_ECI = orbital_perturb +acc_in 
    #compute rsw axis 
    ur = r_ECI / np.linalg.norm(r_ECI) # radial direction
    cross = np.cross(r_ECI, v_ECI)
    un = ( cross/ np.linalg.norm(cross)) # normal direction
    cross2 = np.cross(cross, r_ECI)
    ut = (cross2 / (np.linalg.norm(cross) * np.linalg.norm(r_ECI))) # tangential direction

    RSW_matrix = np.array([[ur[0], ut[0], un[0]], [ur[1], ut[1], un[1]], [ur[2], ut[2], un[2]]])

    RSW_inverse = np.linalg.inv(RSW_matrix)

    vec_a_rsw = np.matmul(RSW_inverse, a_perturb_ECI)




    #osculating elements
    dp = np.sqrt(p/mu) *((2*p)/(w)) *vec_a_rsw[1]

    df = np.sqrt(p/mu) *(vec_a_rsw[0] * np.sin(L) + ( (w+1) * np.cos(L) + f ) *( vec_a_rsw[1]/w ) - (h * np.sin(L) - k * np.cos(L)) * (g/w) * (vec_a_rsw[2]))

    dg = np.sqrt(p/mu) *(-vec_a_rsw[0] * np.cos(L) + ( (w+1) * np.sin(L) + g ) * (vec_a_rsw[1]/w) + (h*np.sin(L) - k*np.cos(L))*(f/w) * vec_a_rsw[2])

    dh = np.sqrt(p/mu) *((s) * np.cos(L)/(2*w)) * vec_a_rsw[2]

    dk = np.sqrt(p/mu) *((s) * np.sin(L)/(2*w)) * vec_a_rsw[2]

    dL = np.sqrt(mu * p) * ((w/p)**2) + (1/w) * np.sqrt(p/mu) * (h * np.sin(L) - k *np.cos(L)) * vec_a_rsw[2]
    
    osculating = np.array([dp, df, dg, dh, dk, dL])

    
    
    return osculating


