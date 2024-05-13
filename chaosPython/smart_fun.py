#python code to hold smart functions to control the attitude / reduce the angular velocity.

"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


smart_fun.py




This Python code defines functions that interact with the driver ansd integrator.
Most of these functions are designed to work as event functions inside an integrator.

NOTE: The arguments passed to events functions will be (t, y, *args) as defined by solve_ivp.
"""

import numpy as np
from .transformations import *
from .sensors import signalToNoise
from .support_fun import circlesIntersection







def pixel_fuel(t, satellite_state, *args):
    """    
    function to return 0 when the pixel has no more fuel. 
    This even function is terminal
    """

    return satellite_state[13]






def r_stop(t, satellite_state, *args):
    """
    Event function to stop the simulation when the perigee altitude reaches 150 km.

    Args:
        t (float): Current simulation time (ignored in this function).
        satellite_state (numpy.ndarray): Satellite state vector containing orbital elements.
        *args: Additional arguments passed to the function (ignored in this function).

    Returns:
        float:
            - A value less than or equal to zero if the perigee altitude is less than or equal to 150 km,
              indicating the stopping condition is met.
            - A positive value otherwise.
    """



    #function to stop simulation when perigee reaches 150 km.
    #extract elements
    p, f, g, h, k, L = satellite_state[:6]

    #get keplerian elements:
    a, e, inc, RAAN, wp, theta = eq_to_kep(p, f, g, h, k, L)

    #get Cartesian position vector
    r_ECI = equi_to_r_cart(p, f, g, h, k, L)

    #compute perigee (height):
    r_p = a*(1 - e) - 6378

    #check if stop condition is met:
    event = 150 - r_p
    return event


def sensorMeasurement(t, satellite_state, thrust, torque_thruster, satellite, grid, selected_grid, controlSystem, environment, sensor, manager, iteration):

    """
    Simulates sensor measurements of the spacecraft's angular velocity.

    This function models a Micro-Electro-Mechanical System (MEMS) gyroscope
    sensor that measures the spacecraft's angular velocity with noise.
    This is an event function, so that the noise is updated at 
    every step of the integration.


    Args:
        t (float): Current simulation time.
        satellite_state (numpy.ndarray): Satellite state vector.
        thrust (float): Thruster thrust value.
        torque_thruster (numpy.ndarray): Thruster torque.
        satellite (object): Satellite object.
        grid (object): Grid object.
        selected_grid (object): Selected grid object.
        controlSystem (object): Control system object.
        environment (object): Environment object.
        sensor (object): Sensor object.
        iteration (int): Current simulation iteration.

    Returns:
        int: Always returns 1, so that solve_ivp will not terminate the integration.
    """
    #get requried variables:
    omega_state = satellite_state[6:]   #quate, ang vel
    omega = omega_state[4:7]            #ang vel


    #if no previous data, measure angle:

    #Read the previous bias instability value if any
    if not hasattr(sensorMeasurement, "biais_prev"):
        sensorMeasurement.biais_prev = sensor.biais_n           #measure of how the error should grow

    #Read the previous measurement time if any
    if not hasattr(sensorMeasurement, 't_prev'):
        sensorMeasurement.t_prev = 0
    
    #Set the current time
    if not hasattr(sensorMeasurement, 't_n'):
        sensorMeasurement.t_n = t


    #after each successful step:
    if t > sensorMeasurement.t_n :

        #update previous time:
        sensorMeasurement.t_prev = sensorMeasurement.t_n
        #update new time:
        sensorMeasurement.t_n = t
        



        #draw new measurement of sensor:
        sensorMeasurement.biais_n, sensorMeasurement.omega_n = sensor.MEMSgyro_AngularVelocity(omega, sensorMeasurement.t_n, sensorMeasurement.t_prev, sensorMeasurement.biais_prev)
        
        #update biais in MEMS and update class
        sensorMeasurement.biais_prev = sensorMeasurement.biais_n
        sensor.MEMSomega = sensorMeasurement.omega_n

        #Write sensor readings to a text file
        file3 = open('output/run{}/MEMS'.format(iteration) + '.txt', 'a')
        file3.write('{:.3f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(sensorMeasurement.t_n, omega[0], omega[1], omega[2], sensor.MEMSomega[0], sensor.MEMSomega[1], sensor.MEMSomega[2]))
        file3.close()

        


    return 1        # -- we don't need to return anything, we just want to update the class at every tiem step




def changeGridState(satellite_state, grid, CS, environment, sensor):
    """
    Updates the firing state of a grid based on the expected signal-to-noise ratio (SN).

    Args:
        satellite_state (numpy.ndarray): Satellite state vector.
        grid (object): Grid object .
        CS (object): Control system object .
        environment (object): Environment object.
        sensor (object): Sensor object.
    """

    #Determines firing decision based on expected SN 
    p, f, g, h, k, L = satellite_state[:6]              #equinoctial element
    vec_v = equi_to_v_cart(satellite_state[:6])         #Velocity vector
    vec_r = equi_to_r_cart(p, f, g, h, k, L)            #position vector
    v = np.linalg.norm(vec_v) #km/s

    #angle relative to velocity direction
    pointing_angle = CS.velocity_tracking(satellite_state, grid)
    
    #compute ion density
    rho_ions = environment.ionDensity(vec_r)

    # compute cone angle from SN 
    #intersection area inside the faraday cup
    A = circlesIntersection(sensor.FCradius, sensor.FCradius, sensor.FCheight, pointing_angle)
    
    #compute signal to noise
    SN = signalToNoise(rho_ions, A, v*1000, sensor.FCepsilon)       #velocity must be in SI units
    
    #determine if thrusting should happen
    if SN >= sensor.FC_SN:
        grid.state = 1*CS.switch        #Do not turn on if CS.switch is zero
        
    elif SN < sensor.FC_SN:
        grid.state = 0


def computeGridState(satellite_state, grid, CS, environment, sensor):
    """
    Computes the firing state of a grid based on the expected signal-to-noise ratio (SN).

    Args:
        satellite_state (numpy.ndarray): Satellite state vector.
        grid (object): Grid object .
        CS (object): Control system object .
        environment (object): Environment object.
        sensor (object): Sensor object.
    """
    #prep data
    p, f, g, h, k, L = satellite_state[:6]
    vec_v = equi_to_v_cart(satellite_state[:6])
    vec_r = equi_to_r_cart(p, f, g, h, k, L)
    v = np.linalg.norm(vec_v) #km/s
    pointing_angle = CS.velocity_tracking(satellite_state, grid)
    
    
    #compute ion density
    rho_ions = environment.ionDensity(vec_r)

    # compute cone angle from SN 
    A = circlesIntersection(sensor.FCradius, sensor.FCradius, sensor.FCheight, pointing_angle)
    
    #compute signal to noise
    SN = signalToNoise(rho_ions, A, v*1000, sensor.FCepsilon)
    
    #determine if thrusting should happen
    if SN >= sensor.FC_SN:
        # grid.state = 1*CS.switch
        gridState = 1
    elif SN < sensor.FC_SN:
        # grid.state = 0
        gridState = 0
    return gridState


def assessGridState2(t, satellite_state, thrust,  torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, manager, iteration):
    """
    Event function to assess the firing state of the last grid.

    This function checks the expected signal-to-noise ratio (SN) for the last thruster grid
    using the `computeGridState` function. It then returns a value based on the grid state
    (1 for firing, 0 for not firing),  scaled by 2 and subtracting 1 to create an event function that crosses zero.

    Returns:
        float:
            - A value of 1  if the SN of the last grid is sufficient for firing.
            - A value of -1 otherwise.
    """
    gridState = computeGridState(satellite_state, grids[-1], controlSystem, environment, sensor)

    return 2* gridState - 1



def assessGridState(t, satellite_state, thrust,  torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, manager, iteration):
    """
    Event function to assess the firing state of the first grid.

    This function checks the expected signal-to-noise ratio (SN) for the last thruster grid
    using the `computeGridState` function. It then returns a value based on the grid state
    (1 for firing, 0 for not firing),  scaled by 2 and subtracting 1 to create an event function that crosses zero.

    Returns:
        float:
            - A value of 1  if the SN of the last grid is sufficient for firing.
            - A value of -1 otherwise.
    """
    gridState = computeGridState(satellite_state, grids[0], controlSystem, environment, sensor)
    return 2* gridState - 1



def torqueProfile(t, satellite_state, thrust, torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, manager, iteration):
    """
    Calculates and logs the combined torque acting on the spacecraft.

    Returns:
        int: Always returns 0 to indicate successful execution.
    """
    force, torque = satellite.envPerturbations(t)
    torquer = np.append([t], torque)
    f = open('run{}/torques.txt'.format(iteration), 'a')
    f.write(''.join('%s\t' % item for item in torquer))
    f.write('\n')
    f.close()
    return 0



