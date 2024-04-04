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
from .sensors import effectiveCupAngle, signalToNoise
from .support_fun import selectGrid, circlesIntersection







def pixel_fuel(t, satellite_state, *args):# t_stop, thrust, nextPixel, torque_thruster, satellite, grids, selected_grids,controlSystem, environment, sensor, iteration):
    #function to return 0 when the pixel has no more fuel. 
    # This even function is terminal

    return satellite_state[16]






def r_stop(t, satellite_state, *args):#t_stop, thrust, nextPixel, torque_thruster, satellite, grid, selected_grid, controlSystem, environment,sensor, iteration):
    #function to stop simulation when perigee reaches 150 km.
    #extract elements
    p, f, g, h, k, L = satellite_state[:6]

    #get keplerian elements:
    a, e, inc, RAAN, wp, theta = eq_to_kep(p, f, g, h, k, L)

    #get Cartesian position vector
    r_ECI = equi_to_r_cart(p, f, g, h, k, L)

    #cpmpute perigee (height):
    r_p = a*(1 - e) - 6378
    # h = np.linalg.norm(r_ECI) - 6378
    # print('h',h)
    #check if stop condition is met:
    event = 150 - r_p
    # print('event', event)
    return event


def sensorMeasurement(t, satellite_state, t_stop, thrust, nextPixel, torque_thruster, satellite, grid, selected_grid, controlSystem, environment, sensor, iteration):


    #get requried variables:
    #true angle relative to velocity
    # true_alpha = controlSystem.velocity_tracking(satellite_state, grid)
    omega_state = satellite_state[6:] #quate, ang vel
    omega = omega_state[4:7] #ang vel
    #get magnitude of spacecraft velocity and position
    # p, f, g, h, k, L = satellite_state[:6]
    # vec_v = equi_to_v_cart(satellite_state[:6])
    # vec_r = equi_to_r_cart(p, f, g, h, k, L)
    # v = np.linalg.norm(vec_v)
    # r = np.linalg.norm(vec_r)

    #get ion density
    # rho_ions = environment.ionDensity(vec_r)
    # print(rho_ions, v)
    #if no previous data, measure angle:
    if not hasattr(sensorMeasurement, "biais_prev"):
        # sensorMeasurement.alpha_prev = sensor.coarseAttitude(true_alpha, v, rho_ions)
        sensorMeasurement.biais_prev = sensor.biais_n

    # if not hasattr(sensorMeasurement, 'biais_n'):
    #     # sensorMeasurement.alpha_n = sensorMeasurement.alpha_prev
    #     sensorMeasurement.biais_prev = sensorMeasurement.biais_n
    #if no previous data, assign time
    if not hasattr(sensorMeasurement, 't_prev'):
        sensorMeasurement.t_prev = 0
    
    if not hasattr(sensorMeasurement, 't_n'):
        sensorMeasurement.t_n = t


    #in root finder, for any time t between t_prev and t_n:
    # if t <= sensorMeasurement.t_n:
    # #obtain a smooth function between alpha_prev and alpha_n
    #     f_alpha = interp1d(np.array([sensorMeasurement.t_prev, sensorMeasurement.t_n]), np.array([sensorMeasurement.alpha_prev, sensorMeasurement.alpha_n]), kind='linear')
    #     alpha_t = f_alpha(t)



    #after each successful step:
    if t > sensorMeasurement.t_n :#and sensorMeasurement.t_n  -sensorMeasurement.t_prev > 0
        # print('t', t, 't_n', sensorMeasurement.t_n, 't_prev', sensorMeasurement.t_prev)

        #update previous time:
        sensorMeasurement.t_prev = sensorMeasurement.t_n
        #update new time:
        sensorMeasurement.t_n = t
        

        #update previous measurements
        # sensorMeasurement.alpha_prev = sensorMeasurement.alpha_n

        #draw new measurement of sensor:
        # sensorMeasurement.alpha_n = sensor.coarseAttitude(true_alpha, v, rho_ions)
        sensorMeasurement.biais_n, sensorMeasurement.omega_n = sensor.MEMSgyro_AngularVelocity(omega, sensorMeasurement.t_n, sensorMeasurement.t_prev, sensorMeasurement.biais_prev)
        # print('draw', sensorMeasurement.omega_n*180/np.pi)
        
        #update biais in MEMS and update class
        sensorMeasurement.biais_prev = sensorMeasurement.biais_n
        sensor.MEMSomega = sensorMeasurement.omega_n

        #get output:
        # alpha_t = sensorMeasurement.alpha_n


        # #internal test
        # dist_boundary = alpha_t - controlSystem.thrustingConeAngle
        # #controlSystem.quadrantControl_truth(omega, grid, satellite)
        # if dist_boundary <= 0:
        #     state = 1
        # elif dist_boundary > 0:
        #     state = 0

        # true_state = grid_state_truth(grid, controlSystem, satellite_state)
        # file = open('run{}/faraday'.format(iteration) + '.txt', 'a')
        # file.write('{}\t{:.7f}\t{:.7f}\n'.format(sensorMeasurement.t_n, true_alpha, alpha_t))
        # file.close()
        
        # file2 = open('run{}/thruster_state'.format(iteration) + '.txt', 'a')
        # file2.write('{:.3f}\t{}\t{}\n'.format(sensorMeasurement.t_n, true_state, state))
        # file2.close()

        file3 = open('output/run{}/MEMS'.format(iteration) + '.txt', 'a')
        file3.write('{:.3f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(sensorMeasurement.t_n, omega[0], omega[1], omega[2], sensor.MEMSomega[0], sensor.MEMSomega[1], sensor.MEMSomega[2]))
        file3.close()

        
        # print('t_step', sensorMeasurement.t_n - sensorMeasurement.t_prev)
    #store this angle into the sensor class for use elsewhere
    # sensor.measuredAngle = alpha_t 


    return 1#alpha_t -- we don't need tor eturn anythin, we just want to update the class at every tiem step




def changeGridState(satellite_state, grid, CS, environment, sensor):

    ###
    #Determines firing decision based on expected SN 
    #prep data
    p, f, g, h, k, L = satellite_state[:6]
    vec_v = equi_to_v_cart(satellite_state[:6])
    vec_r = equi_to_r_cart(p, f, g, h, k, L)
    v = np.linalg.norm(vec_v) #km/s
    pointing_angle = CS.velocity_tracking(satellite_state, grid)
    #compute ion density
    rho_ions = environment.ionDensity(vec_r)

    # compute cone angle from SN 
    # thrusting_cone_angle = effectiveCupAngle(sensor.FC_SN, sensor.FCepsilon, rho_ions, sensor.FCaperture, v*1e3)
    A = circlesIntersection(sensor.FCradius, sensor.FCradius, sensor.FCheight, pointing_angle)
    #compute signal to noise
    SN = signalToNoise(rho_ions, A, v*1000, sensor.FCepsilon)
    #determine if thrusting should happen
    # print('nominal SN', sensor.FC_SN, 'current SN', SN, 'A', A)
    # print('changer', grid.name, grid.state,'poitning', CS.velocity_tracking(satellite_state, grid)*180/np.pi)

    if SN >= sensor.FC_SN:
        grid.state = 1*CS.switch
        
    elif SN < sensor.FC_SN:
        grid.state = 0


def computeGridState(satellite_state, grid, CS, environment, sensor):

    ###
    #Determines firing decision based on expected SN 
    #prep data
    p, f, g, h, k, L = satellite_state[:6]
    vec_v = equi_to_v_cart(satellite_state[:6])
    vec_r = equi_to_r_cart(p, f, g, h, k, L)
    v = np.linalg.norm(vec_v) #km/s
    pointing_angle = CS.velocity_tracking(satellite_state, grid)
    #compute ion density
    rho_ions = environment.ionDensity(vec_r)

    # compute cone angle from SN 
    # thrusting_cone_angle = effectiveCupAngle(sensor.FC_SN, sensor.FCepsilon, rho_ions, sensor.FCaperture, v*1e3)
    A = circlesIntersection(sensor.FCradius, sensor.FCradius, sensor.FCheight, pointing_angle)
    #compute signal to noise
    SN = signalToNoise(rho_ions, A, v*1000, sensor.FCepsilon)
    #determine if thrusting should happen
    # print('nominal SN', sensor.FC_SN, 'current SN', SN, 'A', A)
    # print('EVENT', grid.name, grid.state,'poitning', CS.velocity_tracking(satellite_state, grid)*180/np.pi)

    if SN >= sensor.FC_SN:
        # grid.state = 1*CS.switch
        gridState = 1
    elif SN < sensor.FC_SN:
        # grid.state = 0
        gridState = 0
    return gridState


def assessGridState2(t, satellite_state, t_stop, thrust, nextPixel, torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, iteration):
    
    gridState = computeGridState(satellite_state, grids[-1], controlSystem, environment, sensor)

    return 2* gridState - 1



def assessGridState(t, satellite_state, t_stop, thrust, nextPixel, torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, iteration):
    
    gridState = computeGridState(satellite_state, grids[0], controlSystem, environment, sensor)
    return 2* gridState - 1



def torqueProfile(t, satellite_state, t_stop, thrust, nextPixel, torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, iteration):

    force, torque = satellite.envPerturbations(t)
    torquer = np.append([t], torque)
    f = open('run{}/torques.txt'.format(iteration), 'a')
    f.write(''.join('%s\t' % item for item in torquer))
    f.write('\n')
    f.close()
    return 0



