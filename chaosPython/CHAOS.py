#python code to drive the C.H.A.O.S. code

import os
import numpy as np
import scipy.integrate as integrate
import spiceypy 

from .attitude_ODE import analytical
from .orbital_driver import modified_equinoctial
from .class_grid import Grid
from .class_panel import Panel
from .class_sensor import Sensor
from .support_fun import gridsBurntime, selectGrid
from .class_dataset import Dataset
from .class_satellite import Satellite
from .class_control_system import control_system
from .class_shapeFunction import *
from .class_perturbations import Environment
# from plotter import centerQuadrantPlot
from .smart_fun import  pixel_fuel, r_stop, sensorMeasurement, changeGridState, assessGridState, assessGridState2


'''
Python script to drive the C.H.A.O.S. code. It will handle calling the correct variables and methods,and ensure
the correct torques and force models are used. 
This code will also provide a recap interface and provide ease-of-use for easy simulation setting.
'''


def coupledMotion(t, satellite_state, thrust, torque_thruster, satellite, grids, selected_grid, CS, environment, sensor, it):
    #extract for clarity
    equinoctial = satellite_state[:6] # orbital data
    omega_state = satellite_state[6:] # attitude data
    #update the class first!!
    satellite.modEquinoctialElements = equinoctial
    satellite.quaternion = omega_state[:4]
    satellite.angularVel = omega_state[4:7]
    
    # #check if thruster points in correct direction. If yes, turn grid on. If not, Cube de ALPS is turned off
    for grid in range(len(grids)):
    #     #changes the state of each grid
    #     effectiveGridState(satellite_state, grids[grid], CS, environment, sensor) #uses MEASURED attitude!
        omega_state[11+grid] = grids[grid].state
        # print('grid state', grids[grid].name, grids[grid].state, CS.velocity_tracking(satellite_state, grids[grid])*180/np.pi)
        # print('selected grid', selected_grid.name)
    
    #the selected grid is whichever grid is turned on
    #selected_grid = selectGrid(grids)
    # omega_state[11] = grids[0].state #assign the state as integration variable for tracking

    #evaluate the force and torques acting on the spacecraft
    orbital_perturb, torque_env = satellite.envPerturbations(t)

    # compute osculating elements
    #--The relevant grid is passed
    #--fuel updated in selected_grid
    osculating = modified_equinoctial(t, satellite_state, thrust, orbital_perturb, satellite, selected_grid, environment)

    # comptue derivative of rotational state
    rotational_state = analytical(t, omega_state, torque_thruster, torque_env, satellite,  selected_grid)
    # combine arrays
    state_derivative = np.append(osculating, rotational_state)


    # print('thruster angle', CS.velocity_tracking(omega_state, satellite, grid)*180/np.pi)
    # dist_boundary = thrusting_cone_boundary(t, satellite_state, t_stop, thrust, nextPixel, torque_thruster, satellite, grid, CS, environment, it)
    # # print('ODE function call')
    # # print('event fun value', dist_boundary)
    # # print('breaker in func', dataset.breaker)
    # f = open('run{}/test.txt'.format(it), 'a')
    # f.write('{}\t{}\t{}\t{}\n'.format(t, CS.velocity_tracking(satellite_state, grid), dist_boundary, grid.state*CS.switch))
    # f.close()
    # exit()
    return state_derivative



def CHAOS(satellite, environment, grids, controlSystem, sensor, dataset, t_stop, writepath, iteration, integrator='DOP853', tols=1e-12):

    #extract initial conditions:

    if dataset.y0 is not None: #initial conditions from checkpoint file
        satellite_state = dataset.y0

        #update class
        satellite.modEquinoctialElements = satellite_state[0:6]
        satellite.quaternion = satellite_state[6:10]
        satellite.angularVel = satellite_state[10:13]
        satellite.hysteresisInduction = satellite_state[13:16]


    else: #initial conditions from set-up
        #orbital:
        equinoctial = satellite.modEquinoctialElements

        #attitude:
        #randomise initial quaternion and induction:
        # satellite.initial_induction() #deprecated? 
        quat = satellite.quaternion
        ang_vel = satellite.angularVel
        hyst_induction = satellite.hysteresisInduction
        
        #create initial attitude_state:
        satellite_state = np.zeros((17,))


        satellite_state[0:6] = equinoctial # orbital equinoctial element
        satellite_state[6:10] = quat #attitude quaternion
        satellite_state[10:13] = ang_vel # attitude ang vel
        satellite_state[13:16] = hyst_induction # magnetic hysterisis induction
        satellite_state[16] = 0 #fuel for specified pixel
        # satellite_state[17] = 0 # selected grid state
        for n in range(len(grids)):
            satellite_state = np.append(satellite_state, 0)

    
    satellite.it = iteration
 


    #drive the Cubesat Hybrid Attitude and Orbital Simulator
    #make folders:
    run_folder = 'run{}'.format(iteration)
    writepath = writepath + '/' + run_folder
    if os.path.exists(writepath) != True:
            os.mkdir(writepath)

    if os.path.exists(writepath + '/timetrack') != True:
            os.mkdir(writepath + '/timetrack')

    if os.path.exists(writepath + '/orbital') != True:
        os.mkdir(writepath + '/orbital')


    if os.path.exists(writepath + '/attitude') != True:
        os.mkdir(writepath + '/attitude')

    if os.path.exists(writepath + '/objects') != True:
        os.mkdir(writepath + '/objects')


    # open output files
    writeTimetrack = open(writepath+'/timetrack/timetrack{}.txt'.format(iteration), 'a')
    writeGrid1 = open(writepath+'/objects/grid1_data{}.txt'.format(iteration), 'a')
    writeGrid2 = open(writepath+'/objects/grid2_data{}.txt'.format(iteration), 'a')
    writeQuaterion = open(writepath+'/attitude/quaternion{}.txt'.format(iteration), 'a')
    writeEquinoctial = open(writepath+'/orbital/equinoctial{}.txt'.format(iteration), 'a')

    #calls constants into shorter notation variables
    I = satellite.inertiaMatrix 

    #initialise start time and breaker variable
    t_start = dataset.t_start
    top_breaker = 0#to break infinite loop if they occur
    

    #store information about a pixel firing time
    # dataset.firingTime = grid.pixelDict

    altitude = -1e17
    p = None
    
    #Compute total number of pixels / thruster Heads
    totalNumHeads = 0
    for i in grids:
        totalNumHeads += i.numHeads
    #while there is more than 1 second of thrusting left in each pixel and if we haven't reached the max time t_stop
    #and if we are above 
    while gridsBurntime(grids) > totalNumHeads and t_start < t_stop and altitude < -1: 


        #check if thruster points in correct direction. If yes, turn grid on. If not, grid is turned off
        for grid in range(len(grids)):
            #changes the state of each grid
            changeGridState(satellite_state, grids[grid], controlSystem, environment, sensor) #uses MEASURED attitude!
            # print('CHAOS grid state', grids[grid].name, grids[grid].state)

        #determine which grid the pixel will be selected from
        #-- quadrant selection is grid-agnostic, so not too important, 
        #-- But we need to pass 1 grid to the CS        
        selected_grid = selectGrid(grids)


        #predict next pixel to use. Includes re-ignition condition in the choice /returns None if no pixel chosen
        nextPixel = controlSystem.decidePixel(satellite_state[6:13], selected_grid, satellite, sensor) 

        # if no more pixels are available, check other grids
        if nextPixel == None:
            break            
        
        #compute which quadrant to use if no sensors were used
        omega_state = satellite_state[6:] #quate, ang vel
        omega = omega_state[4:7] #ang vel
        controlSystem.quadrantControl_truth(omega, selected_grid, satellite)


        
         
        #set burn time of the chosen pixel--- this is the integration time
        #burn for Quadrant Ignition Time, or for remaining fuel time, whichever is lower:
       
        #Update the initial condition of the fuel: Use fuel of correct pixel and correct grid:
        satellite_state[16] = selected_grid.burnTime[nextPixel]


        #define at which time the pixel should stop firing: start time + burn time
        t_end = selected_grid.firingTime[nextPixel] + t_start  #default QIT firing time

 
        #ensure t_end does not go over our max time t_stop
        if t_stop < t_end:
            t_end = t_stop

        #if the thrusters are off, i.e.. no pixels fired, integration time is t_stop. 
        #-- this is faster
        if controlSystem.switch == 0:
            t_end = t_stop
            print('thruster off--intergation will not be split')
            
            # Thrust and torques from thrusters are null.
            thrust = np.array([0, 0, 0])
            torque_thruster = np.array([0, 0, 0])
        
        
        #if the grid is on
        elif controlSystem.switch == 1:
            #fetch thrust and torque from the pixel
            thrust_val, torque_thruster_val = selected_grid.generateTorque(nextPixel) 
            thrust = thrust_val * selected_grid.state
            torque_thruster = torque_thruster_val * selected_grid.state



            

        #set initial condition for integration, i.e. the attitude state
        init_cond = satellite_state

        #set number of evaluation points in integration interval
        # eval_time = np.linspace(t_start, t_end, 8640)  #temporary measure.


        #integrate the ODE / propagate the motion with analytical form of eq.
        #args passed as complementary variables to the ODE function
                #tolerances and integrator defined as inputs
        solver = integrate.solve_ivp(coupledMotion, t_span=(t_start, t_end), y0=init_cond, args=(thrust, torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, iteration), events=controlSystem.cut_off_function, dense_output=False, method=integrator, t_eval=None, rtol=tols, atol=tols) 
        

        #if grid is on, check whether the pixel was at all ignited 
        #ensure the final time point of the integration is different from the starting time stamp
        if t_start != solver.t[-1] and controlSystem.switch == 1:

            #update the burn time and the ignition log for data storage
            # grid.burnTime[nextPixel] -= (solver.t[-1] - t_start)  #less fuel after firing for FiringTime 
            selected_grid.ignitionLog[nextPixel] += len(solver.t_events[1])   #how many times a pixel has ignited
            
            #store the firing time of the ignited pixel.
            # dataset.firingTime[nextPixel] = np.append(dataset.firingTime[nextPixel], [solver.t[-1] - t_start])
            #reset breaker --- used to stop simulation if no pixel can be ignited
            dataset.breaker = 0
            top_breaker = 0
            
            
            #store the sequence of pixel firing:
            # dataset.sequence = np.append(dataset.sequence, nextPixel)

        #if the integration did not actually use a pixel
        elif t_start == solver.t[-1]:
            print('\t\t\t\t no advance in time has been made.')
            dataset.breaker +=1 # add to the breaker counter
            top_breaker += 1


        #update state for next integration leg
        satellite_state =  solver.y[:, -1] # velocity in rotated frame,
        #normalise quaternion if needed
        equinoctial = satellite_state[:6]
        quaternion = satellite_state[6:10]
        rest = satellite_state[10:]
        quat_norm = 1#np.linalg.norm(quaternion)
        quaternion = quaternion / quat_norm

        #update fuel --take from the state to update fuel grid
        selected_grid.burnTime[nextPixel] = satellite_state[16] #pixel fuel update

        satellite_state = np.append(equinoctial, quaternion)
        satellite_state = np.append(satellite_state, rest)
        

        # store value of initial start point:
        t_org = t_start
        #update time for next integration leg
        t_start = solver.t[-1]
        
        # progress bar 
        print('torg', t_org, 'tf_integration', t_start)
        print('grid burn progress', t_start / gridsBurntime(grids))
        print('Time progress', t_start / t_stop)
        print('msg', solver.message, 'y:', solver.y[:, -1], 't', solver.t[-1])
        print('number of function eval:', solver.nfev)
        #update time over which the data is outputed
        # eval_time = np.linspace(t_org, solver.t[-1], 8640)
        # # eval_time = np.linspace(t_org, solver.t[-1], len(solver.t))
        # print('eval time len', len(eval_time))
        print('tstop', t_stop, 'final time', solver.t[-1],  '\n\t\t\t\n//////////////////')


        #store data 
        #calls the interpolation to get data at desired time stamps
        # z = solver.sol(eval_time) #dense_output

        #get number of points the integrator used for solving the ODE
        p =  len(solver.t)

        omega_norm =  (np.square(solver.y[10]) + np.square(solver.y[11]) + np.square(solver.y[12]))**(0.5)
        #tracks the quaternion and omega components with respect to time. stored in the dataset class.
        for i in range(0, len(solver.t)):# int(len(eval_time)/25)+1
  

            #if using solver determined output---debug statement:   

            dataset.write_value(writeTimetrack, solver.t[i], solver.y[10, i], solver.y[11, i], solver.y[12, i], omega_norm[i])
            dataset.write_value(writeGrid1, solver.t[i], grids[0].name, solver.y[17, i], sum(grids[0].burnTime.values()), gridsBurntime(grids))
            dataset.write_value(writeGrid2, solver.t[i], grids[-1].name, solver.y[-1, i], sum(grids[-1].burnTime.values()), gridsBurntime(grids))
            # dataset.write_value(writepath+'/gauss/RateData{}'.format(iteration), [dataset.maxIgnition], [dataset.totalFuel], [dataset.maxVelocity], [dataset.absoluteMinFiringTime], [dataset.maxVelocity])
            dataset.write_value(writeQuaterion, solver.t[i], solver.y[6, i], solver.y[7, i], solver.y[8, i], solver.y[9, i]) 
            dataset.write_value(writeEquinoctial, solver.t[i], solver.y[0, i], solver.y[1, i], solver.y[2, i], solver.y[3, i], solver.y[4, i], solver.y[5, i])
            

            
        #if the driver has looped 3 times without progressing in time, break
        if dataset.breaker == 3:
            
            print('breaker condition reached: {} intergation without t moving.'.format(dataset.breaker))
            break


        #store quadrant selection precision data
        quadName = selected_grid.idealQuadrant
        file2 = open(writepath + '/quadrantSelection'.format(iteration) + '.txt', 'a')
        file2.write('{:.3f}\t{}\t{}\n'.format(solver.t[-1], selected_grid.predictedQuadrant, quadName))
        file2.close()

        #recompute altitude
        altitude = r_stop(0, satellite_state, t_stop, thrust, nextPixel, torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, iteration)


        #prep for checkpointing
        fuelGrid = []
        for i in grids:
            fuelGrid.append(i.burnTime)
        #checkpoint the last value of time, integral
        np.savez(writepath + '/checkpoint'.format(iteration), t=t_start, y0=satellite_state, fuel=fuelGrid, sensor_biais=sensor.biais_n)
        print('simulated days', solver.t[-1] / 86400)
        
    #END OF WHILE LOOP 
    # one of the following events occurred:
    #- t_stop has been reached
    #- no fuel left on pixels
    #- 1000 pixel in a row have been selected but couldn't fire

    #close output files
    writeTimetrack.close()
    writeEquinoctial.close()
    writeGrid1.close()
    writeGrid2.close()
    writeQuaterion.close()

    return p    
    ################## END ###############

    
    





