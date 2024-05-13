#python code to drive the CHAOS code

'''
"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


CHAOS.py



This Python code defines the `Cubesat High-fidelity Attitude and Orbit Simulator.
It manages the interaction between various modules:

- `Satellite` class: Represents the spacecraft properties (orbital elements, attitude, etc.)
- `Environment` class: Models the environmental effects on the spacecraft.
- `Grid` class (multiple instances): Represents thruster grids used for attitude & orbit effect.
- `ControlSystem` class: Implements the control logic for thruster firing decisions.
- `Sensor` class: Models the sensor readings used for feedback in the control system.
- `dataset` object: Stores and manages simulation data

The `CHAOS` function performs the following key tasks:

1. **Initialization:**
    - Extracts initial conditions (orbital and attitude states) from either a checkpoint file or user-defined values.
    - Creates the complete state vector combining orbital and attitude data.
    - Sets up folders for data storage.
    - Initializes variables for tracking simulation progress and stopping criteria.

2. **Main Simulation Loop:**
    - Iterates until:
        - There's no fuel left for thrusters across all grids.
        - The maximum simulation time (`t_stop`) is reached.
        - The spacecraft falls below 150 km.


    - Within each iteration:
        - Selects the appropriate grid and pixel for thruster firing based on the control system and sensor readings.
        - Updates the state vector with the selected pixel's fuel consumption.
        - Defines the integration time based on the pixel's firing duration.
        - Integrates the coupled equations of motion (orbital and attitude dynamics) using a chosen numerical integrator.
        - Stores simulation data (time, orbital elements, attitude, grid states) in text files.
        - Updates the state vector for the next integration step.
        - Checks for specific breaking conditions
        - Checkpoints the simulation state for potential restart.

3. **Finalization:**
    - Closes data output files.

'''





import os
import numpy as np
import scipy.integrate as integrate
from .attitude_ODE import analytical
from .smart_fun import r_stop, changeGridState
from .orbital_driver import modified_equinoctial
from .support_fun import gridsBurntime, selectGrid






def coupledMotion(t, satellite_state, thrust, torque_thruster, satellite, grids, selected_grid, CS, environment, sensor, manager, it):
    """
    Calculates the combined derivatives of the spacecraft's orbital and attitude states.

    Args:
        t (float): Current simulation time (seconds).
        satellite_state (numpy.ndarray): State vector containing orbital and attitude data.
            - satellite_state[:6] (orbital): Osculating orbital elements.
            - satellite_state[6:] (attitude): Quaternion and angular velocity and pixel data.
        thrust (numpy.ndarray): Thrust vector acting on the spacecraft (body frame).
        torque_thruster (numpy.ndarray): Torque generated by thrusters (body frame).
        satellite (class object): Instance of the `Satellite` class representing the spacecraft.
        grids (list): List of `Grid` class objects representing the thruster grids.
        selected_grid (class object): The currently selected thruster grid for firing.
        CS (class object): Instance of the `ControlSystem` class implementing the control logic. 
        environment (class object): Instance of the `Environment` class representing the environment.
        sensor (class object): Instance of the `Sensor` class modeling the sensor readings. 
        manager (class object): perturbationManager instance to call perturbation models
        it (int): Current simulation iteration number. 

    Returns:
        numpy.ndarray: Combined state derivatives (orbital and attitude).
    """
    
    
    
    
    
    #extract for clarity
    equinoctial = satellite_state[:6]                   # orbital data
    omega_state = satellite_state[6:]                   # attitude data
    
    #update the class first!!
    satellite.modEquinoctialElements = equinoctial
    satellite.quaternion = omega_state[:4]
    satellite.angularVel = omega_state[4:7]
    



    # #check if thruster points in correct direction. 
    #If yes, turn grid on. If not, Cube de ALPS is turned off
    for grid in range(len(grids)):

        #changes the state of each grid
        omega_state[8+grid] = grids[grid].state




    #evaluate the forces and torques acting on the spacecraft
    orbital_perturb, torque_env = manager.callForceModels(t, satellite)



    # compute osculating elements
    osculating = modified_equinoctial(t, satellite_state, thrust, orbital_perturb, satellite, selected_grid, manager)



    #comptue derivative of rotational state
    rotational_state = analytical(t, omega_state, torque_thruster, torque_env, satellite,  selected_grid)
    


    # combine arrays
    state_derivative = np.append(osculating, rotational_state)


    return state_derivative











def CHAOS(satellite, manager, environment, grids, controlSystem, sensor, dataset, t_stop, writepath, iteration, integrator='DOP853', tols=1e-6):

    """
    Drives the CHAOS spacecraft simulation.

    Args:
        satellite (class object): Instance of the `Satellite` class representing the spacecraft.
        manager (class object): Instace of the `perturbationManager` class handling the force models
        environment (class object): Instance of the `Environment` class representing the environment.
        grids (list): List of `Grid` class objects representing the thruster grids.
        controlSystem (class object): Instance of the `ControlSystem` class implementing the control logic.
        sensor (class object): Instance of the `Sensor` class modeling the sensor readings.
        dataset (object): Object for managing simulation data
        t_stop (float): Maximum simulation time (seconds).
        writepath (str): Path to the directory for storing simulation data.
        iteration (int): Current simulation name.
        integrator (str, optional): Name of the numerical integration method. Defaults to 'DOP853'.
        tols (float, optional): Tolerance levels for the numerical integration. Defaults to 1e-6.

    Returns:
        None
    """



    #########################
    #INITIALISATION
    ########################


    #extract initial conditions:
    #Checks whether a checkpointing file was loaded
    if dataset.y0 is not None:                      
        satellite_state = dataset.y0                    #initial conditions from checkpoint file

        #update class
        satellite.modEquinoctialElements = satellite_state[0:6]
        satellite.quaternion = satellite_state[6:10]
        satellite.angularVel = satellite_state[10:13]
        # satellite.hysteresisInduction = satellite_state[13:16]


    #If no checkpoint files, take the initial conditions as specified 
    # in set_up.py file
    else: 
        #orbital:
        equinoctial = satellite.modEquinoctialElements

        #attitude:
        quat = satellite.quaternion
        ang_vel = satellite.angularVel
        


        #create complete state vector:
        satellite_state = np.zeros((14,))

        #Assign values
        satellite_state[0:6] = equinoctial                      # orbital equinoctial element
        satellite_state[6:10] = quat                            #attitude quaternion
        satellite_state[10:13] = ang_vel                        # attitude ang vel
        satellite_state[13] = 0                                 #fuel for specified pixel -- value will be updated later
        
        # Add the state of the grid for tracking:
        for n in range(len(grids)):
            satellite_state = np.append(satellite_state, 0)     #value will be updated later

    
    #Store the simulation name
    # satellite.it = iteration
 

    #make folders for data storage:
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



    #initialise start time and breaker variable
    t_start = dataset.t_start
    top_breaker = 0                                     #to break infinite loop if they occur
    

    #stopping condition:
    altitude = -1e17                                    #Value updated in simulation loop
                                            
    


    #Compute total number of pixels across all grids
    totalNumHeads = 0
    for i in grids:
        totalNumHeads += i.numHeads





    #########################
    #START CHAOS SIMULATION 
    ########################

    #while there is more than 1 second of thrusting left in each pixel
    # and if we haven't reached the max time t_stop
    #and if we are above the minimum altitude (150 km)
    while gridsBurntime(grids) > totalNumHeads and t_start < t_stop and altitude < -1: 

        ############################
        #GRID AND QUADRANT SELECTION
        ############################

        #check if thruster points in correct direction. 
        #If yes, turn grid on. If not, grid is turned off
        for grid in range(len(grids)):
            #changes the state of each grid
            #uses MEASURED attitude!
            changeGridState(satellite_state, grids[grid], controlSystem, environment, sensor) 




        #determine which grid the pixel will be selected from
        #-- This grid will be passed to the ControlSystem and the ODEs      
        selected_grid = selectGrid(grids)


        #Selects wich pixel to use. Includes re-ignition condition in the choice /returns None if no pixel chosen
        nextPixel = controlSystem.decidePixel(satellite_state[6:13], selected_grid, satellite, sensor) 

        # if no more pixels are available, check other grids
        if nextPixel == None:
            break            
        

        ##FOR COMPARISON'S SAKE
        #compute which quadrant would be used if no sensors were used
        omega_state = satellite_state[6:]                               #quate, ang vel
        omega = omega_state[4:7]                                        #ang vel
        #Seelcts quadrant based on real angular velocity, NOT the masured AngVel
        controlSystem.quadrantControl_truth(omega, selected_grid, satellite)


        



        ############################
        #PREP FOR INTEGRATION
        ############################
        #Update the initial condition of the fuel in the state vector: Use fuel of correct pixel and correct grid:
        satellite_state[13] = selected_grid.burnTime[nextPixel]

         
        #set burn time of the chosen pixel--- this is the integration time
        #define at which time the pixel should stop firing: start time + burn time
        t_end = selected_grid.firingTime[nextPixel] + t_start                   #default QIT firing time

        #ensure t_end does not go over our max time t_stop
        if t_stop < t_end:
            t_end = t_stop

        #if the thrusters are off, i.e.. no pixels fired, integration time is t_stop. 
        #-- this is faster
        if controlSystem.switch == 0:
            t_end = t_stop
            print('thrusters are off--intergation will not be split')
            
            # Thrust and torques from thrusters are null.
            thrust = np.array([0, 0, 0])
            torque_thruster = np.array([0, 0, 0])
        
        
        #if the grid is on
        elif controlSystem.switch == 1:
            #fetch thrust and torque from the pixel
            #Using the correct pixel:
            thrust_val, torque_thruster_val = selected_grid.generateTorque(nextPixel) 

            #If grid is off, set torques and thrusts to 0
            thrust = thrust_val * selected_grid.state
            torque_thruster = torque_thruster_val * selected_grid.state


        #set initial condition for integration, i.e. the state vector
        init_cond = satellite_state





        #####################
        # INTEGRATION
        #####################

        #integrate the ODE / propagate the motion
        #args passed as complementary variables to the ODE function
        #tolerances and integrator defined as inputs
        solver = integrate.solve_ivp(coupledMotion,                                 #coupeld ODE
                                    t_span=(t_start, t_end),                        #start, stop time
                                    y0=init_cond,                                   #initial conditions
                                    args=(thrust, torque_thruster,                  #thrust and torque from the thruster
                                          satellite, grids,                         #List of thruster grids
                                          selected_grid,                            #Selected grid for firing
                                          controlSystem,                            #Control system object
                                          environment,                              #Enviornment object
                                          sensor,                                   #Sensor object
                                          manager,                                  #perturbation manager 
                                          iteration),                               #Simulation name
                                    events=controlSystem.cut_off_function,          #Terminal functions to stop integration
                                    dense_output=False,                                
                                    method=integrator,                              #specify integration method
                                    t_eval=None, 
                                    rtol=tols, atol=tols)                           #specify tolerances
        







        ##############
        # BREAKER CODE
        ##############
        #if grid is on
        #ensure the final time point of the integration is different from the starting time point
        if t_start != solver.t[-1] and controlSystem.switch == 1:

            #update the burn time and the ignition log for data storage
            # grid.burnTime[nextPixel] -= (solver.t[-1] - t_start)  #less fuel after firing for FiringTime 
            # selected_grid.ignitionLog[nextPixel] += len(solver.t_events[1])   #how many times a pixel has ignited
            

            #reset breaker --- used to stop simulation if no pixel are ignited
            dataset.breaker = 0
            top_breaker = 0
            

        #if the integration did not actually propagate:
        elif t_start == solver.t[-1]:
            print('\t\t\t\t no advance in time has been made.')
            dataset.breaker +=1 # add to the breaker counter
            top_breaker += 1

 



        ##################################
        # POST-INTEGRATION DATA MANAGEMENT
        ##################################      
        #update state for next integration leg
        satellite_state =  solver.y[:, -1]                      #upadte state vector
        
        #update values
        equinoctial = satellite_state[:6]
        quaternion = satellite_state[6:10]
        rest = satellite_state[10:]
        
        #normalise quaternion if needed
        quat_norm = 1                                           # to normalise: q = np.linalg.norm(quaternion)
        quaternion = quaternion / quat_norm


        #update fuel --take from the state to update fuel grid
        selected_grid.burnTime[nextPixel] = satellite_state[13] #pixel fuel update

        #reform the state vector (quaternion normalised now)
        satellite_state = np.append(equinoctial, quaternion)
        satellite_state = np.append(satellite_state, rest)
        

        # store value of initial start point:
        t_org = t_start
        #update the start time for next integration leg
        t_start = solver.t[-1]
        


        






        ############################   
        #DATA WRITING
        ############################

        #compute norm of angular velocity throught integration:
        omega_norm = (np.square(solver.y[10]) + np.square(solver.y[11]) + np.square(solver.y[12]))**(0.5)


        #write data to textfile:
        for i in range(0, len(solver.t)): 

            #Angular velocity
            dataset.write_value(writeTimetrack, solver.t[i], solver.y[10, i], solver.y[11, i], solver.y[12, i], omega_norm[i])
            
            #Grid1 data
            dataset.write_value(writeGrid1, solver.t[i], grids[0].name, solver.y[14, i], sum(grids[0].burnTime.values()), gridsBurntime(grids))
            
            #Grid2 data
            dataset.write_value(writeGrid2, solver.t[i], grids[-1].name, solver.y[-1, i], sum(grids[-1].burnTime.values()), gridsBurntime(grids))
            
            #Quaternion data
            dataset.write_value(writeQuaterion, solver.t[i], solver.y[6, i], solver.y[7, i], solver.y[8, i], solver.y[9, i]) 
            
            #Orbital data
            dataset.write_value(writeEquinoctial, solver.t[i], solver.y[0, i], solver.y[1, i], solver.y[2, i], solver.y[3, i], solver.y[4, i], solver.y[5, i])
            

        



        ############################   
        #BREAKER CODE
        ############################

        #if the driver has looped 3 times without the itnegration moving forward in time, break
        if dataset.breaker == 3:
            
            print('breaker condition reached: {} integration without t moving.'.format(dataset.breaker))
            print('Something most likely went wrong')
            break


        #store quadrant selection data
        quadName = selected_grid.idealQuadrant
        file2 = open(writepath + '/quadrantSelection'.format(iteration) + '.txt', 'a')
        file2.write('{:.3f}\t{}\t{}\n'.format(solver.t[-1], selected_grid.predictedQuadrant, quadName))
        file2.close()

        #recompute altitude event function. 
        # This number is negative while the pacecraft is above 150 km
        altitude = r_stop(0, satellite_state, t_stop, thrust, nextPixel, torque_thruster, satellite, grids, selected_grid, controlSystem, environment, sensor, iteration)




        ############################   
        #CHECKPOINTING
        ############################

        #prep for checkpointing
        fuelGrid = []
        for i in grids:
            fuelGrid.append(i.burnTime)

        #checkpoint the last value of time and state vector, and grids values and sensor values
        np.savez(writepath + '/checkpoint'.format(iteration), t=t_start, y0=satellite_state, fuel=fuelGrid, sensor_biais=sensor.biais_n)
        print('simulated time [days]', solver.t[-1] / 86400)



    ############################   
    #END OF WHILE LOOP 
    ############################
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

    return None
 
    ################## END ###############

    
    





