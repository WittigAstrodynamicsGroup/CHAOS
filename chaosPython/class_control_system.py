#Python code to hold the class defining the control system/ definition. First step towards a general re-structuring of the code.

'''
Python code where the control_system class is defined. The class effectively will be the cl1 function, re-implemented as a class. 
Through if statement, it will allow for variation / updates of the control system without having to change the source code, which was the case when dealing 
with the cl1 function. 

The control_system class holds the following attributes:

    expected kwargs inputs:
        -closed_loop:   defines whether the a closed loop control system is activated
        -cut_off_function:   Holds the function name that defines the cut-off conditions. Always terminal
        -omegaMax:   Holds the maximum angular velocity allowed for the satellite
        -multipleIgnition:  Defines whether a thruster pixel can be ignited many times or not
        
The methods of the control system are the following:

        -   closedLoop: The method returns the next pixel that should be ignited, based on the current velocity of the satellite.
                        The method assumes that the angular velocity, the fuel left on each pixel, and whether a pixel has been ignited
                        or not are all known variables at any point in time.
                        To predict which pixel should be ignited in order to reduce the anular velocity of the spacecraft, the CL 
                        assumes a perfectly symmetrical body using the nominal thrust at the location of a pixel. The pixel which reduces
                        the velocity the most is chosen. The prediction of the effect includes the amount of fuel on an individual pixel
                        and the re-usability condition. 



'''


from .quaternions import * 
from .support_fun import findCenter
import numpy as np
import random

from .transformations import equi_to_v_cart

class control_system:

    def __init__(self, **kwargs):

        self.closed_loop = kwargs.get('closed_loop', True) #default CS set to use closed loop control
        self.cut_off_function = kwargs.get('cut_off_function', None) # default CS set to not use the cut-off method.INput is the cut-off function needed.
        self.omegaMax = kwargs.get('omegaMax', None)
        self.multipleIgnition = kwargs.get('multipleIgnition', True)
        self.quadrant = kwargs.get('quadrant', False)
        self.ASD = kwargs.get('ASD', None)
        self.BPQS = kwargs.get('BPQS', False)
        self.naive = kwargs.get('naive', False)
        # self.thrustingConeAngle = kwargs.get('thrustingConeAngle', np.pi/2)
        self.ASDcounter = 0
        self.switch = 0

        # self.closed_loop = closed_loop 
        # self.omegaMax = omegaMax
        # self.multipleIgnition = multipleIgnition
        # self.cut_off_function = cut_off_function
        # self.quadrant = quadrant
        # self.ASD = ASD
        
        # self.BPQS = BPQS

        print('create interface to avoid clashes between control systems: quadrant = 1  means cl = 0...')

    #create methods here

    def closedLoop(self, omega_state, grid, satellite):
        '''
        Method to determine the next pixel to be fired, based on 
            -current angular velocity
            -Fuel left on each pixel
            -re-ignition conditions
            -deprecated?
        '''
        omega = omega_state[4:]
        thrust_nominal = grid.thrustNominal
        inertia_inv = satellite.InverseinertiaMatrix
        pixelPositionCopy = grid.pixelPosition.copy() # assign a copy so that the original object is unchanged

        
        prediction_list = []
        r_pos = []
        to_del_ignition = []
        to_del_fuel = []

        ###check re-ignition and fuel conditions
        #check re-igniton conditions:
        if self.multipleIgnition == False:
            for pix in pixelPositionCopy.keys():
                if grid.ignitionLog[pix] > 0: #if pixel has been ignited at all
                    to_del_ignition.append(pix)
                    
        
      
        #delete pixels that have been ignited
        for item in to_del_ignition:
            del pixelPositionCopy[item]
    
        #check fuel left
        for pix in pixelPositionCopy.keys():
            if grid.burnTime[pix] < 1: #if less than 1 sec of burn time remaining,
                to_del_fuel.append(pix) 

        #delete pixels with not enough fuel
        for item in to_del_fuel:
            del pixelPositionCopy[item]
            # print('pixel deleted fuel:', item)

        if len(pixelPositionCopy.keys()) != 0:  #if there are any pixel left

            # #for every pixel on the grid
            for pix in pixelPositionCopy.keys():
                #uses the remaining burn time of a specific pixel
                t_burn = grid.burnTime[pix] 

                estimated_torque = np.cross(pixelPositionCopy[pix], [thrust_nominal, 0, 0]) #assume torque from nominal thrust in x
            
                estimated_omega_effect = np.matmul(inertia_inv, estimated_torque) * t_burn #pure effect of the thruster
                
                estimated_omega = omega + estimated_omega_effect #estimated total velocity of the cubesat
                
                prediction_list.append(np.linalg.norm(estimated_omega)) # list of all predicted outcomes
                r_pos.append(pix) # thruster pixel order
    

            prediction = min(prediction_list) # find the lowest outcome possible
        
            pix_key = r_pos[prediction_list.index(prediction)] # find which pixel corresponds to the lowest omega.
        else:
            pix_key = None
            prediction = None #if no pixels can be chosen, returns None
        return(pix_key)


    def quadrantControl_truth(self, omega_state, grid, satellite):
        
        '''
        Method to determine the next pixel to be selected in a quadrant configuration
        based on
        -fuel left
        -current angular velocity
        -re-ignition is assumed true
        '''
        omega = omega_state
        
        
        thrust_nominal = grid.thrustNominal
        axis = grid.thrustAxis
        t_burn = grid.QIT #TO CHANGE
        inertia_inv = satellite.InverseinertiaMatrix
        pixelPositionCopy = grid.pixelPosition.copy() # assign a copy so that the original object is unchanged
        #01, 11, 10, 00
        quadrant_centers = [] #[np.array([0.05, 0.03, 0.03]), np.array([0.05, -0.03, 0.03]), np.array([0.05, -0.03, -0.03]), np.array([0.05, 0.03, -0.03])]
        
        prediction_list = []
        to_del_fuel = []
        q_pos = []
        q_name = []

        #check fuel left
        for pix in pixelPositionCopy.keys():
            if grid.burnTime[pix] < 1: #if less than 1 sec of burn time remaining,
                to_del_fuel.append(pix) 

        #delete pixels with not enough fuel
        for item in to_del_fuel:
            del pixelPositionCopy[item]
            
        if len(pixelPositionCopy.keys()) != 0:  #if there are any pixel left
            # split pixel with fuel in quadrants
            quadrantDicts = grid.quadrant(pixelPositionCopy)
            quadrants = list(quadrantDicts)

            #get quadrant centre and names:
            name = []
            for quad in range(len(quadrants)):
                quadrant_centers.append(findCenter(quadrants[quad]))

                name.append(str(quad))

            #Assess the effect of a quadrant if it has fuel left:
            for i in range(len(quadrants)):
                #if there is fuel left
                if len(quadrants[i]) > 0:
                    
                    estimated_torque = np.cross(quadrant_centers[i], thrust_nominal * axis)
                                
                    estimated_omega_effect = np.matmul(inertia_inv, estimated_torque) * t_burn #pure effect of the thruster

                    estimated_omega = omega + estimated_omega_effect #estimated total velocity of the cubesat
                    
                    prediction_list.append(np.linalg.norm(estimated_omega)) # list of all predicted outcomes
                    q_pos.append(quadrants[i])
                    q_name.append(name[i])
            #find the quadrant with minimum effect
            prediction = min(prediction_list)
            index = prediction_list.index(prediction)
            min_effect_quadrant = q_pos[index]
            grid.idealQuadrant = q_name[index]
            pixKey = q_name[index]
        else:
            pixKey = None
        return pixKey




    def quadrantControl(self, omega_state, grid, satellite, sensor):
        
        '''
        Method to determine the next pixel to be selected in a quadrant configuration
        based on
        -fuel left
        -current angular velocity
        -re-ignition is assumed true
        sensor input
        '''
        #read sensor measurement, default to true reading if none
        omega = omega_state[4:]
        # omega = np.array([0, 0, 0])
        if hasattr(sensor, 'MEMSomega'):
            omega = sensor.MEMSomega
        
        thrust_nominal = grid.thrustNominal
        axis = grid.thrustAxis
        # print('axis', axis)
        t_burn = grid.QIT #TO CHANGE
        inertia_inv = satellite.InverseinertiaMatrix
        pixelPositionCopy = grid.pixelPosition.copy() # assign a copy so that the original object is unchanged
        #01, 11, 10, 00
        quadrant_centers = []# [np.array([0.05, 0.03, 0.03]), np.array([0.05, -0.03, 0.03]), np.array([0.05, -0.03, -0.03]), np.array([0.05, 0.03, -0.03])]
        
        prediction_list = []
        to_del_fuel = []
        q_pos = []
        q_name = []

        #check fuel left
        for pix in pixelPositionCopy.keys():
            if grid.burnTime[pix] < 1: #if less than 1 sec of burn time remaining,
                to_del_fuel.append(pix) 

        #delete pixels with not enough fuel
        for item in to_del_fuel:
            del pixelPositionCopy[item]
            
        if len(pixelPositionCopy.keys()) != 0:  #if there are any pixel left
            # split pixel with fuel in quadrants
            quadrantDicts = grid.quadrant(pixelPositionCopy)
            quadrants = list(quadrantDicts)
            #get quadrant centre & names:
            name = []
            for quad in range(len(quadrants)):
                quadrant_centers.append(findCenter(quadrants[quad]))
                name.append(str(quad))

            #Assess the effect of a quadrant if it has fuel left:
            for i in range(len(quadrants)):
                #if there is fuel left
                if len(quadrants[i]) > 0:
                    
                    estimated_torque = np.cross(quadrant_centers[i], thrust_nominal * axis)
                                
                    estimated_omega_effect = np.matmul(inertia_inv, estimated_torque) * t_burn #pure effect of the thruster

                    estimated_omega = omega + estimated_omega_effect #estimated total velocity of the cubesat
                    
                    prediction_list.append(np.linalg.norm(estimated_omega)) # list of all predicted outcomes
                    q_pos.append(quadrants[i])
                    q_name.append(name[i])
            #find the quadrant with minimum effect
            prediction = min(prediction_list)
            index = prediction_list.index(prediction)
            min_effect_quadrant = q_pos[index]
            grid.predictedQuadrant = q_name[index]
            
            # randomly choose pixel to ignite
            pix_key = random.choice(list(min_effect_quadrant.keys()))


        else:
            pix_key = None

        return pix_key
            
    def bestPixelQuadrant(self, omega_state, grid, satellite):
        omega = omega_state[4:]
        thrust_nominal = grid.thrustNominal
        axis = grid.thrustAxis
        t_burn = grid.QIT #uses QIT for pixel burn time
        inertia_inv = satellite.InverseinertiaMatrix
        pixelPositionCopy = grid.pixelPosition.copy() # assign a copy so that the original object is unchanged
        
        prediction_list = []
        to_del_fuel = []
        r_pos = []
        

        #check fuel left
        for pix in pixelPositionCopy.keys():
            if grid.burnTime[pix] < 1: #if less than 1 sec of burn time remaining,
                to_del_fuel.append(pix) 

        #delete pixels with not enough fuel
        for item in to_del_fuel:
            del pixelPositionCopy[item]
            
        if len(pixelPositionCopy.keys()) != 0:  #if there are any pixel left
            # split pixel with fuel in quandrants
            I, II, III, IV = grid.quadrant(pixelPositionCopy)
            quadrants = [I, II, III, IV]

            #Assess the effect of each remaining pixel:
            if len(pixelPositionCopy.keys()) != 0:  #if there are any pixel left

                # #for every pixel on the grid
                for pix in pixelPositionCopy.keys():
                    
                    estimated_torque = np.cross(pixelPositionCopy[pix], thrust_nominal * axis) #assume torque from nominal thrust in x
                
                    estimated_omega_effect = np.matmul(inertia_inv, estimated_torque) * t_burn #pure effect of the thruster
                    
                    estimated_omega = omega + estimated_omega_effect #estimated total velocity of the cubesat
                    
                    prediction_list.append(np.linalg.norm(estimated_omega)) # list of all predicted outcomes
                    r_pos.append(pix) # thruster pixel order
        
                prediction = min(prediction_list) # find the lowest outcome possible
                index = prediction_list.index(prediction) #posiiton in list
                min_effect_pixel = r_pos[index] #pixel key for min effect
                
                #determine which quadrant the pixel is in
                for Q in quadrants:
                    if min_effect_pixel in Q:
                        min_effect_quadrant = Q
                # randomly choose pixel to ignite
            pix_key = random.choice(list(min_effect_quadrant.keys()))

            if prediction > np.linalg.norm(omega):
                self.ASDcounter += 1
            else:
                self.ASDcounter -= 1

            #Anti-Spin Degeneration algorithm
            if self.ASD != None:
                # print('ASD', self.ASDcounter)
                # If omega has increased for ASD amount of times in a row, terminate
                if self.ASDcounter > self.ASD:
                    print('ASD termination')
                    pix_key = None
        else:
            pix_key = None

        return pix_key


    def naiveControl(self, omega_state, grid, satellite):

        """
        Naive approach to operating the thruster. No quadrant separation, instead, the 
        pixel is simply randomly ignited, and fires for X amount of time.
        """
        pixelPositionCopy = grid.pixelPosition.copy() # assign a copy so that the original object is unchanged
        
        to_del_fuel = []

        

        #check fuel left
        for pix in pixelPositionCopy.keys():
            if grid.burnTime[pix] < 1: #if less than 1 sec of burn time remaining,
                to_del_fuel.append(pix) 

        #delete pixels with not enough fuel
        for item in to_del_fuel:
            del pixelPositionCopy[item]
            
        if len(pixelPositionCopy.keys()) != 0: 
            
            #Randomly choses a pixel from the allowable pixels (fuel left)
            pix_key = random.choice(list(pixelPositionCopy.keys()))
            
        else:
            pix_key = None

        return pix_key

    def velocity_tracking(self, satellite_state, grid ):
        # function to return the angle relative to the orbital velocity direction
        axis = grid.thrusterOrientation
        omega_state = satellite_state[6:]
        quat = omega_state[:4]
        #rotate thruster axis in inertial frame 
        axis_in  = body_to_inertial(axis, quat)
        axis_in  = axis_in / np.linalg.norm(axis_in)
        #compute velocity direction
        v = equi_to_v_cart(satellite_state[:6])
        v_unit = v / np.linalg.norm(v)

        #compute thrust / velocity angle 
        coss = np.dot(v_unit, axis_in)
        angle = np.arccos(coss) 

        return angle

    def decidePixel(self, *args):
        '''
        Function called by wrapper. The decidePixel function simply calls other function based on the control system parameters
        The args to be passed are the argument to the function called by the decidePixel function
        '''
        if self.closed_loop == True:
            pix_key = self.closedLoop(args[0], args[1], args[2])
            return pix_key

        elif self.quadrant == True:
            self.multipleIgnition = True
            
            #set max angular velocity 
            #OR t_avg to max angular velocity #in grid
            pix_key = self.quadrantControl(args[0], args[1], args[2], args[3])
            return pix_key
        
        elif self.BPQS == True:
            self.multipleIgnition = True

            #et the control loop
            pix_key = self.bestPixelQuadrant(args[0], args[1], args[2])
            return pix_key


        elif self.naive == True:
            self.multipleIgnition = True
            pix_key = self.naiveControl(args[0], args[1], args[2])
            return pix_key