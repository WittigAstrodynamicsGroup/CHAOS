#Python code to hold the class defining the control system/ definition. First step towards a general re-structuring of the code.

"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_control_system.py


**Control System Module**

This module defines the `control_system` class, which implements 
thruster control functionalities for a spacecraft simulation. 

The class provides methods for:

* **Quadrant control strategies:**
    * 'quadrantControl_truth': Determines 
      the optimal quadrant to fire based on the true state for maximum 
      control effectiveness.

    * 'quadrantControl': Similar to `quadrantControl_truth`, 
    but uses sensor readings (if available) 
    instead of the true state.

* **Velocity tracking:**
    * `velocity_tracking`: Calculates the angle 
      between the thruster axis and the spacecraft's orbital velocity 
      direction.

* **Control function selection:**
    * `decidePixel(self, *args)`: Calls the user-defined control function 
      provided during initialization.

**External Dependencies:**

* `quaternions.py` for quaternion operations.
* `support_fun.py` for utility functions (findCenter).
* `numpy` for numerical computations.
* `random` for random number generation (quadrantControl).
* `transformations.py` for state vector conversion (equi_to_v_cart).

"""


from .quaternions import * 
from .support_fun import findCenter
import numpy as np
import random

from .transformations import equi_to_v_cart

class control_system:

    def __init__(self, **kwargs):

        self.CSfunction = kwargs.get("CSfunction", None)            #Specify the Control system function
        self.switch = 0                                             #Cube-de-ALPS switch. 
                                                                    #---if 0, no thruster firing will occur. 
        self.cut_off_function = kwargs.get('cut_off_function', None)





    #####################
    #DEFINE METHODS HERE
    #####################


    def quadrantControl_truth(self, omega_state, grid, satellite):
        
        """
        Determines the quadrant to be fired in a quadrant control strategy, 
        based on the true state of the spacecraft.

        This function determines the most appropriate pixel to fire next within 
        a thruster configuration divided into quadrants. It considers the following
        factors:

        * **Remaining fuel:** Pixels with less than 1 second of burn time remaining 
        are excluded.
        * **Current angular velocity (`omega_state`):** The function uses the 
        actual state to determine the relevant quadrant for control.
        * **Re-ignition:** Re-ignition is assumed to be always possible.

        The function performs the following steps:

        1. **Fuel check:** It identifies and removes pixels from consideration 
        that have insufficient remaining fuel.

        2. **Quadrant division:** It splits the remaining pixels with fuel into 
        their respective quadrants based on the thruster configuration (`grid`).

        3. **Quadrant center and name extraction:** It calculates the center 
        position and assigns a name (string representation of the quadrant 
        number) for each quadrant with remaining pixels.

        4. **Effect prediction:** For each quadrant with fuel, it estimates the 
        resulting angular velocity after firing the entire quadrant using a 
        simplified dynamic model. This estimation considers the nominal thrust, 
        axis of thrust, quadrant center, inertia of the spacecraft, and 
        nominal burn time.

        5. **Optimal quadrant selection:** It selects the quadrant that 
        minimizes the predicted final angular velocity magnitude. This 
        corresponds to the quadrant that would have the most significant 
        corrective effect on the current state.

        6. **name assignement:** It assigns the name of the selected quadrant 
        to the control system class

        Args:
            self: Reference to the `Satellite` object.
            omega_state: Current angular velocity of the spacecraft (rad/s).
            grid: Reference to the thruster configuration object.
            satellite: Reference to the `Satellite` object (used for inertia).

        Returns:
            str: The key (name) of the selected quadrant, or `None` if no pixels 
                have sufficient fuel.
        """

        
        omega = omega_state                                 #Use real state to determine quadrant
        
        
        thrust_nominal = grid.thrustNominal                 #nominal thrust  
        axis = grid.thrustAxis                              #thrust direction
        t_burn = grid.QIT                                   #Nominal time for which pixel is fired
        inertia_inv = satellite.InverseinertiaMatrix        #Inverse inertia matrix
        pixelPositionCopy = grid.pixelPosition.copy()       # assign a copy so that the original object is unchanged

        #define useful lists
        quadrant_centers = []   #will hold position of centres for each quadrant
        prediction_list = []    #will hold predicted final angular vel after firing 
        to_del_fuel = []        #holds the pixels that have run out of fuel
        q_pos = []              #Holds the quadrant position in loop
        q_name = []             #holds the quadrant name 



        #check fuel left
        for pix in pixelPositionCopy.keys():
            if grid.burnTime[pix] < 1:              #if less than 1 sec of burn time remaining, consider pixel as burnt out
                to_del_fuel.append(pix) 

        #delete pixels with not enough fuel
        for item in to_del_fuel:
            del pixelPositionCopy[item]
            

        #if there are any pixel left
        if len(pixelPositionCopy.keys()) != 0:  

            # split pixel with fuel in quadrants
            quadrantDicts = grid.quadrant(pixelPositionCopy)
            quadrants = list(quadrantDicts)                             #make into a list

            #get quadrant centres and names:
            name = []
            for quad in range(len(quadrants)):                          #for each quadrant
                quadrant_centers.append(findCenter(quadrants[quad]))    #find the centre   

                name.append(str(quad))                                  #get name in a string 



            #Assess the effect of a quadrant if it has fuel left:
            for i in range(len(quadrants)):

                #if there is fuel left
                if len(quadrants[i]) > 0:                               
                    

                    ###simplified dynamics to estimate the angular velocity after firing
                    #assume nominal thrust
                    #firing from centre of quadrant
                    

                    #estimated torque
                    estimated_torque = np.cross(quadrant_centers[i], thrust_nominal * axis)

                    #estimated delta omega            
                    estimated_omega_effect = np.matmul(inertia_inv, estimated_torque) * t_burn #pure effect of the thruster


                    #estimated final omega
                    estimated_omega = omega + estimated_omega_effect #estimated total velocity of the cubesat
                    

                    #append to list 
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
        
        """
        Determines the quadrant to be fired in a quadrant control strategy, 
        based on the sensor reading of the state of the spacecraft.

        This function determines the most appropriate pixel to fire next within 
        a thruster configuration divided into quadrants. It considers the following
        factors:

        * **Remaining fuel:** Pixels with less than 1 second of burn time remaining 
        are excluded.
        * **Current angular velocity (`omega_state`):** The function uses the 
        actual state to determine the relevant quadrant for control.
        * **Re-ignition:** Re-ignition is assumed to be always possible.

        The function performs the following steps:

        1. **Fuel check:** It identifies and removes pixels from consideration 
        that have insufficient remaining fuel.

        2. **Quadrant division:** It splits the remaining pixels with fuel into 
        their respective quadrants based on the thruster configuration (`grid`).

        3. **Quadrant center and name extraction:** It calculates the center 
        position and assigns a name (string representation of the quadrant 
        number) for each quadrant with remaining pixels.

        4. **Effect prediction:** For each quadrant with fuel, it estimates the 
        resulting angular velocity after firing the entire quadrant using a 
        simplified dynamic model. This estimation considers the nominal thrust, 
        axis of thrust, quadrant center, inertia of the spacecraft, and 
        nominal burn time.

        5. **Optimal quadrant selection:** It selects the quadrant that 
        minimizes the predicted final angular velocity magnitude. This 
        corresponds to the quadrant that would have the most significant 
        corrective effect on the current state.

        6. **Pixel selection:** It selects a random pixel within the selected quadrant

        Args:
            self: Reference to the `Satellite` object.
            omega_state: Current angular velocity of the spacecraft (rad/s).
            grid: Reference to the thruster configuration object.
            satellite: Reference to the `Satellite` object (used for inertia).

        Returns:
            str: The key (name) of the selected pixel, or `None` if no pixels 
                have sufficient fuel.
        """
        #read sensor measurement, default to true state if no measurement available
        omega = omega_state[4:]
        # omega = np.array([0, 0, 0])
        if hasattr(sensor, 'MEMSomega'):
            omega = sensor.MEMSomega
        
        thrust_nominal = grid.thrustNominal                 #nominal thrust  
        axis = grid.thrustAxis                              #thrust direction
        t_burn = grid.QIT                                   #Nominal time for which pixel is fired
        inertia_inv = satellite.InverseinertiaMatrix        #Inverse inertia matrix
        pixelPositionCopy = grid.pixelPosition.copy()       # assign a copy so that the original object is unchanged

        #define useful lists
        quadrant_centers = []   #will hold position of centres for each quadrant
        prediction_list = []    #will hold predicted final angular vel after firing 
        to_del_fuel = []        #holds the pixels that have run out of fuel
        q_pos = []              #Holds the quadrant position in loop
        q_name = []             #holds the quadrant name 



        #check fuel left
        for pix in pixelPositionCopy.keys():
            if grid.burnTime[pix] < 1:              #if less than 1 sec of burn time remaining, consider pixel as burnt out
                to_del_fuel.append(pix) 

        #delete pixels with not enough fuel
        for item in to_del_fuel:
            del pixelPositionCopy[item]
            


         #if there are any pixel left
        if len(pixelPositionCopy.keys()) != 0:  

            # split pixel with fuel in quadrants
            quadrantDicts = grid.quadrant(pixelPositionCopy)
            quadrants = list(quadrantDicts)                             #make into a list

            #get quadrant centres and names:
            name = []
            for quad in range(len(quadrants)):                          #for each quadrant
                quadrant_centers.append(findCenter(quadrants[quad]))    #find the centre   

                name.append(str(quad))                                  #get name in a string 



            #Assess the effect of a quadrant if it has fuel left:
            for i in range(len(quadrants)):

                #if there is fuel left
                if len(quadrants[i]) > 0:                               
                    

                    ###simplified dynamics to estimate the angular velocity after firing
                    #assume nominal thrust
                    #firing from centre of quadrant
                    

                    #estimated torque
                    estimated_torque = np.cross(quadrant_centers[i], thrust_nominal * axis)

                    #estimated delta omega            
                    estimated_omega_effect = np.matmul(inertia_inv, estimated_torque) * t_burn #pure effect of the thruster


                    #estimated final omega
                    estimated_omega = omega + estimated_omega_effect #estimated total velocity of the cubesat
                    

                    #append to list 
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
            



    

    def velocity_tracking(self, satellite_state, grid ):


        """
        Calculates angle between thruster axis and orbital velocity.


        Args:
            satellite_state: Spacecraft state vector (position, velocity).
            grid: Thruster configuration object.

        Returns:
            float: Angle between thruster and orbital velocity (rad).
        """


        # Extract thruster axis
        axis = grid.thrusterOrientation

        # Extract state
        omega_state = satellite_state[6:]
        quat = omega_state[:4]

        #rotate thruster axis in inertial frame 
        axis_in  = body_to_inertial(axis, quat)
        axis_in  = axis_in / np.linalg.norm(axis_in)

        #compute velocity direction
        v = equi_to_v_cart(satellite_state[:6])
        v_unit = v / np.linalg.norm(v)

        #compute thrust - orbital velocity angle 
        coss = np.dot(v_unit, axis_in)
        angle = np.arccos(coss) 

        return angle




    def decidePixel(self, *args):
        '''
        Function called by wrapper. The decidePixel function simply calls the user-defined control function.
        The args to be passed are the argument to the function called by the decidePixel function
        '''

        pix_key = self.CSfunction(args[0], args[1], args[2], args[3])

        
        return pix_key
        
