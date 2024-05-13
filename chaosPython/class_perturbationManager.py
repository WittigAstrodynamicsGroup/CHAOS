"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_perturbationManager.py



This Python module defines the `perturbationManager` class, which manages 
environmental perturbations acting on a spacecraft during orbit simulations.

The `perturbationManager` class:

- Stores central body information (`centralBody`) and reference epoch (`epoch`).
- Maintains a list of registered force models (`forceModelsList`).
- Provides methods for:
    - Initializing the simulation environment by propagating the central body 
      and epoch information to all registered force models (`initialSetUp`).
    - Calculating the total environmental perturbations acting on the spacecraft 
      by combining the effects from various force models (`callForceModels`). 
      This function iterates through the force models, accumulates their 
      contributions to acceleration and torque, and returns the total 
      force and torque vectors.

**Key Concepts:**

- **Force Models:** Specific classes representing environmental perturbations, using the interface 
    from perturbationInterface (class_perturbationInterface.py)
"""

import numpy as np





class perturbationManager:





    ####################
    # INSTANTIATION
    ####################

    def __init__(self):



        self.epoch = np.array([1, 1, 2000])                 #simulation epoch
        self.centralBody = 'earth'                          #simulation cenrtal body             
        self.mu = 398600.4418                               #gravitational parameter, km-5-vallado
        self.bodyRadius = 6378.1363                         #Primary body radius

        self.forceModelsList = []                           #list where the perturbation class are stored

        

        #constants
        ##
        #planetary constants
        # self.R_E = 6378.1363 # km
        # self.R_Sun = 696000 #km
        # self.mu_sun = 1.327122e11#1.32712428e11 #km^3 /s^2 gmat:132712440017.99
        # self.mu_moon = 4.902800305555e3#4902.799 #km^3 /s2 GMAT: 4902.8005821478





    ####################
    # METHODS
    ####################

    def initialSetUp(self): 

        """
        Propagates central body and epoch information to all registered force models.

        This function ensures consistency in the environmental conditions 
        experienced by the spacecraft during the simulation. It iterates through 
        the list of force models (`self.forceModelsList`) and sets the epoch and centralBody atributes.
        """

        for cls in self.forceModelsList:
            cls.centralBody = self.centralBody
            cls.epoch = self.epoch
            cls.R_E = self.bodyRadius



    def callForceModels(self, time, satelliteClass):

        """
        Calculates the total environmental perturbations acting on the spacecraft by combining forces from various models.

        This function iterates through a list of force models (`self.forceModelsList`) 
        and accumulates the environmental accelerations and torques acting on the 
        spacecraft for the given simulation time (`time`) and spacecraft properties 
        (`satelliteClass`).


        **Returns:**
            tuple(numpy.ndarray): A tuple containing two NumPy ndarrays representing:
                - accVec (km/s^2): Total force vector acting on the spacecraft in the simulation frame due to all environmental perturbations considered by the provided force models.
                - torqueVec (Nm): Total torque vector acting on the spacecraft in the simulation frame due to all environmental perturbations considered by the provided force models.
        """

        #define acceleration and torque vectors
        accVec = np.array([0., 0., 0.])
        torqueVec = np.array([0., 0., 0.])

        #for each provided perturbation
        for f in self.forceModelsList:
            acc, torque = f.computePerturb(time, satelliteClass)

            accVec += acc
            torqueVec += torque

        return accVec, torqueVec
