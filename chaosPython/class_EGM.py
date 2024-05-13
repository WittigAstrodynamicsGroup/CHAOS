"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_EGM.py




This Python module defines the `EarthGravityModel` class, 
which implements an Earth Gravity Model (EGM) for calculating 
gravitational perturbations acting on a spacecraft.

The `EarthGravityModel` class inherits from the `perturbationInterface` 
class, adhering to the common interface for calculating environmental 
forces and torques.

Key functionalities include:

- **Gravitational Acceleration Calculation:**
  The `computePerturb` method utilizes the provided `egm08` function 
  to calculate the gravitational acceleration vector acting on the 
  spacecraft based on its position and the EGM coefficients.

- **EGM Coefficients Loading:**
  The class constructor (`__init__`) handles loading the EGM 
  coefficients (C and S) from data files located in a specified 
  data directory.
"""

import numpy as np
import os
from .class_perturbationInterface import perturbationInterface
from .EGM_V import matricise, egm08

class EarthGravityModel(perturbationInterface):

    def __init__(self):


        #degree and order of the Model
        self.degree = 2
        self.order = self.degree

        #epoch of the perturbation: 
        self.epoch = None                               ##overriden by perturbationManager


        ###################
        # LOAD COEFFICIENTS
        ###################

        #get data path
        self.dataFilePath = os.path.dirname(__file__) + '/Data'

        ##upload coefficient
        self.C = matricise(2, filename = self.dataFilePath + '/Environment/Gravity/EGM96_n36') #EGM coefficient
        self.S = matricise(3, filename = self.dataFilePath + '/Environment/Gravity/EGM96_n36') #EGM coefficients






    def computePerturb(self, time, satelliteClass):
        """
        Calculates the environmental perturbations acting on the spacecraft for the EGM08.

        **Arguments:**
            - self: Reference to the current object instance.
            - time: Simulation time (s).
            - satelliteClass: Reference to a class holding information about the spacecraft 
                            (e.g., mass, position, coefficients).

        **Returns:**
            tuple(numpy.ndarray): A tuple containing two NumPy ndarrays representing:
                - accVector (km/s^2): Total force vector acting on the spacecraft in the simulation frame due to the modeled environmental effect.
                - torqueVector (km^2/s^2): Total torque vector acting on the spacecraft in the simulation frame (currently set to zero).
            """
        
        torqueVector = np.array([0, 0, 0])                     ##no torque

        accVector =  egm08(time, satelliteClass.eci_position(), self.C, self.S, self.epoch, self.degree, self.order )

        return accVector, torqueVector