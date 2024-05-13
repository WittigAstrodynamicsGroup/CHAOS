"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_nBody.py




This Python module defines the `nBody` class, which implements a 
point-mass gravitational perturbation model for a spacecraft simulation.

The `nBody` class inherits from the `perturbationInterface` class, 
adhering to a common interface for calculating environmental forces 
and torques acting on a spacecraft. It specifically models the 
gravitational attraction of another celestial body (represented 
as a point mass) on the spacecraft.

**Key functionalities include:**

- **Point-Mass Gravity Calculation:**
  The `computePerturb` method calculates the gravitational acceleration 
  vector acting on the spacecraft due to the point mass's gravity. 
  It utilizes the `nGravity` function and considers the spacecraft's 
  position, the point mass's position obtained using `spiceTargPos` 
  (SPICE library), and the point mass's gravitational 
  parameter (`self.mu`).
"""

import numpy as np


from .class_perturbationInterface import perturbationInterface
from .perturbationModelling import nGravity
from .SPICEroutine import spiceTargPos




class nBody(perturbationInterface):

    ####################
    # INSTANTIATION
    ####################

    def __init__(self, gravitationalParameter, SPICEname):
        

        #inputs
        self.mu = gravitationalParameter                # Gravitational parameter
        self.name = SPICEname                           # String for the name of the body
        
        #handled by perturbationManager
        self.centralBody = "earth"                            #defaults to EARTH
        self.epoch = None                               #overridden by perturbationManager        



    ####################
    # METHODS
    ####################



    def computePerturb(self, time, satelliteClass):

        """
        Calculates the environmental perturbations acting on the spacecraft due to another point mass.

        This function implements the perturbation model for this derived class 
        of `perturbationInterface`. It calculates the total force vector acting 
        on the spacecraft due to the gravitational attraction of another point mass.

        **Arguments:**
            - self: Reference to the current object instance.
            - time: Simulation time (specific format depends on the `spiceTargPos` function).
            - satelliteClass: Reference to a class holding information about the spacecraft 
                                (e.g., mass, position, coefficients).

        **Returns:**
            tuple(numpy.ndarray): A tuple containing two NumPy ndarrays representing:
                - accVector (km/s^2): Total force vector acting on the spacecraft in the simulation centralBody due to the modeled environmental effect.
                - torqueVector (Nm): Total torque vector acting on the spacecraft in the simulation centralBody (currently set to zero).
        """
        
        torqueVector = np.array([0, 0, 0])              #no torque
        
        
        #compute position of body:
        vecBody = spiceTargPos(time, self.epoch, self.name, self.centralBody)

        #compute acceleration vector
        accVector = nGravity(satelliteClass.eci_position(), vecBody, self.mu)

        return accVector, torqueVector