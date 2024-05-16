
"""


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_SRP.py





This Python module defines two classes for modeling solar radiation pressure (SRP) 
perturbations on a spacecraft during simulations:

- **LowFidelitySRP:**
  This class implements a low-fidelity SRP model. It calculates the 
  acceleration acting on the spacecraft due to SRP based on averaged and 
  constant satellite properties (e.g., area, coefficient of reflection). 
  The `LFsolarRadiationPressure` function is used for this calculation.

- **HighFidelitySRP:**
  This class implements a high-fidelity SRP model. It calculates the 
  total force and torque vectors acting on the spacecraft by considering 
  its detailed geometry (using facet data from an STL file). 
  The `solarRadiationPressure` function is used for this 
  computation, which sums the forces acting on each facet to obtain 
  the total force and torque vectors.

Both classes inherit from the `perturbationInterface` class, 
ensuring a consistent approach for handling environmental perturbations 
within the simulation framework.

**Assumptions**

- Both models assume the Sun is the primary source of SRP.


"""







import numpy as np
from .class_perturbationInterface import perturbationInterface
from .perturbationModelling import solarRadiationPressure, LFsolarRadiationPressure
from .SPICEroutine import spiceTargPos






class LowFidelitySRP(perturbationInterface):

        
    """
    Implementation of the low fidelity SRP perturbation model, using averaged and constant satellite properties
    """



    ####################
    # INSTANTIATION
    ####################

    def __init__(self):
        


        #handled by perturbationManager:
        self.epoch = None                               #overridden by perturbationManager
        self.centralBody = 'earth'                      #defaults to earth      
        self.R_E = 6378.1363                            #Primary body radius


        #constants
        self.R_Sun = 696000                             # sun radius, km  
        self.lightSpeed = 299792458                     # speed of light, m/s

        
        
    ####################
    # METHODS
    ####################


    def computePerturb(self, time, satelliteClass):
        """
        Calculates the environmental perturbations acting on the spacecraft due to solar radiation pressure (SRP).

        **Arguments:**
            - self: Reference to the current object instance.
            - time: Simulation time (seconds).
            - satelliteClass: Reference to a class holding information about the spacecraft 
                            (e.g., mass, position, coefficients).

        **Returns:**
            tuple(numpy.ndarray): A tuple containing two NumPy ndarrays representing:
                - accVector (km/s^2): Total force vector acting on the spacecraft in the simulation frame due to solar radiation pressure.
                - torqueVector (km^2/s^2): Total torque vector acting on the spacecraft in the simulation frame (currently set to zero).
        """




        torqueVector = np.array([0, 0, 0])              #no torques 


        #compute sun position
        vec_sun = spiceTargPos(time, self.epoch, 'sun', self.centralBody)


        #compute SRP acceleration
        accVector = LFsolarRadiationPressure(satelliteClass.eci_position(),
                                             vec_sun, 
                                             satelliteClass.area,
                                             satelliteClass.mass,
                                             satelliteClass.Cr,             #coefficient of reflection
                                             self.lightSpeed, 
                                             self.R_E, 
                                             self.R_Sun)



        return accVector, torqueVector








class HighFidelitySRP(perturbationInterface):

        
    """
    Implementation of the High fidelity SRP perturbation model, taking the satellite geometry into account
    """



    ####################
    # INSTANTIATION
    ####################

    def __init__(self):
        


        #handled by perturbationManager:
        self.epoch = None                               #overridden by perturbationManager
        self.centralBody = 'earth'                      #defaults to earth
        self.R_E = 6378.1363                            #Primary body radius




        #constants
        self.R_Sun = 696000                             # sun radius, km  
        self.lightSpeed = 299792458                     # speed of light, m/s
        
        
        
    ####################
    # METHODS
    ####################


    def computePerturb(self, time, satelliteClass):

        """
        Calculates the environmental perturbations acting on the spacecraft due to solar radiation pressure (SRP).

        This function implements a high-fidelity SRP model for this derived class 
        of `perturbationInterface`. It calculates the total force and torque vectors 
        acting on the spacecraft's gemoetry due to solar radiation pressure.

     
        **Arguments:**
            - self: Reference to the current object instance.
            - time: Simulation time (seconds).
            - satelliteClass: Reference to a class holding information about the spacecraft 
                            (e.g., mass, position, coefficients).

        **Returns:**
            tuple(numpy.ndarray): A tuple containing two NumPy ndarrays representing:
                - accVector (km/s^2): Total force vector acting on the spacecraft in the simulation frame due to solar radiation pressure.
                - torqueVector (Nm): Total torque vector acting on the spacecraft in the simulation frame due to solar radiation pressure.
        """



        #get position and velocity vectors
        r_eci = satelliteClass.eci_position()

        #get satellite data:
        #STL data --each facet
        arr_area = satelliteClass.areas
        arr_normal = satelliteClass.normals
        arr_Cr = satelliteClass.Crs
        COMpos =  satelliteClass.positions - satelliteClass.centreOfMass      

        #properties:
        m = satelliteClass.mass
        quaternion = satelliteClass.quaternion

        #compute sun position
        vec_sun = spiceTargPos(time, self.epoch, 'sun', self.centralBody)


        #compute SRP acceleration
        accVector, torqueVector = solarRadiationPressure(r_eci, 
                                           vec_sun, 
                                           m, 
                                           quaternion, 
                                           arr_area, 
                                           arr_normal, 
                                           COMpos, 
                                           arr_Cr, 
                                           self.lightSpeed,
                                           self.R_E,
                                           self.R_Sun)



        return accVector, torqueVector
