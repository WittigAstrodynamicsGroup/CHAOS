"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_atmosphericDrag.py







This Python module defines two classes for modeling atmospheric drag 
perturbations on a spacecraft during simulations:

- **LowFidelityDrag:**
  This class implements a low-fidelity drag model. It calculates the 
  drag acceleration acting on the spacecraft based on averaged and 
  constant satellite properties (e.g., area, drag coefficient). 
  The `LFdrag` function is used for this calculation.

- **HighFidelityDrag:**
  This class implements a high-fidelity drag model. It calculates the 
  aerodynamic acceleration and torques acting on the spacecraft by 
  considering its detailed geometry (using facet data from 
  an STL file). The `drag` function is used for this 
  computation, which sums the forces acting on each facet to 
  obtain the total force and torque vectors.

Both classes inherit from the `perturbationInterface` class, 
ensuring a consistent approach for handling environmental perturbations 
within the simulation framework.


**Common Assumptions and Considerations:**

- Both models assume the Earth's atmosphere is the primary source 
  of drag.


"""

from .perturbationModelling import LFdrag, drag, density
from .class_perturbationInterface import perturbationInterface
from .EnvironmentModelling import interpF10Index_ext, atmTemperature
from .transformations import julianDate
from .SPICEroutine import spiceTargPos
import numpy as np
import os



class LowFidelityDrag(perturbationInterface):

    """
    Low fidelity drag model, computing drag acceleration based on averaged and constant satellite properties.
    """
    ####################
    # INSTANTIATION
    ####################

    def __init__(self):

        ###################
        # LOAD COEFFICIENTS
        ###################

        #get data path
        self.dataFilePath = os.path.dirname(__file__) + '/Data'
        
        #space weather file
        self.interpf10, self.interpctr81 = interpF10Index_ext(self.dataFilePath + '/Environment/Atmosphere/SpaceWeather-test.txt')


        #handled by perturbationManager:
        self.epoch = None                               #overridden by perturbation manager
        self.centralBody = 'earth'                      #defaults to earth
        self.R_E = 6378.1363                            #defaults to Earth

        #constants:
        self.EarthFlattening = 1 / 298.257                  #flattening coefficeint of the Earth (Vallado)
        self.omega_E = np.array([0, 0, 0.0000729211585530]) #earth rotation vector, rad/s


    ####################
    # METHODS
    ####################

    def computePerturb(self, time, satelliteClass):   
        

        torqueVector = np.array([0, 0, 0])              #no torques


        #get position and velocity vectors
        r_eci = satelliteClass.eci_position()
        v_eci = satelliteClass.eci_velocity()


        #transform into epoch to julian date
        fractionalTime =  time / 86400
        day, month, year = self.epoch
        JD = julianDate(fractionalTime, day, month, year)
        

        #extract solar activity
        F10 = self.interpf10(JD)
        F10_avg = self.interpctr81(JD)

        #sun position:
        vec_sun = spiceTargPos(time, self.epoch, 'sun', self.centralBody)

        #atmospheric temperature
        T = atmTemperature(F10, F10_avg, vec_sun, r_eci, self.EarthFlattening)


        #call density functions
        r = np.linalg.norm(r_eci)
        h = (r - self.R_E) # km
        rho = density(h, T)

        #compute acceleration vector
        accVector = LFdrag(r_eci, v_eci, rho, satelliteClass.area, satelliteClass.Cd, satelliteClass.mass, self.omega_E)


        return accVector, torqueVector












class HighFidelityDrag(perturbationInterface):
    """
    High fidelity drag model, computing aerodynamic acceleration and torques based on spacecraft geometry
    """

    ####################
    # INSTANTIATION
    ####################

    def __init__(self):

        ###################
        # LOAD COEFFICIENTS
        ###################

        #get data path
        self.dataFilePath = os.path.dirname(__file__) + '/Data'
        
        #space weather file
        self.interpf10, self.interpctr81 = interpF10Index_ext(self.dataFilePath + '/Environment/Atmosphere/SpaceWeather-test.txt')


        #handled by perturbationManager:
        self.epoch = None                               #overridden by perturbation manager
        self.centralBody = 'earth'                      #defaults to earth


        #constants:
        self.EarthFlattening = 1 / 298.257                  #flattening coefficeint of the Earth (Vallado)
        self.omega_E = np.array([0, 0, 0.0000729211585530]) #earth rotation vector, rad/s
        self.R_E = 6378.1363 # km






    ####################
    # METHODS
    ####################

    def computePerturb(self, time, satelliteClass):   
        
        """
        Calculates the environmental perturbations acting on the spacecraft due to atmospheric drag.

        This function implements the perturbation model for this derived class 
        of `perturbationInterface`. It calculates the total force vector acting 
        on the spacecraft due to aerodynamic drag forces caused by the atmosphere.

        **Process:**

        1. **Satellite Data Acquisition:**
        - Retrieves satellite properties from `satelliteClass`:


        2. **Time Conversion and Solar Activity Data:**
        - Converts simulation time to Julian Date (JD) to extract solar activity data (F10 and F10bar)

        3. **Sun Position and Atmospheric Temperature:**
        - Determines the Sun's position in the ECI frame using `spiceTargPos`.
        - Calculates atmospheric temperature based on solar activity, Sun's 


        4. **Drag Force Calculation:**
        - Calls the `drag` function to compute the total acceleration vector 
            due to atmospheric drag. This function calculates 
            the force acting on each facet and sums them to obtain the 
            total force & torque vector. Density is calculated within the `drag` 
            function itself.



        **Arguments:**
            - self: Reference to the current object instance.
            - time: Simulation time (seconds).
            - satelliteClass: Reference to a class holding information about the spacecraft 
                            (e.g., mass, position, coefficients).

        **Returns:**
            tuple(numpy.ndarray): A tuple containing two NumPy ndarrays representing:
                - accVector (km/s^2): Total force vector acting on the spacecraft in the simulation frame due to atmospheric drag.
                - torqueVector (Nm): Total torque vector acting on the spacecraft in the simulation frame (currently set to zero).
        """


        #get position and velocity vectors
        r_eci = satelliteClass.eci_position()
        v_eci = satelliteClass.eci_velocity()

        #get satellite data:
        #STL data --each facet
        arr_area = satelliteClass.areas
        arr_normal = satelliteClass.normals
        arr_Tw = satelliteClass.Tws
        arr_alpha_coeff = satelliteClass.alpha_coeffs
        COMpos =  satelliteClass.positions - satelliteClass.centreOfMass      


        #properties:
        m = satelliteClass.mass
        quaternion = satelliteClass.quaternion


        #transform into epoch to julian date
        fractionalTime =  time / 86400
        day, month, year = self.epoch
        JD = julianDate(fractionalTime, day, month, year)
        

        #extract solar activity
        F10 = self.interpf10(JD)
        F10_avg = self.interpctr81(JD)

        #sun position:
        vec_sun = spiceTargPos(time, self.epoch, 'sun', self.centralBody)

        #atmospheric temperature
        T = atmTemperature(F10, F10_avg, vec_sun, r_eci, self.EarthFlattening)


        #compute acceleration and torque vector
        ##density computed inside the drag function
        accVector, torqueVector = drag(r_eci,
                        v_eci,
                        m,
                        self.omega_E,
                        self.R_E, 
                        arr_area, 
                        T, 
                        arr_Tw, 
                        quaternion, 
                        arr_normal, 
                        COMpos, 
                        arr_alpha_coeff)


        return accVector, torqueVector