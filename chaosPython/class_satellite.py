"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_satellite.py



This module defines the `Satellite` class, which represents a spacecraft 
in a simulation environment. 


**State initialization:**

* The class allows for various initialization options using keyword arguments. 
  Default values are provided for a standard 1U CubeSat.

* Users can specify initial orbital elements, state vectors (position and 
  velocity), or both. The class performs the necessary conversions to ensure 
  consistent internal representation.

**Geometry:**

* The class supports An STL file defining the spacecraft's geometry.

* When using an STL file, the class calculates relevant properties 
  (areas, normals, facet positions) based on the mesh data, considering 
  the provided scale factor to convert units from the STL file to meters.

**Methods:**

* The class defines several methods, including:


    * `eci_position(self)`: Computes the spacecraft's ECI position vector 
      based on its modified equinoctial elements using an external 
      `equi_to_r_cart` function.

    * `eci_velocity(self)`: Computes the spacecraft's ECI velocity vector 
      based on its modified equinoctial elements using an external 
      `equi_to_v_cart` function.
    """




import numpy as np
from stl import mesh
from .transformations import kep_to_eq, rv2kpl, kpl2rv, equi_to_r_cart, equi_to_v_cart
from .support_fun import gridsBurntime
import warnings


class Satellite:
    #keyword input argument--default values for standard 1U CubeSat.
    def __init__(self,  **kwargs):
        
        


        #####################
        # SATELLTIE STATE
        #####################
        #attitude
        self.quaternion = kwargs.get('quaternion', None)            #initial quaternion
        self.angularVel = kwargs.get('angularVel', None)            #Initial angular velocity [rad]

        #orbital cartesian
        self.velocity = kwargs.get('velocity', None)                #Initial orbital velocity [km/s]
        self.position = kwargs.get('position', None)                #initial orbital position [km]

        #orbital element
        #defaults to ISS orbit
        self.keplerianElements = kwargs.get('keplerianElements', np.array([6781.36,         #SMA
                                                                            0.0021168,      #ECC
                                                                            51.6*np.pi/180, #INC
                                                                            90*np.pi/180,   #RAAN
                                                                            0,              #WP
                                                                            0]))            #TA
        #Modified equinoctial elements
        self.modEquinoctialElements = kwargs.get('equinoctial', None)




        #####################
        # SATELLITE DATA
        #####################
        #defaults to a 1U CubeSat
        self.mass = kwargs.get('mass', 1.33)                #mass [kg]
        self.Cd = kwargs.get('Cd', 2.2)                     #drag coefficient
        self.Cr = kwargs.get('Cr', 0.8)                     #reflection coefficient, for SRP
        self.area = kwargs.get('area', 0.015)               # average area, for low-fidelity models
        self.Tw = kwargs.get('Tw', 300) #k                  #Satellite wall temperature
        self.alpha_coeff = kwargs.get('alpha_coeff', 0.9)   #Alpha coefficient (Sentman Model)




        #####################
        # GEOMETRY
        #####################
        # self.panelList = kwargs.get('panelList', [])
        self.STLfile = kwargs.get('STLfile', None)          #Stores the location of the STL file
        self.scale_factor = kwargs.get('scale_factor', 1e6) #how many units of STL 
                                                            #file in 1 m2?
                                                            #----Typically, STL files uses millimeters^2 as unit


        self.inertiaMatrix = kwargs.get('inertiaMatrix', np.diag([0.00216667]*3))   #Inertia matrix of the satellite
        self.InverseinertiaMatrix = np.linalg.inv(self.inertiaMatrix)               #Inverse of the inerta matrix
        self.centreOfMass = kwargs.get('centreOfMass', np.array([0, 0, 0]))         #Position of the satellite CoM

        

        #####################
        # CONSTANTS
        #####################
        self.mu = 398600.4418  #environment.mu                  # Earth gravity parameter #398600.4418   #km^5- from vallado
        self.mu_E = self.mu * 1e9                               # Earth gravity parameter- SI




        #####################
        # INITIALISATION
        #####################

        #Based on which variable is given, transform to modifiedEquinoctial
        if self.position is not None and self.velocity is not None:
            #get keplerian elements from r,v in km & km/s
            self.keplerianElements = rv2kpl(self.position, self.velocity)

        else:
            # get eci r, v in km & km/s
            self.position, self.velocity = kpl2rv(self.keplerianElements[0],        #SMA
                                                    self.keplerianElements[1],      #ECC
                                                    self.keplerianElements[2],      #INC 
                                                    self.keplerianElements[3],      #RAAN 
                                                    self.keplerianElements[4],      #AOP 
                                                    self.keplerianElements[5],      #THETA
                                                    self.mu)                        #MU--in km^3/s^2
        
        #transform into modEq
        a, e, inc, raan, wp, theta = self.keplerianElements

        #assign to attributes
        self.modEquinoctialElements = kep_to_eq(a, e, inc, raan, wp, theta)



        ###########################
        # INITIALISATION OF GEOMETRY
        ###########################

        if self.STLfile is not None:

            #read stl file
            HFstl = mesh.Mesh.from_file(self.STLfile)           #Load STL file


            #Convert the STL units (typically mm^2) to m^2
            self.areas = np.array([HFstl.areas[i][0] / self.scale_factor for i in range(len(HFstl.areas))])


            #get normals to the facets
            self.normals = HFstl.get_unit_normals()

            #get position of CoM of each facet
            self.positions = [np.array([np.sum(HFstl.x[i])/3, np.sum(HFstl.y[i])/3, np.sum(HFstl.z[i])/3])/np.sqrt(self.scale_factor) for i in range(len(HFstl.x))]
            
            #For each facet, provide Temperature, reflection coefficient, and alpha coefficient
            self.Tws = [self.Tw for i in range(len(HFstl.areas))]                       #Facet Wall temperature
            self.Crs = [self.Cr for i in range(len(HFstl.areas))]                       #Facet reflection coefficient
            self.alpha_coeffs = [self.alpha_coeff for i in range(len(HFstl.areas))]     #Facet alpha coefficient


            #If you have many facets, pop warning
            if len(HFstl.areas) > 30:
                warnings.warn("WARNING: {} facets detected. Is this normal? Slow run time can be expected.".format(len(HFstl.areas)))


            
    
    
    #####################
    # DEFINE METHODS HERE
    #####################


    


    def eci_position(self):

        """
        Calculates the spacecraft's ECI position vector.

        This function computes the ECI (Earth-Centered Inertial) position vector of the 
        spacecraft based on its modified equinoctial elements. It uses the 
        `equi_to_r_cart` function (written in`transformations.py`)
        to perform the conversion. The result is stored in both a class 
        attribute (`self.position`) and returned.

        Args:
            self: Reference to the spacecraft object.

        Returns:
            numpy.array: The ECI position vector of the spacecraft (km).
        """


        #compute eci position from modified equinoctial elements
       
        p, f, g, h, k, L = self.modEquinoctialElements          ##Get the equinoctial elements
        
        r_ECI = equi_to_r_cart(p, f, g, h, k, L)                ##Function in "transformations.py"
        
        self.position = r_ECI                                   ##Update the class attribute

        return r_ECI






    def eci_velocity(self):
        """
        Calculates the spacecraft's ECI velocity vector.

        This function computes the ECI (Earth-Centered Inertial) velocity vector of the 
        spacecraft based on its modified equinoctial elements. The result is stored 
        in both a class attribute (`self.velocity`) and returned.

        Args:
            self: Reference to the spacecraft object.

        Returns:
            numpy.array: The ECI velocity vector of the spacecraft (km/s).
        """

        #compute eci velocity from modified equinoctial elements

        v_ECI = equi_to_v_cart(self.modEquinoctialElements)     #function in `transformations.py`

        self.velocity = v_ECI                                   ##sUpdate the class attribute

        return v_ECI






   

   