#python script to hold the sensor class

"""
This python scripts holds the code to the sensor class. Few definition will be held here as a trial for 
"true oop", so that inputs can also be user-defined functions
"""

from .sensors import faradayCupAngle, MEMSgyro
import numpy as np
from numpy.random import default_rng
class Sensor:

    def __init__(self, **kwargs):
        #define initial attributes
        #Faraday Cup
        self.FCepsilon = kwargs.get('currentError', None)
        self.FCaperture = kwargs.get('apertureArea', None)
        self.FCheight = kwargs.get('FCheight', 10*1e-3)
        self.FCradius = kwargs.get('FCradius', 15*1e-3)
        self.FC_SN = kwargs.get('SN', 1)


        #MEMS gyro
        self.biais_n = kwargs.get('MEMSbiais', np.array([0, 0, 0]))
        self.sigma_n = kwargs.get('sigmaNoise', (0.15*np.pi/180)/60)
        self.sigma_b = kwargs.get('sigmaBiais', (0.3*np.pi/180)/3600)
        
        
        #define random number generator
        self.seed = kwargs.get('seed', None) #for testing purposes
        self.seed2 = kwargs.get('seed2', None)
        self.rng = default_rng(seed=self.seed)
        self.rng2 = default_rng(seed=self.seed2)
    #define methods here



    def coarseAttitude(self, alpha, v, rho):
        v = v*1e3
        angle = faradayCupAngle(alpha, self.FCepsilon, v, self.FCaperture, rho, self.rng)

        return angle

    def MEMSgyro_AngularVelocity(self, omega, t, t_prev, biais_n):
            delta_t = t - t_prev
            biais_new, ang_vel = MEMSgyro(omega, delta_t, biais_n, self.sigma_n, self.sigma_b, self.rng, self.rng2)

            return biais_new, ang_vel


