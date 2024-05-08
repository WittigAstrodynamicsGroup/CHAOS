
"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_sensor.py


This Python module defines the `Sensor` class, which simulates the behavior 
of two satellite sensors: a Faraday cup and a MEMS gyroscope. The class 
provides methods to calculate coarse attitude measurements from the Faraday 
cup and noisy angular velocity measurements from the MEMS gyro.
"""

from .sensors import faradayCupAngle, MEMSgyro
import numpy as np
from numpy.random import default_rng
class Sensor:

    def __init__(self, **kwargs):
        #define initial attributes
        #Faraday Cup
        self.FCepsilon = kwargs.get('currentError', None)                           # Current measurement error
        self.FCaperture = kwargs.get('apertureArea', None)                          # Aperture of Faraday cup                         
        self.FCheight = kwargs.get('FCheight', 10*1e-3)                             # Faraday cup height
        self.FCradius = kwargs.get('FCradius', 15*1e-3)                             # Faraday cup radius
        self.FC_SN = kwargs.get('SN', 1)                                            # Signal-to-noise ratio threshold 


        #MEMS gyro
        self.biais_n = kwargs.get('MEMSbiais', np.array([0, 0, 0]))                 # Biais instability of MEMS gyro
        self.sigma_n = kwargs.get('sigmaNoise', (0.15*np.pi/180)/60)                # MEMS gyro noise
        self.sigma_b = kwargs.get('sigmaBiais', (0.3*np.pi/180)/3600)               # MEMS gyros biais
        
        
        #define random number generator
        self.seed = kwargs.get('seed', None)                                        #Fix the seeds to an int to be able to repeat the ramdon sequence
        self.seed2 = kwargs.get('seed2', None)
        self.rng = default_rng(seed=self.seed)
        self.rng2 = default_rng(seed=self.seed2)
    #define methods here


        #Simulates the measured angle based on the error of the Faraday cup
    def coarseAttitude(self, alpha, v, rho):
        v = v*1e3
        angle = faradayCupAngle(alpha, self.FCepsilon, v, self.FCaperture, rho, self.rng)

        return angle



    #Simulates the measurement of a MEMS gyroscope, including bias drift and random noise.
    def MEMSgyro_AngularVelocity(self, omega, t, t_prev, biais_n):
            delta_t = t - t_prev
            biais_new, ang_vel = MEMSgyro(omega, delta_t, biais_n, self.sigma_n, self.sigma_b, self.rng, self.rng2)

            return biais_new, ang_vel


