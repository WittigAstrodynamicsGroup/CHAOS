#Python script to hold the Panel class

"""
Panel Class definition.
Holds the mechanical properties of the panel (Cd, Cl, Cr, etc...) and calls functions to compute
the force and torque acting on itself. 

NOTE: Not all forces and torques can be computed from panel level!! ex: EGM, LS, GG....
"""

import numpy as np
# 
class Panel:

        def __init__(self, **kwargs):
            
                # self.env = Environment
                #Surface properties
                # self.Cd = kwargs.get('Cd', 2.2)
                self.Cr = kwargs.get('Cr', 0.8)

                self.area = kwargs.get('area', 0.01) #m^2
                self.surfNormal = kwargs.get('normal', None)
                self.surfPos = kwargs.get('pos', None)

                self.length = kwargs.get('length', 0.1)
                self.width = kwargs.get('width', 0.1)

                self.Tw = kwargs.get('Tw', 300) # K
                self.alpha = kwargs.get('alpha', 0.9)
  

    ##call methods here

        # def Cd(self, angle, Ta, Vorb):

        #         Cd, Cl = MoeSentman(angle, self.area, 1, Ta, self.Tw, Vorb, self.alpha) 

        #         return Cd, Cl