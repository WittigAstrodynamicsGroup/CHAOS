#python code to hold the grid class


'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_grid.py


This Python module defines the `Grid` class, which represents a grid of 
thruster heads on a CubeSat surface. The grid object holds data about the 
thruster configuration and provides methods to manage thruster operations.
'''


import numpy as np

class Grid:

    def __init__(self,  **kwargs):

        self.numHeads = kwargs.get('numHeads', None)                                 # Number of pixels on the thruster grid
        self.quadrantDivision = kwargs.get('quadrantDivision', None)                 # Function to divide pixels into quadrants
        self.totalburnTime = kwargs.get('totalburnTime', 365*24*3600)                # Total firing time of the thruster
        self.fuelMass = kwargs.get('fuelMass', 0.1) #100g in kg                      # total fuel mass on the grid
        self.pixelPosition = kwargs.get('pixelPosition', None)                       # Dict containing the position of the pixels, 
                                                                                     # in the format {key: positionVector}
        self.thrustNominal = kwargs.get('thrustNominal', 17.6e-6)                    # Nominal thrust to be used
        self.pulseVariation = kwargs.get('pulseVariation', 0.1)                      # Variation of the thrust at each pulse [0, 1]
        self.pixelPulseRate = kwargs.get('pixelPulseRate', 20)                       # number of Pulses per seconds
        self.QIT = kwargs.get('QIT', 100)                                            # Quadrant Ignition Time: Max. time a thruster
                                                                                     # fires before the control systems is called again
        self.thrusterOrientation = kwargs.get('thrusterOrientation',                 # Body-fixed direction in which the thruster fires
                                               np.array([1, 0, 0]))
        self.thrustAxis = kwargs.get('thrustAxis', - self.thrusterOrientation)       # The direction in which the acceleration 
                                                                                     # due to the thrust will act 

                                                    
        #Use Central Limit Theorem to compute variation of the average thrust:
        self.thrustVariation = (self.thrustNominal * self.pulseVariation) / (np.sqrt(self.QIT * self.pixelPulseRate))
        
        #Initialise random number generator for thrust randomisation
        self.rng = np.random.default_rng()

    

        ##Create dictionaries to store remaining fuel information for each pixel
        #set the fuel allocation and ignition list
        self.burnTime = {}                      #evolving dict where the remaining burntime of pixel is stored
        self.firingTime = {}                    #default time for which pixel fires
        
        
        if self.numHeads == len(self.pixelPosition.keys()):
            #for each pixel
            for i in self.pixelPosition.keys():
                self.burnTime[i] = (self.totalburnTime / self.numHeads)     #each pixel is allocated the same burntime 
                self.firingTime[i] = (self.QIT)                             #A pixel fires for QIT  seconds before CS is called


    
        else:
            raise ValueError('The number of heads does not match the dict key length. i.e. the number of heads  {}\t is different than the amount of head positions.{} '.format(self.numHeads, len(self.pixelPosition.keys())))
        

    #define methods here


    def updateFiringTime(self, time):
        '''
        Updates the list of firing time
        '''
        for i in self.pixelPosition.keys():
            self.firingTime[i] = (time)


    def generateTorque(self, chosenPixel):
        '''
        generates torque and thrust for a given pixel position. 
        Uses the object's inbuilt data.
        '''
        
        #generate a random thrust in the thrust axis direction
        #according to a gaussian centered of ThrustNominal 
        thrust =  self.rng.normal(self.thrustNominal, self.thrustVariation) * (self.thrustAxis)

        #torque vector, cross product of r and thrust vectors
        torque = np.cross(self.pixelPosition[chosenPixel], thrust) 
        return(thrust, torque)
    


    def quadrant(self, r_dict):
        """
        divides the grid of pixel into quadrants, based on the user-defined function.
        """

        dicts = self.quadrantDivision(r_dict)

        return dicts
        

