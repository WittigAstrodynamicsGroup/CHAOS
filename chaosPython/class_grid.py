#python code to hold the grid class


'''
Grid class definition. The class holds the grid data and the methods to represents to actions of the grid. Subject to modification as 
the interaction with the wrapper / other classes gets implemented.
The grid attributes are the following

    args expected attributes:
    The first argument to be passed is a *args containing the arguments for the shapeFunction
        -numHeads: The number of thruster heads on the thruster / information on the thruster head number. This is to be combined with a generating function
        -columnLength: Number of pixel in a column or column equivalent array
        -rowLength: Number of pixel in a row of row equivalent array
    kwargs expected attributes
    -shapeFunction: Function generating the position of the thruster heads based on numHeads and the desired pattern / shape / distribution.
    -position: Holds the full position of each thruster heads. NOTE: not a direct input
    -burnTime: Holds the fuel mass of each thruster head, linked to their position. NOTE: could be burn time instead
    -thrustNominal: nominal thrust of the thruster heads
    -thrustVariation: thrust variation of the thruster heads

The grid methods are:

    -ignite(): Reduces the fuel / burn time of a thruster head (associated to a position). This should be coupled with the integration of the ODE at this position
    -gradient(): Allocate fuel based on the position of the thruster heads.
    -gradientPlot(): Plot the allocation of fuel for a visual understanding of the distribution. 


'''


import numpy as np

class Grid:

    def __init__(self,  **kwargs):

        numHeads = kwargs.get('numHeads', None)
        columnLength = kwargs.get('columnLength', None)
        rowLength = kwargs.get('rowLength', None)
        
        quadrantDivision = kwargs.get('quadrantDivision', None)
        totalburnTime = kwargs.get('totalburnTime', 365*24*3600) 
        self.fuelMass = kwargs.get('fuelMass', 0.1) #100g in kg
        pixelPosition = kwargs.get('pixelPosition', None)
        thrustNominal = kwargs.get('thrustNominal', 17.6e-6)
        pulseVariation = kwargs.get('pulseVariation', 0.1)
        pixelPulseRate = kwargs.get('pixelPulseRate', 20)
        QIT = kwargs.get('QIT', (totalburnTime/numHeads))
        self.thrusterOrientation = kwargs.get('thrusterOrientation', np.array([1, 0, 0]))
        self.thrustAxis = kwargs.get('thrustAxis', - self.thrusterOrientation)
        
        self.numHeads = numHeads
        self.columnLength = columnLength
        self.rowLength = rowLength
        self.quadrantDivision = quadrantDivision
        self.totalburnTime = totalburnTime
        self.pixelPosition = pixelPosition
        self.thrustNominal = thrustNominal
        self.pulseVariation = pulseVariation
        self.pixelPulseRate = pixelPulseRate
        self.QIT = QIT
        self.thrustVariation = (self.thrustNominal * self.pulseVariation) / (np.sqrt(self.QIT * self.pixelPulseRate))
        self.thrust_std = self.pulseVariation / (np.sqrt(self.QIT * self.pixelPulseRate))
        self.state = 0
        self.rng = np.random.default_rng()
        print('numhead interface required. Ensure the ifnormation is given, and not clashing with other data')
        '''
            #set the amount of pixels
            if self.numHeads == None and self.rowLength == None and self.columnLength == None:
                raise ValueError('Please input either the total number of heads of the amount of heads in row /  column')

            elif self.numHeads != None and self.columnLength != None and self.rowLength != None:
                if self.numHeads == self.rowLength * self.columnLength:
                    pass
                else: 
                    raise warnings.warn('The total number of pixel does not match, the grid is not square symmetrical')

            elif self.numHeads == None and self.columnLength == None and self.rowLength != None:
                self.numHeads = self.rowLength**2

            elif self.numHeads == None and self.rowLength == None and self.columnLength != None:
                self.numHeads =  self.columnLength**2

            elif self.numHeads == None and self.rowLength != None and self.columnLength != None:
                self.numHeads = self.rowLength * self.columnLength
            else: 
                self.numHeads = self.numHeads
        '''
        #set the position vecotrs

        #ensure not too much or too little data is given, and that the r_dict doesn't clash with shapefunction 
        # if self.shapeFunction == None and self.pixelPosition == None:
        #     raise ValueError('Need either an input function or manually input the position dict.')
        
        # elif self.shapeFunction != None and self.pixelPosition != None:
        #     raise ValueError('Cannot specify both a generating function and a manual input for the position. Set one of them to None')

        # elif self.pixelPosition ==  None and self.shapeFunction != None:
            
        #     self.pixelPosition = self.shapeFunction

        # elif self.pixelPosition != None and self.shapeFunction == None:
        #     self.pixelPosition = self.pixelPosition
        

        #assign the pixel layout:
        

        #set the fuel allocation and ignition list
        self.burnTime = {} #evolving dict where the remaining burntime of pixel is stored
        self.ignitionLog = {}
        self.firingTime = {} #default time for which pixel fires
        self.pixelDict = {}
        if self.numHeads == len(self.pixelPosition.keys()):
            #create dict, each pixel poisition is allocated a fuel mass (uniform)
            for i in self.pixelPosition.keys():
                self.burnTime[i] = (self.totalburnTime / self.numHeads)
                self.firingTime[i] = (self.QIT)
                self.ignitionLog[i] = 0
                self.pixelDict[i] = []

        
        
        else:
            raise ValueError('The number of heads does not match the dict key length. i.e. the number of heads  {}\t is different than the amount of head positions.{} '.format(self.numHeads, len(self.pixelPosition.keys())))
        

    #define methods here


    def updateFiringTime(self, time):
        '''
            Updates the list of firing time, to mimic non-uniform fuel distribution 
        '''
        for i in self.pixelPosition.keys():
            self.firingTime[i] = (time)


    def generateTorque(self, chosenPixel):
        '''
            generates torque for a given pixel position. Uses he object's inbuilt data.
        '''
        
        #generate a random thrust in the thrust axis direction
        thrust =  self.rng.normal(self.thrustNominal, self.thrustVariation) * (self.thrustAxis)

        #torque vector, cross product of r and thrust vectors
        torque = np.cross(self.pixelPosition[chosenPixel], thrust) 
        return(thrust, torque)
    


    def quadrant(self, r_dict):
        ##divides the grid of pixel into quadrants, based on the user-defined function.
        dicts = self.quadrantDivision(r_dict)

        return dicts
        

if __name__ =='__main__':
    from class_shapeFunction import shapeFunction
    import random
    uniform = shapeFunction() 
    
    grid1 = Grid(10, shapeFunction = uniform.uniformGrid, numHeads=100)
    I, II, III, IV =  grid1.quadrant(grid1.pixelPosition)
    # print(I)
    xI = []
    yI = []
    yII = []
    xII = []
    xIII = []
    yIII = []
    xIV = []
    yIV = []
    
    import matplotlib.pyplot as plt
    for i in I.keys():
        xI.append(I[i][1])
        yI.append(I[i][2])

    for i in II.keys():
        xII.append(II[i][1])
        yII.append(II[i][2])

    for i in III.keys():
        xIII.append(III[i][1])
        yIII.append(III[i][2])

    for i in IV.keys():
        xIV.append(IV[i][1])
        yIV.append(IV[i][2])

    plt.scatter(xI, yI, c='red')
    # plt.scatter(xII, yII, c='blue')
    # plt.scatter(xIII, yIII, c='yellow')
    # plt.scatter(xIV, yIV, c= 'green')
    plt.ylim((-0.04, 0.04))
    plt.xlim((-0.04, 0.04))
    plt.grid(True, linestyle='dashed')
    plt.savefig('localPlot/quadrants.png')
    # plt.show()
    mags = [np.linalg.norm(i[1:]) for i in I.values()]
    print(mags)
    print('val', (((17.6e-6 * max(mags))/0.00182)*100)*180/np.pi)
    print('computed std', np.sqrt(((np.max(mags)-np.min(mags))**2)/12), 'computed mean', ((np.max(mags)+np.min(mags)))/2)
    
    
    
    iterations = 500000
    positions = []
    for i in range(iterations):

        pos = random.choice(mags)
        positions.append(pos)
    
    print('test mean', np.mean(positions))
    print('test std', np.std(positions))

    print(np.mean(mags))
    print(np.std(mags))
    # print('position:', grid1.pixelPosition["'2 2'"])
    # print('thrust, torque:', grid1.generateTorque("'2 2'"), grid1.thrustNominal)
    # print(grid1.burnTime)
    # print(grid1.totalburnTime)

