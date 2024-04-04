#python code to hold dataset class


import numpy as np


'''
Storage class that handles storing and parsing data and writing it to a textfile.
NOTE: The data is stored only in its instance. The data stored here is wiped after each iteration, unless written into a file.
For consistency, everything to be written should be stored within the object. 
The methods available are the following:

        -parse:         Handle the basic parsing of data from simulator. This include splitting arrays into components and else
        -write_data:    Write 5 different arrays to a textfile. Each column is a data array. 
        -write_value:   Write 5 different values(1 number) to a textfile. Each column is a value

'''
class Dataset:
    #No input argument when instantiated
    def __init__(self):
        #create a bunch of list, dicts, etc for storage
        self.omega_x = []
        self.omega_y = []
        self.omega_z = []
        self.omega_norm = []
        
        self.quaternion_x = []
        self.quaternion_y = []
        self.quaternion_z = []
        self.quaternion_r = []
        # self.quaternion_norm = []
        self.UnitQuaternion_x = []
        self.UnitQuaternion_y = []
        self.UnitQuaternion_z = []
        self.UnitQuaternion_r = []

        self.time_dump = []

        self.thrust = []
        self.thrust_x = []
        self.thrust_y = []
        self.thrust_z = []

        self.torque = []
        self.torque_x = []
        self.torque_y = []
        self.torque_z = []
        self.torque_norm = []

        self.sequence = []

        self.InitialFuel = 0
        self.remainingFuel = []
        self.totalFuel = 0
        self.stdFuel = 0

        self.totalImpulse = []
        self.pixelImpulse = []

        self.ignitionLog = []
        self.maxIgnition = 0
        self.meanIgnition = 0
        self.stdIgnition = 0


        self.maxVelocity = []

        self.firingTime = {}
        self.meanFiringTime = {}
        self.minFiringTime = {}
        self.absoluteMinFiringTime = 0
        
        self.breaker = 0
        self.t_start = 0
        self.y0 = None
    
    # define methods here




    def parse(self):
        '''
        Method to automatically parse all the data when called
        '''
        #set time in days
        self.time_dump = self.time_dump #have the time in seconds

        #compute norm of omega
        self.omega_norm = np.sqrt(self.omega_x**2 + self.omega_y**2 + self.omega_z**2)
        self.maxVelocity = max(self.omega_norm)

        #split torque components:
        for i in range(len(self.torque)):
            self.torque_x = np.append(self.torque_x, self.torque[i][0])
            self.torque_y = np.append(self.torque_y, self.torque[i][1])
            self.torque_z = np.append(self.torque_z, self.torque[i][2])
            self.torque_norm = np.sqrt(self.torque_x**2 + self.torque_y**2 + self.torque_z**2)
        
        
        #split thrust components:
        for i in range(len(self.torque)):
            self.thrust_x = np.append(self.thrust_x, self.thrust[i][0])
            self.thrust_y = np.append(self.thrust_y, self.thrust[i][1])
            self.thrust_z = np.append(self.thrust_z, self.thrust[i][2])

        #calculate total impulse:
        self.totalImpulse = sum(self.pixelImpulse)

        #calculate average firing time:
        abis = []

        for i in self.firingTime.keys():
            abis = np.append(abis, self.firingTime[i])
            if len(self.firingTime[i]) != 0:
                self.meanFiringTime[i] = np.mean(self.firingTime[i])
            # if len(self.minFiringTime[i]) != None:
                self.minFiringTime[i] = np.min(self.firingTime[i])
        
        if len(abis) != 0:
            self.absoluteMinFiringTime = min(abis)

        #calculate the average ignition per pixel:
        self.meanIgnition = np.mean(list(self.ignitionLog.values()))
        self.maxIgnition = max(list(self.ignitionLog.values()))

        #calculate the total remaining fuel:
        self.totalFuel = 100 * (self.InitialFuel - sum(self.remainingFuel.values())) / (self.InitialFuel)
        
        #unitise quaternions before writing. Divide each element of a quaternion by the norm of said quaternion.
        #compute norm
        self.quaternion_norm = np.sqrt(np.square(self.quaternion_x) + np.square(self.quaternion_y)  + np.square(self.quaternion_z) + np.square(self.quaternion_r))
        # quaternion_norm = 1

        #divide each element by the norm
        self.UnitQuaternion_x = self.quaternion_x 
        self.UnitQuaternion_y = self.quaternion_y 
        self.UnitQuaternion_z = self.quaternion_z 
        self.UnitQuaternion_r = self.quaternion_r 



    def write_data(self, filename, data1, data2, data3, data4, data5):
        '''
        method to write  5 data arrays into a file 
        '''
        
        f = open(filename + '.txt', 'a')
        # f = filename
        for i in range(len(data1)):
            f.write("{}\t{}\t{}\t{}\t{}\n".format(data1[i], data2[i], data3[i], data4[i], data5[i]))
        f.close()

    
    def write_value( self, filename, *args):
        # f = open(filename + '.txt', 'a')
        f = filename
        f.write(''.join('%s\t' % item for item in args))
        f.write('\n')
        # f.close()

def write_value( filename, *args):
    # f = open(filename + '.txt', 'a')
    f = filename
    f.write(''.join('%s\t' % item for item in args))
    f.write('\n')
    # f.close()


# def write_value( filename, *args):
#     # f = open(filename + '.txt', 'a')
#     f = filename
#     f.write(''.join('%s\t' % item for item in args))
#     f.write('\n')
#     # f.close()
'''
test code
'''

if __name__ == '__main__':
    t=[]
    x=[]
    y=[]
    z=[]
    norm=[]
    # timetrack = open('timetrack0.txt', 'r')
    # for line in timetrack.readlines():
    #     item = line.split('\t')
    #     t.append(float(item[0]))
    #     x.append(float(item[1]))
    #     y.append(float(item[2]))
    #     z.append(float(item[3]))
    #     norm.append(float(item[4]))
    # print(type(z[0]))
    abis = []
    a = {1:[10], 2:[5, 1], 3:[8, 0, 4]}
    for i in a.keys():
        abis = np.append(abis, a[i])
        a[i] = np.mean(a[i])
    print(a)
    print(min(abis))



    import os
    print(os.getcwd())







    
    

