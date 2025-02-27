"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


set_up.py


This script sets up the configuration and parameters for running a (optional: parallel) simulation using the CHAOS module.
It initializes the environment, defines the satellite parameters, control system settings, and other relevant parameters.
The simulation is run using the chaosPython module, which handles the dynamics of the satellite and environmental interactions.

The script also includes functionality for checkpointing to resume simulations in case of interruption and visualization of simulation results.

Dependencies:
- spiceypy (https://spiceypy.readthedocs.io/en/main/installation.html#how-to-install-from-source-for-bleeding-edge-updates)
- NAIF SPICE Kernels (https://naif.jpl.nasa.gov/pub/naif/generic_kernels/)
- scipy
- numpy
- mpi4py (if parallel simulations required)
- chaosPython (custom module)

Usage:
    python3 set_up.py

"""

import os 
import numpy as np
from mpi4py import MPI

import chaosPython as cp

"""SET-UP FOR SIMULATIONS"""


######################################
#initial keplerian orbital parameters
######################################
sma = 6378 + 350                    #semi-major axis
ecc = 0                             #eccentricity
inc = 98.7  *np.pi/180              #inclination
raan = 0    *np.pi/180              #Right Ascension
wp = 0      *np.pi/180              #argument of perigee
trueAnomaly = 0 *np.pi/180          #True Anomaly




######################################
#Simulations Identifiers
#   Set a unique identifier to the simulation
######################################

#For parallel simulations
# comm = MPI.COMM_WORLD               #Initialise the MPI communicator
# rank = comm.Get_rank()              #get the rank of this current process

#If you want only one simulation
rank  = 'test'                         #The number here will be the simulation ID   

#If you need to run a few specific simulation IDs 
# to_run = [1, 2, 3]                #This example runs simulation 1, 2, 3 in parallel
# rankCPU = comm.Get_rank()
# rank = to_run[rankCPU]



######################################
#Define the force models to be used
######################################

#Call the perturbation classes
#Each of them has a "computePerturb(t, satelliteClass)" method, 
#and returns a tuple (accelerationVector, TorqueVector) [km/s^2, Nm].
#To create your own model, create a class derived from perturbationInterface (cp.perturbationInterface)

#EGM
egm = cp.EarthGravityModel()

#Moon /Sun gravity
MoonGravity = cp.nBody(gravitationalParameter=4.902800305555e3,                 #GMAT value: 4902.8005821478
                        SPICEname='moon') 
SunGravity = cp.nBody(gravitationalParameter=1.327122e11,                       #GMAT value: 132712440017.99
                       SPICEname='sun')      


#SRP -- low fidelity or high fidelity to choose from:
LFsrp = cp.LowFidelitySRP()
HFsrp = cp.HighFidelitySRP()


#Drag-- Low fidelity or high fidelity to choose from:
LFdrag = cp.LowFidelityDrag()
HFdrag = cp.HighFidelityDrag()






######################################
#Call on a manager class to handle
#force models
######################################


#All the perturbation classes are handled by the perturbation Manager class:
manager = cp.perturbationManager()

#give the perturbation desired for the simulation --orbital and attitude perturbation models!!
## --- empty list for keplerian motion
manager.forceModelsList = [egm, MoonGravity, SunGravity, HFdrag, HFsrp]

#set simulation details:
manager.epoch = np.array([1, 7, 2014])                              #Simulation start time: i.e. 1st of July 2014

#set central body--
#These are passed to the 
#perturbation classes
manager.centralBody = 'earth'                                       #Central body of simulation 
manager.mu = 398600.4418                                            #associated gravitational parameter    
manager.bodyRadius = 6378.1363                                      #Primary body radius


#IMPORTANT:
#This updates the epoch and central body in 
#all perturbation classes appended in forceModelsList
manager.initialSetUp()                                              #call this at the end of the manager set-up!                                        

######################################
#Initialise the environment model
##in this case it's ionic density only
######################################

LEO = cp.Environment()


######################################
#initialise grid of Cube de alps 
# -- Assume copper is used as fuel source
######################################

fuel_burn_time=365*86400                                    #Total fuel burn time is seconds

############
# Thruster 1
# Set on face +x 
############


face1PixelPosition = cp.AA_2FC_grid()                        #Dictionary of pixel positions with 
                                                            #format {PixelID: pixelPositionVector}


#instantiate the thruster grid class
cube_de_alps = cp.Grid(numHeads = 44,                       #number of pixels
                       pixelPosition = face1PixelPosition,    #dictionary containing the positions of the pixels
                       pulseVariation=0.1,                  #standard deviation of each pulse
                       QIT = 100,                           #Quadrant ignition time
                       totalburnTime=fuel_burn_time/2,      #Fuel [in seconds] for this grid
                       thrusterOrientation = np.array([1, 0, 0]))        #Thrust direction for this grid.

#set the fuel mass
cube_de_alps.fuelMass = 0.0801/2                            #fuel mass in kg ---80.1g
cube_de_alps.state = 0                                       #on/off switch. Inconsequential, as CHAOS 
                                                                #will evaluate when this should be on or off
cube_de_alps.quadrantDivision = cp.basicQuadrant             #Function that divides the pixels into the desired quadrants
cube_de_alps.name = 'FaceA'                                  #Give a name to the grid




############
# Thruster 2
# Set on face -x 
############


oppositeFacePixelPosition = cp.oppositeFace(cp.AA_2FC_grid())   #Dictionary of pixel positions with 
                                                                #format {PixelID: pixelPositionVector}

#instantiate the thruster grid class
cube_de_alps2 = cp.Grid(numHeads = 44,                       #number of pixels
                       pixelPosition = oppositeFacePixelPosition,    #dictionary containing the positions of the pixels
                       pulseVariation=0.1,                  #standard deviation of each pulse
                       QIT = 100,                           #Quadrant ignition time
                       totalburnTime=fuel_burn_time/2,      #Fuel [in seconds] for this grid
                       thrusterOrientation = np.array([-1, 0, 0]))        #Thrust direction for this grid.

#set the fuel mass
cube_de_alps2.fuelMass = 0.0801/2                            #fuel mass in kg ---80.1g
cube_de_alps2.state = 0                                       #on/off switch. Inconsequential, as CHAOS 
                                                                #will evaluate when this should be on or off
cube_de_alps2.quadrantDivision = cp.basicQuadrant             #Function that divides the pixels into the desired quadrants
cube_de_alps2.name = 'FaceB'                                  #Give a name to the grid

###########################
#initialise Satellite 
#which holds Cube-de-ALPS 
###########################

thruster_list = [cube_de_alps, cube_de_alps2]                 #Create a list of all the thruster objects






###############################
# initialise the satellite class
###############################

inertia = np.array([[0.00182,0,0], [0,0.00185,0], [0,0,0.00220]])           #1U CubeSat--SLUCUBE inertia matrix

#initial position
kep = np.array([sma, ecc, inc, raan, wp, trueAnomaly])                      #orbital parameters 


cubesat  = cp.Satellite(inertiaMatrix=inertia,                  #Inertia matrix of the cubesat
                        keplerianElements=kep,                  #Passes the keplerian element of the CubeSat. this 
                                                                #is taken into account only as the initial position.
                                                                #If a checkpointing file exist, this value in inconsequential

                        mass=1.2,                               #Total wet mass of the satellite
                        STLfile= cp.moduleDirName + '/Data/STLs/1UCubesat.STL')     #Feed an STL 
                                                                                    #file holding the geometry of the satellite
                                                                                    #Needed for full the high-fidelity perturbation modelling.
                                                                                    

###If you want the dry mass, remove the fuel!!!
# cubesat.mass = cubesat.mass - cube_de_alps.fuelMass - cube_de_alps2.fuelMass


##################################################
#Finalise setting initial condition for satellite
#Attitude 
##################################################

cubesat.quaternion = np.array([0.5, 0.5, 0.5, 0.5])     ##Quaternion giving the initial pointing --norm should be 1

#angular velocity --radians
cubesat.angularVel = np.array([1, 0.1, 0.2])*np.pi/180     ##Initial angular velocity vector, in radians




###############################################
###initialise Sensor class
# Holds Faraday cup information 
# HOlds gyroscope information
###############################################


#Faraday cup
sensor = cp.Sensor( apertureArea=7.08e-4,           #Sets the aperture area of the Faraday cup
                    currentError=1e-9,              #Sets the expected background current error reading
                    SN=5 )                          # Sets the signal to noise ratio above which the thruster is triggered


#Gyroscope: STIM277H---An aerospace-grade 3-axis angular velocity sensor
##Random Walk parameters
sensor.sigma_b = 0*(0.3/3600)*np.pi/180     #Bias instability: Secular growth of the standard deviation over time   
sensor.sigma_n = 0*(0.15/60)*np.pi/180      #Angular Random Walk: noise at each step





###############################################
###Initialise control system class
# Holds information about the terminal events
# Controls the firing decision
###############################################



#terminal event functions 
cp.pixel_fuel.terminal = True                       ##Stops the integration if the selected pixel runs out of fuel
cp.r_stop.terminal = True                           ## Stops the integration if the lower altitude 
                                                    # limit of 150 km has been reached

cp.assessGridState.terminal = True                  ##Stops the integration if the first grid crosses the cone boundary
cp.assessGridState2.terminal = True                 ##Stops the integration if the second grid crosses the cone boundary

#list of event functions
events = [cp.pixel_fuel, cp.r_stop, cp.sensorMeasurement, cp.assessGridState, cp.assessGridState2]

##Initialise the control system, and specify the event functions
CS = cp.control_system(cut_off_function=events)

CS.CSfunction = CS.quadrantControl                  #Specify the control system function

CS.switch=1                                         #on/off switch -- If off (=0), 
                                                    # the thruster will not fire in the simulation.
                                                    #The integration will be done in one leg, no interruptions.



###############################################
###Initialise Dataset class
##Holds functions to write the data to file
##Holds the start time data
##Holds Restart data if necessary
###############################################


store = cp.Dataset()                                #Dataset class
store.t_start = 0                                   #Time at the start of the simulation
outputPath = 'output'                               #Name of the output path



###############################################
###Checkpointing code
#Checks if a restart file exists
#If a restart file exists, reads data
##              and restarts from there
###############################################

#Checks if a restart file exists
if os.path.exists(outputPath + '/run{}/checkpoint.npz'.format(rank)) == True:
    npzfile = np.load(outputPath + '/run{}/checkpoint.npz'.format(rank), allow_pickle=True)         #Load the restart file
    
    t_start = npzfile['t']                                          ##Read the last saved time stamp              
    y0 = npzfile['y0']                                              #Read the state at the last saved time stamp
    fuelDict = npzfile['fuel'] #list of dicts                       #Read the grid fuel at the last saved ....
    previous_sensor_biais = npzfile['sensor_biais']                 #Read sensor data

    store.t_start = t_start                                         #Overwrites the simulation start time
    store.y0 = y0                                                   #Overwrite the simulatuion initial condition

    for i in range(len(thruster_list)):
        print(i)
        thruster_list[i].burnTime = fuelDict[i]                     ##Overwrites the fuel for each pixel, on each grid
    sensor.biais_n = previous_sensor_biais                          ##The error in the sensor reading is updated





###############################################
###Final set-up
###############################################
    

#max simulated time -- This is a hard limit, no exception will be made
stopper = 86400 * 1  # 1 day in seconds



# list of SPICE kernels to load
    ##Make sure you have downloaded these files from 
    ###     https://naif.jpl.nasa.gov/pub/naif/generic_kernels/
    ###     and store them in 
    #          "./chaosPython/Data/kernels/{}"

kernels_to_load = ['/lsk/naif0012.tls',
                       '/spk/de440.bsp',
                       '/fk/earth_itrf93.tf',
                       '/pck/pck00010.tpc',
                       '/pck/earth_000101_230102_221009.bpc',
                       '/pck/earth_200101_990628_predict.bpc']


#run 
#load spiceypy kernels
cp.loadKernels(kernels_to_load)

p = cp.CHAOS(cubesat,                   #Satellite object  
                manager,                #perturbation manager class    
                LEO,                    #Environment object
                thruster_list,          #List of thruster object
                CS,                     #Control system object
                sensor,                 #Sensor object
                store,                  #Dataset object
                t_stop=stopper,         #Hard simulated time limit
                writepath=outputPath,   #Location of output data file
                iteration=rank,         #Simulation ID (int)     
                integrator='DOP853',    #Integration algorithm used
                tols=1e-6)              #Tolerance for the integrator


#unload kernels after the simulation
cp.unloadKernels(kernels_to_load)




#ploting code

cp.plotData(rank,                           #ID of the simulation to plot
            datapath=outputPath,            #Where to read the data files from
            figpath='./figs/')              #Where to store the figure (pdfs)