#python code to hold the extraction code as a class


import numpy as np

import os



# from .plotter import individual_timetracks, plot_histogram, plot_single_timetrack


# '''
# The Extractor class encapsulates all the procedures, data storage and manipulation required 
# for the analysis of parametric studies regarding the attitude simulator.
# It produces average timetrack plots and for each grid/perf scenarios, 
# and heatmaps to compare the different set-ups. 
# The class must be initialised with the data from the Monte-Carlo study: grid-size, perf-size, and
# number of iteration

# the methods available are:
#                 -read_data:             read textfile data as stored by the wrapper_function file

#                 -iteration_analysis:    handles the statistical analysis of the data, 
#                                         to produce the statistical data textfile

#                 -timetrack:             Produces the average timetrack plots. 
#                                         This mean taking the average value for each time point t, 
#                                         and plotting the point with a std dev     
# '''

# class Extractor:
#     '''
#     Class to handle the data analysis. Initialise with 
#     grid size (array)
#     perf value (array)
#     # iteration (int)
#     '''
#     #initialization upon instantiation
#     def __init__(self, grid_list, perf_list, iteration, **kwargs):
#          #create a bunch of list, dicts, etc for storage // similar to dataset
#         x_labels = kwargs.get('x_labels', None)
#         y_labels = kwargs.get('y_labels', None)

#         self.x_labels = x_labels
#         self.y_labels = y_labels
#         self.grid_list = grid_list
#         self.perf_list = perf_list
#         self.iteration = iteration

#         self.omega_x = []
#         self.omega_y = []
#         self.omega_z = []
#         self.omega_norm = []

#         self.time_dump = []

#         self.torque = []
#         self.torque_x = []
#         self.torque_y = []
#         self.torque_z = []
#         self.torque_norm = []

#         self.hmap_mean = {}
#         self.hmap_peak = {}
#         self.hmap_std = {}
#         self.hmap_impulse = {}
#         self.hmap_ratio = {}
#         self.hmap_ignition = {}
#         self.hmap_fuel = {}
#         self.hmap_firing_time = {}
        
        
#         self.maxTrack = 0

#         self.absolute_path = os.getcwd()
        
    
#     def read_data(self, filename):
#         t=[]
#         x=[]
#         y=[]
#         z=[]
#         norm=[]
#         datafile = open(filename, 'r')
#         for line in datafile.readlines(): #split each line
#             item = line.split('\t') #split each element
#             t.append(float(item[0]))
#             x.append(float(item[1]))
#             y.append(float(item[2]))
#             z.append(float(item[3]))
#             norm.append(float(item[4]))

#         return t, x, y, z, norm

    

#     def iteration_analysis(self):
#         #for each simluation and iteration: gather data:
            
#         for grid_num in range(len(self.grid_list)):
#             grid_name = 'grid{}x{}'.format(self.grid_list[grid_num], self.grid_list[grid_num])

#             #create heatmap entry for eah grid
#             self.hmap_mean['{}'.format(self.grid_list[grid_num])] = []
#             self.hmap_std['{}'.format(self.grid_list[grid_num])] = []
#             self.hmap_impulse['{}'.format(self.grid_list[grid_num])] = []
#             self.hmap_ratio['{}'.format(self.grid_list[grid_num])] = []
#             self.hmap_ignition['{}'.format(self.grid_list[grid_num])] = []
#             self.hmap_firing_time['{}'.format(self.grid_list[grid_num])] = []
#             self.hmap_fuel['{}'.format(self.grid_list[grid_num])] = []
#             self.hmap_peak['{}'.format(self.grid_list[grid_num])] = []

#             for perf_num in range(len(self.perf_list)):
#                 perf_name = grid_name + '/perf{}'.format(self.perf_list[perf_num])
#                 if os.path.exists(perf_name + '/Plots') == False:
#                     os.mkdir(perf_name + '/Plots')

#                 self.top_tx = []
#                 self.top_ty = []
#                 self.top_tz= []
#                 self.top_t_norm = []

#                 self.rotation = []
#                 self.omega_stat = []
#                 self.w_x_stat = []
#                 self.w_y_stat = []
#                 self.w_z_stat = []

#                 self.tx_stat = []
#                 self.ty_stat = []
#                 self.tz_stat = []

#                 self.final_wx = []
#                 self.final_wy = []
#                 self.final_wz = []
#                 self.final_w_norm = []

#                 self.averageIgnition = []
#                 self.averageFiringTime = []
#                 self.averageRemainingFuel = []
#                 self.maxVelocity = []
#                 #for each iteration
#                 for i in range(self.iteration):
#                     self.sequence, self.torque_x, self.torque_y, self.torque_z, self.torque_norm = self.read_data(perf_name + '/stat/StatData{}.txt'.format(i))
#                     self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(i))
#                     self.meanIgnition, self.totalFuel, self.maxTrack, self.meanFiringTime, self.redundant = self.read_data(perf_name + '/gauss/RateData{}.txt'.format(i))
#                     self.top_tx = np.append(self.top_tx, np.sum(self.torque_x))
#                     self.top_ty = np.append(self.top_ty, np.sum(self.torque_y))
#                     self.top_tz = np.append(self.top_tz, np.sum(self.torque_z))
#                     self.top_t_norm = np.append(self.top_t_norm, np.sum(self.torque_norm))

#                     self.final_wx = np.append(self.final_wx, self.omega_x[-1])
#                     self.final_wy = np.append(self.final_wy, self.omega_y[-1])
#                     self.final_wz = np.append(self.final_wz, self.omega_z[-1])
#                     self.final_w_norm = np.append(self.final_w_norm, self.omega_norm[-1])
#                     self.rotation = self.final_w_norm / (2 * np.pi)

#                     self.averageIgnition = np.append(self.averageIgnition, self.meanIgnition)
#                     self.averageFiringTime = np.append(self.averageFiringTime, self.meanFiringTime[0]) #in sec
#                     self.averageRemainingFuel = np.append(self.averageRemainingFuel, self.totalFuel)
#                     self.maxVelocity = np.append(self.maxVelocity, self.maxTrack)

#                 #compute statistical parameters:
#                 self.rot_stat = [np.mean(self.rotation), np.std(self.rotation)]
#                 self.omega_stat =[np.mean(self.final_w_norm), np.std(self.final_w_norm)]
#                 self.w_z_stat = [np.mean(self.final_wz), np.std(self.final_wz)]
#                 self.w_y_stat = [np.mean(self.final_wy), np.std(self.final_wy)]
#                 self.w_x_stat = [np.mean(self.final_wx), np.std(self.final_wx)]
#                 self.tz_stat = [np.mean(self.top_tz), np.std(self.top_tz)]
#                 self.ty_stat = [np.mean(self.top_ty), np.std(self.top_ty)]
#                 self.tx_stat = [np.mean(self.top_tx), np.std(self.top_tx)]
#                 self.ignition_stat = [np.mean(self.averageIgnition), np.std(self.averageIgnition)]
#                 self.firing_stat = [np.mean(self.averageFiringTime), np.std(self.averageFiringTime)]
#                 self.fuel_stat = [np.mean(self.averageRemainingFuel), np.std(self.averageRemainingFuel)]
#                 self.peak_stat = [np.mean(self.maxVelocity), np.std(self.maxVelocity)]

#                 #print data to file:
#                 #delete previous file:
#                 if os.path.exists(perf_name + '/Plots/final velocity statistical data'):
#                     os.remove(perf_name + '/Plots/final velocity statistical data')
                
#                 #concise writing
#                 statistical_data = {'rotation':self.rot_stat, '|omega|':self.omega_stat, 'total wz':self.w_z_stat,'total wy':self.w_y_stat, 'total wx':self.w_x_stat, 'total tz':self.tz_stat, 'total ty':self.ty_stat, 'total tx':self.tx_stat}
#                 parameters = [self.rotation, self.omega_norm, self.omega_z, self.omega_y, self.omega_x, self.torque_z, self.torque_y, self.torque_x]
                
#                 #fill the heatmaps
#                 self.hmap_mean['{}'.format(self.grid_list[grid_num])] = np.append(self.hmap_mean['{}'.format(self.grid_list[grid_num])],  statistical_data['|omega|'][0])
#                 self.hmap_std['{}'.format(self.grid_list[grid_num])] = np.append(self.hmap_std['{}'.format(self.grid_list[grid_num])], statistical_data['|omega|'][1])
#                 self.hmap_ignition['{}'.format(self.grid_list[grid_num])] = np.append(self.hmap_ignition['{}'.format(self.grid_list[grid_num])], self.ignition_stat[0])
#                 self.hmap_fuel['{}'.format(self.grid_list[grid_num])] = np.append(self.hmap_fuel['{}'.format(self.grid_list[grid_num])], self.fuel_stat[0])
#                 self.hmap_firing_time['{}'.format(self.grid_list[grid_num])] = np.append(self.hmap_firing_time['{}'.format(self.grid_list[grid_num])], self.firing_stat[0])
#                 self.hmap_peak['{}'.format(self.grid_list[grid_num])] = np.append(self.hmap_peak['{}'.format(self.grid_list[grid_num])], max(self.maxVelocity))
                
                
#                 #write the statistical data file 
#                 for i in range(len(statistical_data)):
#                     keyname = list(statistical_data.keys())[i]              
#                     f = open(perf_name + '/Plots/final velocity statistical data', 'a')    
#                     f.write('{} \n max value \t {} \n min value \t {} \n mean value \t {} \n standard deviation \t {} \n'.format( keyname, max(parameters[i]), min(parameters[i]), statistical_data[keyname][0], statistical_data[keyname][1] ))
#                     f.close()
                
#                 # create folder for gaussian plots
#                 if os.path.exists(perf_name + '/Plots/gaussian') == False:
#                     os.mkdir(perf_name + '/Plots/gaussian')
                
#                 nbins = 50  #number of bins in histograms
#                 plot_histogram(nbins, self.averageFiringTime, self.firing_stat, 'Min firing time distribution', 'quadrant firing time', 'firing', perf_name + '/Plots/gaussian/', save=0, show=1)        
#                 plot_histogram(nbins, self.averageIgnition, self.ignition_stat, 'Max ignition distribution', '# of ignition', 'ignition',  perf_name + '/Plots/gaussian/', save=0, show=1)        
#                 plot_histogram(nbins, self.averageRemainingFuel, self.fuel_stat, 'Used fuel distribution', '\% of remaining fuel', 'fuel', perf_name + '/Plots/gaussian/', save=0, show=1)        
#                 plot_histogram(nbins, self.maxVelocity, self.peak_stat, 'Peak velocity distribution', 'Peak velocity', 'peak_dist', perf_name + '/Plots/gaussian/', save=0, show=1)


   
#     def timetrack(self):
#         import time
#         #for each grid
#         for grid_num in range(len(self.grid_list)):
#             grid_name = 'grid{}x{}'.format(self.grid_list[grid_num], self.grid_list[grid_num])
#             self.hmap_peak['{}'.format(self.grid_list[grid_num])] = []
#             print('{}'.format(grid_name))
#             time_timetrack_start = time.time()
#             #store this information
            

#             #for each performance
#             for perf_num in range(len(self.perf_list)):
#                     perf_name = grid_name + '/perf{}'.format(self.perf_list[perf_num])
#                     print('\t \t {}'.format(perf_name))

#                     #define dicts for storage and sata manipulation
#                     self.omega_x_mean = []
#                     self.omega_y_mean = []
#                     self.omega_z_mean = []
#                     self.omega_norm_mean = []

#                     self.omega_x_err = []
#                     self.omega_y_err = []
#                     self.omega_z_err = []
#                     self.omega_norm_err = []
#                     self.peakVelocityList = []


#                     x_vel = {}
#                     y_vel = {}
#                     z_vel = {}
#                     counter = 0           

#                     #get number of timetrack files
#                     timefile_name = os.listdir(perf_name + '/timetrack')

#                     #initialise max velocity
#                     self.maxVelocity = [0]
#                     #for each file
#                     for file in timefile_name:

#                         #initialise the array.
#                         x_vel['{}'.format(counter)] = []
#                         y_vel['{}'.format(counter)] = [] 
#                         z_vel['{}'.format(counter)] = []   

#                         #read each file:
#                         self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(str(counter)))

#                         #dict containing each elements x, y, z for each timetrack file
#                         x_vel['{}'.format(counter)] = np.append(x_vel['{}'.format(counter)], self.omega_x) #  every 1for solve_ivp every 100 point for odeint,
#                         y_vel['{}'.format(counter)] = np.append(y_vel['{}'.format(counter)], self.omega_y)
#                         z_vel['{}'.format(counter)] = np.append(z_vel['{}'.format(counter)], self.omega_z)
                        
#                         #move 1 file up
#                         counter += 1
                    
#                     #for each time element:
                    
#                     for time_elmt in range(len(self.time_dump)):
#                         y_velocity_at_time = []
#                         x_velocity_at_time = []
#                         z_velocity_at_time = []
#                         norm_velocity_at_time = []        

#                         for file_elmt in range(counter):
#                             x_velocity_at_time = np.append(x_velocity_at_time, x_vel['{}'.format(file_elmt)][time_elmt]) #vertical x-velocity vaues
#                             y_velocity_at_time = np.append(y_velocity_at_time, y_vel['{}'.format(file_elmt)][time_elmt]) #vertical y-velocity values
#                             z_velocity_at_time = np.append(z_velocity_at_time, z_vel['{}'.format(file_elmt)][time_elmt]) #vertical z-velocity values
#                             at_this_time_norm = np.sqrt(z_vel['{}'.format(file_elmt)][time_elmt]**2 + y_vel['{}'.format(file_elmt)][time_elmt]**2 + x_vel['{}'.format(file_elmt)][time_elmt]**2)
#                             norm_velocity_at_time = np.append(norm_velocity_at_time, at_this_time_norm) #vertical norm-velocity values
#                         #calculate mean 
#                         self.omega_x_mean = np.append(self.omega_x_mean, np.mean(x_velocity_at_time)) 
#                         self.omega_y_mean = np.append(self.omega_y_mean, np.mean(y_velocity_at_time))
#                         self.omega_z_mean = np.append(self.omega_z_mean, np.mean(z_velocity_at_time))
#                         self.omega_norm_mean = np.append(self.omega_norm_mean, np.mean(norm_velocity_at_time))
#                         #calculate peak velocity
#                         # self.peakVelocityList = np.append(self.peakVelocityList, max(norm_velocity_at_time))
#                         # if max(norm_velocity_at_time) > self.maxVelocity:
#                         #     self.maxVelocity = max(norm_velocity_at_time)
#                         # calculate std dev
#                         self.omega_x_err.append(np.std(x_velocity_at_time))
#                         self.omega_y_err.append(np.std(y_velocity_at_time))
#                         self.omega_z_err.append(np.std(z_velocity_at_time))
#                         self.omega_norm_err.append(np.std(norm_velocity_at_time))       

#                     #calculate equivalent average time:
#                     time_average = np.mean(self.omega_norm_mean)
                    
#                     #fill dict for heatmap
#                     # self.hmap_peak['{}'.format(self.grid_list[grid_num])] = np.append(self.hmap_peak['{}'.format(self.grid_list[grid_num])], self.maxVelocity)

#                     #write data in file:
#                     f = open(perf_name + '/Plots/final velocity statistical data', 'a')    
#                     f.write('Max track angular velocity norm \t {} \n Time averaged velocity \t {} \n'.format(self.maxVelocity, time_average))
#                     f.close()

#                     # if not there, create the storing folder:
#                     if os.path.exists(perf_name + '/Plots/timetrack') == False:
#                         os.mkdir(perf_name + '/Plots/timetrack')

#                     #plot the time tracks:
#                     individual_timetracks(self.time_dump, *[self.omega_x_mean, self.omega_x_err, self.omega_y_mean, self.omega_y_err, self.omega_z_mean, self.omega_z_err, self.omega_norm_mean, self.omega_norm_err], save_path=perf_name + '/Plots/timetrack/', save=0, show=0)

#                     #plot distribution of peak velocity
#                     # self.peak_stat = [np.mean(self.peakVelocityList), np.std(self.peakVelocityList)]
#                     # nbins = 50  #number of bins in histograms
#                     # plot_histogram(nbins, self.peakVelocityList, self.peak_stat, 'Peak velocity distribution', 'Peak velocity', 'peak_dist', self.absolute_path + '/Plots/stat/', save=0, show=1)
#             time_timetrack_end = time.time()

#     def typical_time_profile(self, example=0):
#         #for each grid
#         import time
#         for grid_num in range(len(self.grid_list)):
#             grid_name = 'grid{}x{}'.format(self.grid_list[grid_num], self.grid_list[grid_num])
#             print('{}'.format(grid_name))
        
#             time_timetrack_start = time.time()
#             #store this information
            

#             #for each performance
#             for perf_num in range(len(self.perf_list)):
#                     import matplotlib.pyplot as plt

#                     perf_name = grid_name + '/perf{}'.format(self.perf_list[perf_num])
#                     print('\t \t {}'.format(perf_name))
#                     if os.path.exists(perf_name + '/Plots/timetrack') == False:
#                         os.mkdir(perf_name + '/Plots/timetrack')
#                     #initialise counter
#                     counter = 0           
                    
#                     #get number of timetrack files
#                     timefile_name = os.listdir(perf_name + '/timetrack')
#                     #plot simplistic of all timetracks
#                     # for the norm:
#                     #initialise counter
#                     counter = 0  
#                     for file in timefile_name:
#                         self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(str(counter)))
#                         plt.plot(self.time_dump, self.omega_norm, color='#D3D3D3')
#                         counter +=1
#                     self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(example))
#                     plot_single_timetrack(self.time_dump, self.omega_norm, perf_name + '/Plots/timetrack/norm.pdf', 'Time evolution of the angular velocity norm', 'norm', show=0, save=1)

#                     # # for y values
#                     # # initialise counter
#                     # counter = 0  
#                     # for file in timefile_name:
#                     #     self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(str(counter)))
#                     #     plt.plot(self.time_dump[::100], self.omega_y[::100], color='#D3D3D3')
#                     #     counter +=1
                    
#                     # self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(example))
#                     # plot_single_timetrack(self.time_dump, self.omega_y, perf_name + '/Plots/timetrack/omega_y.pdf', 'Time evolution of the angular velocity y-component', 'y-component', show=1, save=1)

#                     # #for z values:
#                     # #initialise counter
#                     # counter = 0  
#                     # for file in timefile_name:
#                     #     self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(str(counter)))
#                     #     plt.plot(self.time_dump[::100], self.omega_z[::100], color='#D3D3D3')
#                     #     counter +=1
                    
#                     # self.time_dump, self.omega_x, self.omega_y, self.omega_z, self.omega_norm = self.read_data(perf_name + '/timetrack/timetrack{}.txt'.format(example))
#                     # plot_single_timetrack(self.time_dump, self.omega_z, perf_name + '/Plots/timetrack/omega_z.pdf', 'Time evolution of the angular velocity z-component', 'z-component', show=1, save=1)
                


def read_data(filename):
    t=[]
    x=[]
    y=[]
    z=[]
    norm=[]
    datafile = open(filename, 'r')
    for line in datafile.readlines(): #split each line
        item = line.split('\t') #split each element
        t.append(float(item[0]))
        x.append(float(item[1]))
        y.append(float(item[2]))
        z.append(float(item[3]))
        norm.append(float(item[4]))

    return t, x, y, z, norm


def read_grid_data(filename):
    t=[]
    x=[]
    y=[]
    z=[]
    norm=[]
    datafile = open(filename, 'r')
    for line in datafile.readlines(): #split each line
        item = line.split('\t') #split each element
        t.append(float(item[0]))
        x.append(item[1])
        y.append(float(item[2]))
        z.append(float(item[3]))
        norm.append(float(item[4]))

    return t, x, y, z, norm

