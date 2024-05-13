#Python code to hold the orbit perturbation class.


'''
Python scrip to hold the environment class definition. The inputs define the activation/ details of certain force models (atmospheric temperature, 
spherical harmonics degree, Epoch of Sun/Moon, etc..)
, and the methods
compute the acceleration/torques due to a given perturbation. 
'''

import numpy as np
from .EnvironmentModelling import atmTemperature, interpF10Index_ext, rhoIonsFunc, MoeSentman
from .J77_coeffs import a, b
from .EGM_V import egm08, matricise
from .SPICEroutine import spiceMoonPos, spiceSunPos
# from celestial_motion import moonAlmanac, sunAlmanac, sunAlmanacVallado
# from misc_function import setJulia
from .quaternions import body_to_inertial, inertial_to_body

# from misc_function import  interpDensity, readSpaceWeatherFile
import os
from .transformations import LatLong, julianDate

class Environment:

    def __init__(self, **kwargs):

        #define initial atttributes
        self.J77 = kwargs.get('J77', True)

        self.srp = kwargs.get('SRP', True)

        #define constants
        self.atmTemperature = 1100 # Kelvin

        self.epoch = np.array([1, 1, 2000])
        self.epoch_num = self.epoch

        day, month, year = self.epoch_num
        self.JD0 = julianDate(0, day, month, year) #needs updating!!
        #constants
        ##
        self.epoch_limit = np.array([20, 4, 2020])
        self.lightSpeed = 299792458 #m/s
        # self.SolarPressure = 4.57e-6 #N/m^2
        #planetary constants
        self.R_E = 6378.1363 # km
        self.R_Sun = 696000 #km
        self.mu = 398600.4418   #km-5-vallado
        self.mu_sun = 1.327122e11#1.32712428e11 #km^3 /s^2 gmat:132712440017.99
        self.mu_moon = 4.902800305555e3#4902.799 #km^3 /s2 GMAT: 4902.8005821478
        self.omega_E = np.array([0, 0, 0.0000729211585530]) #earth rotation
        self.EarthFlattening = 1 / 298.257


        #space weather
        self.interpf10, self.interpctr81 = interpF10Index_ext(self.dataFilePath + '/Environment/Atmosphere/SpaceWeather-test.txt')
        self.interpOions, self.interpHions = rhoIonsFunc()


        self.LFperturbFunction = []
        self.HFperturbFunction = []
    #define methods here


    

    def orbitalPerturbation(self, t, vec_r_eci, vec_v_eci, area, m, Cd, Cr):
        
        #set vectors
        egm_perturb = np.array([0., 0., 0.])
        L_perturb = np.array([0., 0., 0.])
        S_perturb = np.array([0., 0., 0.])




        #pure orbital perturbations

        egm_perturb = self.EGM(t, vec_r_eci, vec_v_eci)

        L_perturb = self.moonGravity( t, vec_r_eci, vec_v_eci, area, m, Cd, Cr)

        S_perturb = self.sunGravity( t, vec_r_eci, vec_v_eci, area, m, Cd, Cr)




        total_perturb = egm_perturb + L_perturb + S_perturb 


        return total_perturb


    def LFpanelPerturbation(self, t, vec_r_eci, vec_v_eci, area, m, Cd, Cr):

        #for GG, inertiaMatrix and Quaternions need to be passed. But then, the ORBITAL_DRIVER function needs it too...
        # Problematic. Create separate function for torques?
        drag_perturb = np.array([0., 0., 0.])
        SRP_perturb = np.array([0., 0., 0.])



        # print('t',t, 'r', vec_r_eci, 'v', vec_v_eci, 'A', area, 'm', m, 'CD', Cd, 'CR', Cr)
        #get the vectors from other classes
        vec_sun = spiceSunPos(t, self.epoch, 'earth') # km 


        ##temp
        r, lat, long = LatLong(vec_r_eci)

        #get F10.7 values
        day, month, year = self.epoch
        self.JD0 = julianDate(0, day, month, year)
        JD_t = self.JD0 + t / 86400
        F10 = self.interpf10(JD_t)
        F10_avg = self.interpctr81(JD_t)
        T = atmTemperature(F10, F10_avg, vec_sun, vec_r_eci, self.EarthFlattening)


        #call density functions
        h = (r - self.R_E) # km
        rho = self.density(h, T)

        #pure orbital perturbations

        if self.srp == True:
            # print('SRP on')
            SRP_perturb = self.solarRadiationPressure2(vec_r_eci, vec_sun, area, m, Cr)
                        
        if self.J77 == True:
            # print('drag on')
            drag_perturb = self.drag2(vec_r_eci, vec_v_eci, rho, area, Cd, m)


        total_perturb = drag_perturb + SRP_perturb


        return total_perturb


    def panelPerturbation(self, t, vec_r_eci, vec_v_eci, m, quaternion, arr_area, arr_normal, arr_Tw, arr_alpha_coeff, arr_Cr, arr_position):
        drag_perturb = np.array([0., 0., 0.])
        SRP_perturb = np.array([0., 0., 0.])
        SRP_torque = np.array([0., 0., 0.])
        aero_torque = np.array([0., 0., 0.])

                
        #get the vectors 
        vec_sun = spiceSunPos(t, self.epoch, 'earth') # km 


        #get F10.7 values
        day, month, year = self.epoch
        self.JD0 = julianDate(0, day, month, year)
        JD_t = self.JD0 + t / 86400
        F10 = self.interpf10(JD_t)
        F10_avg = self.interpctr81(JD_t)
        T = atmTemperature(F10, F10_avg, vec_sun, vec_r_eci, self.EarthFlattening)


        #attitude perturbations
        #Force computation
        if self.srp == True:
            # print('panel srp')
            #array of ACC vector--inertial frame, km/s^2
            arr_SRP_perturb = self.solarRadiationPressure( vec_r_eci, vec_sun, m, quaternion, arr_area, arr_normal, arr_Cr)
            
            # torques: force on each facet times the position of each facet
            #get force in SI units:
            arr_SRP_force = m * arr_SRP_perturb*1000
            #rotate force in body-frame
            body_frame_SRP_force = inertial_to_body(arr_SRP_force, quaternion) 
            #compute body-fixed torques
            arr_SRP_torque = np.cross(arr_position, body_frame_SRP_force)
            
            #sums all the forces to get the total force /torque
            SRP_perturb = np.einsum('ij->j', arr_SRP_perturb)
            SRP_torque = np.einsum('ij->j', arr_SRP_torque)

        if self.J77 == True:
            #Array of acceleration vector, inertial frame, km/s^2
            arr_drag_perturb = self.drag(vec_r_eci, vec_v_eci, m, arr_area, T, arr_Tw, quaternion, arr_normal, arr_alpha_coeff)
            
            # torques: force on each facet times the position of each facet
            #get force in SI units
            arr_aero_force = m* arr_drag_perturb*1000
            #rotate inertial forces into body frame
            body_frame_aero_force = inertial_to_body(arr_aero_force, quaternion) 
            arr_aero_torque = np.cross(arr_position, body_frame_aero_force)
            
            #sums all the forces to get the total force /torque
            drag_perturb = np.einsum('ij->j', arr_drag_perturb)
            aero_torque = np.einsum('ij->j', arr_aero_torque)    
            # print('drag acc, km/s', arr_drag_perturb)
            # print('drag torque', arr_aero_torque)

        total_perturb = drag_perturb + SRP_perturb
        total_torque = aero_torque + SRP_torque


        #returns the total orbital and torque perturbations
        return total_perturb, total_torque


    def ionDensity(self, vec_r):
        #compute altitude
        h = np.linalg.norm(vec_r) - self.R_E

        #compute ion density
        H_density = self.interpHions(h)*1e6 #from cm3 to m3
        O_density = self.interpOions(h)*1e6

        #return total ion density
        return H_density + O_density



