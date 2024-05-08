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
        self.egm = kwargs.get('egm', True)
        self.gg = kwargs.get('gravityGradient', False)
        self.lunarGravity = kwargs.get('lunarGravity', True)
        self.solarGravity = kwargs.get('solarGravity', True)
        self.srp = kwargs.get('SRP', True)

        #define constants
        self.atmTemperature = 1100 # Kelvin
        self.egmDegree = 2
        self.egmOrder = self.egmDegree
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

        #egm
        #get data path
        self.dataFilePath = os.path.dirname(__file__) + '/Data'

        self.C = matricise(2, filename = self.dataFilePath + '/Environment/Gravity/EGM96_n36') #EGM coefficient
        self.S = matricise(3, filename = self.dataFilePath + '/Environment/Gravity/EGM96_n36') #EGM coefficients
        #space weather
        self.interpf10, self.interpctr81 = interpF10Index_ext(self.dataFilePath + '/Environment/Atmosphere/SpaceWeather-test.txt')
        self.interpOions, self.interpHions = rhoIonsFunc()
        #temporary
        # self.interp = interpDensity('J77_spot/GMATresults_atmDensity.txt') # temporary interpolation over the GMAT density

    #define methods here

    def density(self, h, T_inf):
  
    # calculates the atmospheric density based on J77. 
    # h in km, from Earth surface. Not valid below 155km.
    # T_inf in Kelvin, between 650 and 1350 K.
    # inputs: 
    # h       altitude in km 
    # T_inf   temperature in Kelvin
        T_min = 650
        T_max = 1350
        num_partial_atm = 8 #number of partial atmospheres 
        deg  = 9 # degrees and order of the static parameters

        #ensure T_inf is bounded 
        if T_inf < T_min:
            T_inf = T_min

        elif T_inf > T_max:
            T_inf = T_max

        T_unitless = (T_inf - T_min) /  (T_max - T_min)
        H_p = []
        rho_p = []
        rho = 0
        T_matrix = np.array([T_unitless**0, T_unitless**1, T_unitless**2, T_unitless**3,T_unitless**4,T_unitless**5, T_unitless**6, T_unitless**7, T_unitless**8 ])
        
        a_para = np.matmul(a, T_matrix)
        b_para = np.matmul(b, T_matrix)

        H_p = -1/(a_para)
        rho_p = np.exp(b_para)
        
        for i in range(num_partial_atm):
 
            rho += rho_p[i] * np.exp(-h/H_p[i])
            
            
        return rho #kg/m^3 

    def drag(self, r_eci, v_eci,  m, arr_area, Ta, arr_Tw, quaternion, arr_normal, arr_alpha_coeff):
        
        # Function to calculate the acceleration due to drag at a given radius and atmospheric temperature.
        # The function computes the drag acting on a given face, hence area and normal of face mus be provided. 
        # If the face is not facing the correct direction, no drag is generated.
        #  
        # Acceleration given in km.
        # r:             position vector, in km, from Earth centre. 
        # T_inf:         in Kelvin, between 650 and 1350 K.
        # v:             velocity vector, in km 
        
   
        v = v_eci *1000 #in m/s
        vec_v_rel = v - np.cross(self.omega_E, r_eci*1000) # Include rotation of atmosphere
        v_rel = np.linalg.norm(vec_v_rel)
        v_dir = vec_v_rel / v_rel #velocity unit vector

        #set face normal in inertial frame
        arr_normal_in = body_to_inertial(arr_normal, quaternion)
        # print('quat nomr', np.linalg.norm(quaternion))
        
        arr_v_dir_lift = np.cross(np.cross(v_dir, arr_normal_in), v_dir)


        #need to normalise all the normals relative to their own magnitude

        #define inverse of norm of each row
        arr_norm_normal_in = 1 / np.linalg.norm(arr_normal_in, axis=1) 
        #multiply each normal vector by its inverse norm: normalisation
        arr_normal_in = np.einsum('i, ij->ij', arr_norm_normal_in, arr_normal_in) 
        #angle of normal relative to velocity direction
        dott = np.dot(arr_normal_in, v_dir)

        # print('norm of the nomrals', np.linalg.norm(arr_normal_in, axis=1))

        # print('dot', dott)
        # print('v_dir norm', np.linalg.norm(v_dir))
        arr_angle = np.arccos(dott)
        # print('arr_angle', arr_angle*180/np.pi)
        # print('area',area)
        #density
        h = (np.linalg.norm(r_eci) - self.R_E)
        # print('h', h)
        rho = self.density(h, Ta)
        
        #compute drag and lift coefficient, Aref=1
        arr_Cd, arr_Cl = MoeSentman(arr_angle, arr_area, 1, Ta, arr_Tw, v_rel, arr_alpha_coeff )
        
        #drg and lift forces --#area not present as Cd /Cl alerady scaled for it.
        drag = -  (0.5 * (arr_Cd) * rho * (v_rel**2)) / m  #SI
        lift =  (0.5 * (arr_Cl) * rho * (v_rel**2)) / m # m/s

        #in inertial frame
        arr_vec_drag = np.einsum('i, ij->ij', drag, np.array([v_dir])) #drag is same dir for all facets
        arr_vec_lift = np.einsum('i, ij->ij', lift, arr_v_dir_lift) #lift is normal to facet
        # print('angle', arr_angle*180/np.pi, 'Cd', arr_Cd, 'Cl', arr_Cl, 'drag', arr_vec_drag, 'lift', arr_vec_lift)

        #sum of the vectors
        arr_vec_a = arr_vec_drag + arr_vec_lift
        # print('a', arr_area)
        return arr_vec_a/1000  #acceleeration due to drag, in km/s^2
        

    def drag2(self, r_eci, v_eci, rho, area, Cd, m):
        # Function to calculate the acceleration due to drag at a given radius and atmospheric temperature.
        # The function computes the drag acting on a given face, hence area and normal of face must be provided. 
        # If the face is not facing the correct direction, no drag is generated.
        #  
        # Acceleration given in km.
        # r:             position vector, in km, from Earth centre. 
        # T_inf:         in Kelvin, between 650 and 1350 K.
        # v:             velocity vector, in km 
        
        v = v_eci *1000 #in m/s


        vec_v_rel = v - np.cross(self.omega_E, r_eci*1000) # Include rotation of atmosphere
        v_rel = np.linalg.norm(vec_v_rel)
        v_dir = vec_v_rel / v_rel #velocity unit vector


        vec_a = v_dir * (-0.5 * (Cd * area) * rho * (v_rel**2)) / m  #SI


        # print(vec_a/1000)
        return vec_a/1000  #acceleeration due to drag, in km/s^2


    def EGM(self, t, r_eci, v_eci):

        # print('reci', r_eci, 'veci',v_eci)
        a_egm = egm08(t, r_eci, v_eci, self.C, self.S, self.epoch, self.egmDegree, self.egmOrder)

        return a_egm # eci

    
    def moonGravity(self, t, r_sat, r_E_L):
        #Function to provide the perturbing acceleration due to the Moon a


        r_l_sat = -(r_sat - r_E_L)

        acc_Luna  = self.mu_moon * ( (r_l_sat / (np.linalg.norm(r_l_sat)**3)) - (r_E_L / (np.linalg.norm(r_E_L)**3)) )
        

        return acc_Luna

    def sunGravity(self, t, r_sat, r_E_Sun):
        #Function to provide the perturbing acceleration due to the Moon and Sun gravity.
        #The model uses the Astronomical Almanac to produce the position of the Moon and Sun with respect to the 
        #Earth. 


        # r_E_Sun = sunAlmanacVallado(t, self.epoch) #Geocentric position of the Sun
        # r_E_Sun = spiceSunPos(t, self.epoch, 'earth')
        r_s_sat = -(r_sat - r_E_Sun)

        acc_sun  = self.mu_sun * ( (r_s_sat / (np.linalg.norm(r_s_sat)**3)) - (r_E_Sun / (np.linalg.norm(r_E_Sun)**3)) )

        return acc_sun
    



    def percentShadow(self, vec_r, vec_r_sun):
        #usefule vectors

        vec_sat_sun = vec_r_sun - vec_r
        #Apparent radius of Sun and Earth 
        R_E_apparent = np.arcsin(self.R_E / np.linalg.norm(vec_r))
        R_S_apparent = np.arcsin(self.R_Sun / np.linalg.norm(vec_sat_sun))

        # print('err', vec_r_sun, vec_r)
        percentSun = 0
        #Apparent separation

        D_apparent = np.arccos( (- np.dot(vec_r, vec_sat_sun)) / (np.linalg.norm(vec_r)*np.linalg.norm(vec_sat_sun)) )
        
        #Illumination conditions
        if D_apparent >= (R_S_apparent + R_E_apparent):
            p = 0 # IN sun!

        elif D_apparent <= (R_E_apparent - R_S_apparent):
            # if R_E_apparent >= R_S_apparent: #Full occulatation
            p = 100 # IN SHADOW!!
             

        #Partial eclipse / Partial occultation
        elif D_apparent < (R_S_apparent + R_E_apparent) and D_apparent > abs(R_E_apparent - R_S_apparent):

            # Area overlap
            c_1 = ((D_apparent**2) + (R_S_apparent**2) - (R_E_apparent**2)) / (2 * D_apparent)
            c_2 = np.sqrt((R_S_apparent**2) - (c_1**2))
            A = (R_S_apparent**2) * np.arccos(c_1 / R_S_apparent) + (R_E_apparent**2) * np.arccos((D_apparent - c_1)/R_E_apparent) - D_apparent * c_2
            
            p = 100 * (A /(np.pi * (R_S_apparent**2 )))
        
        else: #Annular eclipse
                p = 100 * ((R_E_apparent**2) / (R_S_apparent**2))
        #Transform proportion of shadow in proportion of sun:
        percentSun = 100 - p

        return percentSun  #in percent!!!!


    def solarRadiationPressure(self, vec_r, vec_sun, m, quaternion, arr_area, arr_normal, arr_Cr):

        p = self.percentShadow(vec_r, vec_sun) #in/out of shadow

        sun_dir = vec_sun / np.linalg.norm(vec_sun) # unit vector of sun 
        vec_sun_sat = vec_r - vec_sun #sun to satellite vector, km 
        solarIrradiance = 3.844405613e26 / ((4 * np.pi) * (np.linalg.norm(vec_sun_sat*1000)**2)) #W/m^2         
        SolarPressure = solarIrradiance / self.lightSpeed #N/m^2


        #get panel normal in the inertial frame
        arr_normal_in = body_to_inertial(arr_normal, quaternion)
        dot_prod = np.dot(arr_normal_in, sun_dir) #dot product of 2 unit vectors


        #compute if facet is facing sun direction
        arr_facet_state = dot_prod > 0 #boolean array
        #projected area
        arr_A = arr_area * dot_prod  * arr_facet_state#null if facet is not facing the sun direction.


        #compute solar pressure acceleration
        #magnitude of force for every facet
        arr_srp = (p/100) * SolarPressure *(( arr_Cr * arr_A) / m) #SolarPressure
        #force direction (unit vector)
        vec_srp_dir = (vec_sun_sat / np.linalg.norm(vec_sun_sat))

        #array of force vectors (arr_srp[i] * vec_srp_dir)
        arr_vec_a = np.einsum('i, ij->ij', arr_srp, np.array([vec_srp_dir])) 



        #output is an array of force vectors
        return arr_vec_a/1000 #km/s

    def solarRadiationPressure2(self, vec_r, vec_sun, area, m, Cr):

        p = self.percentShadow(vec_r, vec_sun) #in/out of shadow
 
        vec_sun_sat = vec_r - vec_sun #sun to satellite vector, km 

        # print('shadow model', p)
        r_au = 149597870.691 #GMAT 1AU, KM

        solarIrradiance = 3.844405613e26 / ((4 * np.pi) * (np.linalg.norm(vec_sun_sat*1000)**2)) #W/m^2 
        solarPressure = solarIrradiance / self.lightSpeed #N/m^2


        A = area#* ( abs(np.dot(x_in, sun_dir)) + abs(np.dot(y_in, sun_dir)) + abs(np.dot(z_in, sun_dir)))
        # print('area', A)
        

          #vector between sun and satellite
        vec_a =  (p/100) * ((solarPressure * Cr * A) / m) * ((vec_sun_sat *1000)/ np.linalg.norm(vec_sun_sat*1000))

        # vec_a = (p/100) * solarPressure * Cr * (A / m) *((r_au*1000)**2)* (vec_sun_sat*1000 )/ (np.linalg.norm(vec_sun_sat*1000)**3)


        return vec_a/1000 #km/s


    def gravityGradient(self, vec_r_eci, quaternion, inertiaMatrix):
        radial = np.linalg.norm(vec_r_eci)
        nadir_eci = -vec_r_eci / radial #inertial frame
        nadir = inertial_to_body(nadir_eci, quaternion)
        A = np.matmul(inertiaMatrix, nadir)
        #torque:
        t_gravity = 3 * ( self.mu / (radial)**3 ) * np.cross(nadir, A)
    
        return t_gravity

    

    def orbitalPerturbation(self, t, vec_r_eci, vec_v_eci, area, m, Cd, Cr):
        
        #set vectors
        egm_perturb = np.array([0., 0., 0.])
        L_perturb = np.array([0., 0., 0.])
        S_perturb = np.array([0., 0., 0.])


        # print('t',t, 'r', vec_r_eci, 'v', vec_v_eci, 'A', area, 'm', m, 'CD', Cd, 'CR', Cr)
        #get the vectors from other classes
        vec_sun = spiceSunPos(t, self.epoch, 'earth') # km 
        vec_moon = spiceMoonPos(t, self.epoch, 'earth')



        #pure orbital perturbations
        if self.egm == True:
            # print('egm on')
            egm_perturb = self.EGM(t, vec_r_eci, vec_v_eci)

        if self.lunarGravity == True:
            # print('Luni on')
            L_perturb = self.moonGravity(t, vec_r_eci, vec_moon)

        if self.solarGravity == True:
            # print('Solar on')
            S_perturb = self.sunGravity(t, vec_r_eci, vec_sun)

        # if self.srp == True:
        #     # print('SRP on')
        #     SRP_perturb = self.solarRadiationPressure2(vec_r_eci, vec_sun, area, m, Cr)
                        
        # if self.J77 == True:
        #     # print('drag on')
        #     drag_perturb = self.drag2(vec_r_eci, vec_v_eci, rho, area, Cd, m)


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

        #TEST
        # p = self.percentShadow(vec_r_eci, vec_sun)
        # file = open('test_atm' + '.txt', 'a')
        # file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(t, F10, F10_avg, T, rho, p))
        # file.close()

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
if __name__ == '__main__':
    from class_satellite import Satellite
    import matplotlib.pyplot as plt



    # cubesat = Satellite()
    # cubesat.initial_induction()
    # drag = Environment(cubesat)
    # drag.J77 = 1
    # drag.egm = False
    # T_inf = 1000 # 
    # print(drag.drag())
    # print(drag.EGM(0))
    # print(drag.perturbation(0))
   
    # h = np.linspace(6378 + 100 , 6378 + 2500, 2000)
    # rho = []
    # drag_list1 = []
    # drag_list2 = []
    # drag_list3 = []
    # for i in range(len(h)):
    #     v = np.array([np.sqrt(3.986e5 / h[i]), 0, 0])
    #     # rho = np.append(rho, np.linalg.norm(drag.drag()))
    #     drag_list1 = np.append(drag_list1, np.linalg.norm(drag(np.array([h[i], 0, 0]), T_inf,v, 2.2, 0.015)))
    # #     # drag_list2 = np.append(drag_list2, np.linalg.norm(drag(np.array([h[i], 0, 0]), 700,v, 2.2, 0.015)))
    # #     # drag_list3 = np.append(drag_list3, np.linalg.norm(drag(np.array([h[i], 0, 0]), 1300,v, 2.2, 0.015)))

    
    # plt.plot(drag_list1, h-6378, label='density @ 1000K')
    # plt.xlabel(r'Atmospheric density (J77) $\left[\frac{kg}{m^3}\right]$')
    # plt.ylabel('Altitude [km]')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()