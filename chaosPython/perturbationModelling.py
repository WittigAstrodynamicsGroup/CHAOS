"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


peturbationsModelling.py


 This Python code defines functions for modeling various spacecraft perturbations.

 These functions are designed to estimate the environmental forces and 
 torques acting on a spacecraft orbiting the Earth.


 The provided functions address various perturbation sources:
 - density: calculates atmospheric density using the Jacchia-77 model.
 - LFdrag: estimates low-fidelity drag acceleration due to atmosphere.
 - nGravity: computes gravitational acceleration from a point mass.
 - percentShadow: determines the percentage of the Sun's disk visible from the spacecraft (for SRP calculations).
 - LFsolarRadiationPressure: estimates low-fidelity solar radiation pressure acceleration.
 - drag: calculates total aerodynamic acceleration and torque due to drag and lift (higher fidelity model).
 - solarRadiationPressure: computes total solar radiation pressure acceleration and torque acting on the spacecraft.
"""
import numpy as np
from .J77_coeffs import a, b
from .quaternions import body_to_inertial, inertial_to_body
from .EnvironmentModelling import MoeSentman




def density(h, T_inf):
    """
    Calculates atmospheric density using the Jacchia77 atmospheric model.

    This function estimates the atmospheric density (kg/m^3) at a given 
    altitude (`h` in km) based on the Jacchia77 model. The model assumes a 
    temperature at infinity (`T_inf` in Kelvin) between 650K and 1350K.

    **Note:** The model is not valid below 155 km altitude.

    Args:
        h (float): Altitude above Earth's surface in kilometers.
        T_inf (float): Temperature OF atmosphere in Kelvin (between 650K and 1350K).

    Returns:
        float: Atmospheric density at the specified altitude (kg/m^3).
    """


    #Temperature bounds of the model
    T_min = 650
    T_max = 1350

    #model details
    num_partial_atm = 8 #number of partial atmospheres 


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
    

    #density
    for i in range(num_partial_atm):

        rho += rho_p[i] * np.exp(-h/H_p[i])
        
    #return density
    return rho                                  #kg/m^3 



def LFdrag(r_eci, v_eci, rho, area, Cd, m, omega_E):

    """
    Calculates the spacecraft's acceleration due to atmospheric drag (Low-Fidelity).

    This function implements a simplified model to estimate the drag 
    acceleration acting on a spacecraft at a given position (`r_eci` in km) 
    and velocity (`v_eci` in km/s) in Earth-Centered Inertial (ECI) frame. 

    **Note:** The model assumes a drag force acting on a single face based 
    on the provided area and Cd. It does not account for drag on other 
    surfaces or complex vehicle geometry.

    Args:
        r_eci (numpy.ndarray): Spacecraft position vector in ECI frame (km).
        v_eci (numpy.ndarray): Spacecraft velocity vector in ECI frame (km/s).
        rho (float): Atmospheric density (kg/m^3).
        area (float): Cross-sectional area facing the relative wind (m^2).
        Cd (float): Drag coefficient.
        m (float): Spacecraft mass (kg).
        omega_E (float): Earth's rotation rate (rad/s).

    Returns:
        numpy.ndarray: Spacecraft acceleration due to drag in ECI frame (km/s^2).
    """


    #compute velocity relative to rotating atmosphere
    v = v_eci *1000                                         #in m/s
    vec_v_rel = v - np.cross(omega_E, r_eci*1000)           # Include rotation of atmosphere
    
    v_rel = np.linalg.norm(vec_v_rel)                       #norm
    v_dir = vec_v_rel / v_rel                               #velocity unit vector


    vec_a = v_dir * (-0.5 * (Cd * area) * rho * (v_rel**2)) / m         #acceleraiton due to drag----SI


    return vec_a/1000  #acceleration due to drag, in km/s^2



def nGravity(r_sat, r_E_n, mu_obj):
    """
    Calculates the gravitational acceleration on a satellite due to a point mass.

    This function computes the gravitational acceleration acting on a 
    satellite (at position `r_sat` in km) due to the gravitational 
    attraction of a point mass (located at `r_E_n` in km) with gravitational 
    parameter `mu_obj` (km^3/s^2). It assumes a point mass approximation 
    for the object's gravity.

    Args:
        r_sat (numpy.ndarray): Satellite position vector in km.
        r_E_n (numpy.ndarray): Position vector of the point mass in km.
        mu_obj (float): Gravitational parameter of the point mass (km^3/s^2).

    Returns:
        numpy.ndarray: Gravitational acceleration on the satellite due to 
                        the point mass (km/s^2).
    """

    #satellite to point mass vector
    r_n_sat = -(r_sat - r_E_n)

    ##acceleration computation
    acc_n  = mu_obj * ( (r_n_sat / (np.linalg.norm(r_n_sat)**3)) - (r_E_n / (np.linalg.norm(r_E_n)**3)) )
    

    return acc_n






def percentShadow(vec_r, vec_r_sun, radiusEarth, radiusSun):
    """
    Calculates the apparent percentage of the Sun's disk visible from a satellite position.

    This function considers the geometry of the Sun, Earth, and the satellite 
    to determine the portion of the Sun's disk that is occulted (blocked) by 
    the Earth as seen from the satellite's viewpoint. It accounts for full 
    illumination, complete shadow, partial eclipses (partial occultation 
    and partial illumination), and annular eclipses.

    Args:
        vec_r (numpy.ndarray): Satellite position vector in km.
        vec_r_sun (numpy.ndarray): Sun position vector in km.
        radiusEarth (float): Radius of the Earth in km.
        radiusSun (float): Radius of the Sun in km.

    Returns:
        float: Percentage of the Sun's disk visible from the satellite (%).
                - 0%: Full illumination (no shadow)
                - 100%: Complete shadow
                - Values between 0% and 100% indicate partial visibility 
                  due to eclipses.
    """


    #satellite to sun vector
    vec_sat_sun = vec_r_sun - vec_r


    #Apparent radius of Sun and Earth 
    R_E_apparent = np.arcsin(radiusEarth / np.linalg.norm(vec_r))
    R_S_apparent = np.arcsin(radiusSun / np.linalg.norm(vec_sat_sun))

    percentSun = 0

    #Apparent separation
    D_apparent = np.arccos( (- np.dot(vec_r, vec_sat_sun)) / (np.linalg.norm(vec_r)*np.linalg.norm(vec_sat_sun)) )
    
    #Illumination conditions
    if D_apparent >= (R_S_apparent + R_E_apparent):         #if apparent separation is greater than combied Sun-earth radii
        p = 0 # In Sun!

    elif D_apparent <= (R_E_apparent - R_S_apparent):
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
    percentSun = 100 - p                #"How much of the sun disk is visible?"

    return percentSun  #in percent!!!!




def LFsolarRadiationPressure(vec_r, vec_sun, area, m, Cr, lightSpeed, centralBodyRadius, starRadius):
    """
    Calculates the solar radiation pressure acceleration acting on a satellite.

    This function estimates the force exerted on a satellite due to solar 
    radiation pressure (SRp). 

    **Note:** This is a simplified model and does not account for all 
    complexities of SRP, such as spacecraft geometry and radiation pressure 
    variations.

    Args:
        vec_r (numpy.ndarray): Satellite position vector in km.
        vec_sun (numpy.ndarray): Sun position vector in km.
        area (float): Cross-sectional area facing the Sun (m^2).
        m (float): Satellite mass (kg).
        Cr (float): Solar radiation coefficient.

    Returns:
        numpy.ndarray: Solar radiation pressure acceleration on the satellite 
                        (km/s^2).
    """



    #Percent of sun disk visible
    p = percentShadow(vec_r, vec_sun, centralBodyRadius, starRadius)               #in/out of shadow


    vec_sun_sat = vec_r - vec_sun                   #sun to satellite vector, km 


    #compute solar irradiance at satellite position
    solarIrradiance = 3.844405613e26 / ((4 * np.pi) * (np.linalg.norm(vec_sun_sat*1000)**2)) #W/m^2 
    solarPressure = solarIrradiance / lightSpeed #N/m^2


    A = area
    #vector between sun and satellite
    vec_a =  (p/100) * ((solarPressure * Cr * A) / m) * ((vec_sun_sat *1000)/ np.linalg.norm(vec_sun_sat*1000))

    return vec_a/1000                               #km/s






def drag(r_eci, v_eci,  m, omega_E, R_E, arr_area, Ta, arr_Tw, quaternion, arr_normal, arr_position, arr_alpha_coeff):
    
    """
    Calculates the total aerodynamic acceleration and torque acting on a spacecraft due to atmospheric drag and lift.

    This function implements a more complex model for aerodynamic forces 
    compared to `LFdrag`. It considers multiple facets (surfaces) on the 
    spacecraft and calculates the drag and lift forces acting on each facet 
    based on its orientation relative to the airflow (velocity). The forces 
    are then summed to obtain the total aerodynamic acceleration and torque 
    acting on the spacecraft's center of mass.

    NOTE: Model valid only for convex geometries.

    Args:
        r_eci (numpy.ndarray): Spacecraft position vector in ECI frame (km).
        v_eci (numpy.ndarray): Spacecraft velocity vector in ECI frame (km/s).
        m (float): Spacecraft mass (kg).
        omega_E (float): Earth's rotation rate (rad/s).
        R_E (float): Earth's radius (km).
        arr_area (numpy.ndarray): Array of areas (m^2) for each spacecraft facet.
        Ta (float): Atmospheric temperature (Kelvin).
        arr_Tw (numpy.ndarray, optional): Array of wall temperatures (Kelvin) 
                                          for each facet.
        quaternion (numpy.ndarray): Spacecraft attitude quaternion.
        arr_normal (numpy.ndarray): Array of unit normal vectors in body frame 
                                    for each facet.
        arr_position (numpy.ndarray): Array of positions (meters) of each facet's 
                                      center of mass relative to the spacecraft's 
                                      center of mass.
        arr_alpha_coeff (numpy.ndarray): Array of alpha coefficients for the 
                                         aerodynamic force model.

    Returns:
        tuple:
            - total_aero_acceleration (numpy.ndarray): Total aerodynamic 
                                                      acceleration (km/s^2).
            - total_aero_torque (numpy.ndarray): Total aerodynamic torque 
                                                   (Nm).
    """

    ###################
    #get velocity relative 
    #to rotating amtosphere
    ###################

    v = v_eci *1000                                     #in m/s
    vec_v_rel = v - np.cross(omega_E, r_eci*1000)       # Include rotation of atmosphere
    v_rel = np.linalg.norm(vec_v_rel)                   #norm
    v_dir = vec_v_rel / v_rel                           #velocity unit vector



    ###################
    #Lift direction
    ###################

    #get face normal in inertial frame
    arr_normal_in = body_to_inertial(arr_normal, quaternion)

    #get lift direction for each facet
    arr_v_dir_lift = np.cross(np.cross(v_dir, arr_normal_in), v_dir)


    ###################
    #facet normals
    #relative to velocity 
    #direction
    ###################

    #need to normalise all the normals relative to their own magnitude

    #define inverse of norm of each row
    arr_norm_normal_in = 1 / np.linalg.norm(arr_normal_in, axis=1) 

    #multiply each normal vector by its inverse norm: normalisation
    arr_normal_in = np.einsum('i, ij->ij', arr_norm_normal_in, arr_normal_in) 

    #angle of normal relative to velocity direction
    dott = np.dot(arr_normal_in, v_dir)
    arr_angle = np.arccos(dott)


    ###################
    #Aerodynamic 
    #acceleration
    ###################

    #atmospheric density
    h = (np.linalg.norm(r_eci) - R_E)                   #altitude
    rho = density(h, Ta)                                #density
    


    #compute drag and lift coefficient, scaled for area of facet, Aref=1
    arr_Cd, arr_Cl = MoeSentman(arr_angle, arr_area, 1, Ta, arr_Tw, v_rel, arr_alpha_coeff )
    


    #drg and lift forces --#area not present as Cd /Cl alerady scaled for it.
    drag = -  (0.5 * (arr_Cd) * rho * (v_rel**2)) / m  #SI
    lift =  (0.5 * (arr_Cl) * rho * (v_rel**2)) / m # m/s


    #in inertial frame
    arr_vec_drag = np.einsum('i, ij->ij', drag, np.array([v_dir]))      #drag is same dir for all facets
    arr_vec_lift = np.einsum('i, ij->ij', lift, arr_v_dir_lift)         #lift is normal to drag




    ###################
    #sum of the vectors
    #-- total aerodynamic
    # acceleration on each facet
    ###################

    arr_vec_a = arr_vec_drag + arr_vec_lift

    ##Compute total acceleration and torque acting on satellite as a whole
    # torques: force on each facet times the position of each facet

    #get force in SI units
    arr_aero_force = m* arr_vec_a

    #rotate inertial forces into body frame
    body_frame_aero_force = inertial_to_body(arr_aero_force, quaternion) 
    
    #cross force with facet CoM position to obtain torque
    arr_aero_torque = np.cross(arr_position, body_frame_aero_force)
    
    #sums all the forces to get the total force /torque
    total_aero_acceleration = np.einsum('ij->j', arr_vec_a/1000)    # acceleration in km/s^2
    total_aero_torque = np.einsum('ij->j', arr_aero_torque)         #torque in Nm      


    return total_aero_acceleration, total_aero_torque  #acceleration & torque due to drag 






def solarRadiationPressure( vec_r, vec_sun, m, quaternion, arr_area, arr_normal, arr_position, arr_Cr, lightSpeed, centralBodyRadius, starRadius):
    """
    Calculates the total solar radiation pressure (SRP) acceleration and torque acting on a spacecraft.

    This function implements a model to determine the SRP force acting on 
    multiple facets (surfaces) of a spacecraft. 

    Args:
        vec_r (numpy.ndarray): Spacecraft position vector in ECI frame (km).
        vec_sun (numpy.ndarray): Sun position vector in ECI frame (km).
        m (float): Spacecraft mass (kg).
        quaternion (numpy.ndarray): Spacecraft attitude quaternion.
        arr_area (numpy.ndarray): Array of areas (m^2) for each spacecraft facet.
        arr_normal (numpy.ndarray): Array of unit normal vectors in body frame 
                                    for each facet.
        arr_position (numpy.ndarray): Array of positions (meters) of each facet's 
                                      center of mass relative to the spacecraft's 
                                      center of mass.
        arr_Cr (numpy.ndarray): Array of solar radiation coefficients for each facet.

    Returns:
        tuple:
            - SRP_perturb (numpy.ndarray): Total SRP acceleration (km/s^2).
            - SRP_torque (numpy.ndarray): Total SRP torque (Nm).
    """

    #constants:


    #Percent of sun disk visible
    p = percentShadow(vec_r, vec_sun,centralBodyRadius, starRadius)               #in/out of shadow


    #Useful vectors
    sun_dir = vec_sun / np.linalg.norm(vec_sun)     # unit vector of sun 
    vec_sun_sat = vec_r - vec_sun                   #sun to satellite vector, km 


    ##solar pressure
    solarIrradiance = 3.844405613e26 / ((4 * np.pi) * (np.linalg.norm(vec_sun_sat*1000)**2))    #W/m^2         
    SolarPressure = solarIrradiance / lightSpeed                                                #N/m^2


    #get facet normal in the inertial frame
    arr_normal_in = body_to_inertial(arr_normal, quaternion)




    #compute if facet is facing sun direction
    dot_prod = np.dot(arr_normal_in, sun_dir)               #dot product of 2 unit vectors
    arr_facet_state = dot_prod > 0                          #boolean array
    #projected area
    arr_A = arr_area * dot_prod  * arr_facet_state          #null if facet is not facing the sun direction.




    #compute solar pressure acceleration

    #magnitude of force for every facet
    arr_srp = (p/100) * SolarPressure *(( arr_Cr * arr_A) / m)      #SolarPressure
    #force direction (unit vector)
    vec_srp_dir = (vec_sun_sat / np.linalg.norm(vec_sun_sat))

    #array of acceleration vectors (arr_srp[i] * vec_srp_dir)
    arr_vec_a = np.einsum('i, ij->ij', arr_srp, np.array([vec_srp_dir])) 



    #Compute acceleration and torque acting on whole satellite
    # torques: force on each facet times the position of each facet
    #get force in SI units:
    arr_SRP_force = m * arr_vec_a

    #rotate force in body-frame
    body_frame_SRP_force = inertial_to_body(arr_SRP_force, quaternion) 
    #compute body-fixed torques
    arr_SRP_torque = np.cross(arr_position, body_frame_SRP_force)
    
    #sums all the forces to get the total force /torque
    SRP_perturb = np.einsum('ij->j', arr_vec_a/1000)    # total acceleration---km/s
    SRP_torque = np.einsum('ij->j', arr_SRP_torque)     #total torque, Nm 

    #output is an array of force vectors
    return SRP_perturb, SRP_torque #total acceleration and torque from SRP effects

