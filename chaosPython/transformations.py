# Python code for useful transformation functions

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


transformations.py


This Python code defines functions for performing various transformations, including:
    - Coordinate transformations (e.g., spherical to cartesian)
    - Conversions between element representations (e.g., keplerian to equinoctial)
    - Rotations between reference frames (e.g., inertial to body-fixed)
'''



import warnings
import numpy as np

# mu = 3.986004418e5
# VEC_EARTH_ROTATION = np.array([0, 0, 0.0000729211585530])
# R_E = 6378.1363 # km//
# R_E_polar = R_E#6356.751 #km, polar radius
# ecc_R_E = 0#0.081819221456
'''
Transformation of r,v vectors to keplerian elements
'''
def rv2kpl(VEC_r, VEC_v):
    '''
    This function transforms position (VEC_r) and velocity (VEC_v) vectors from ECI (Earth-centered inertial) frame
    to keplerian elements. Algorithm from Vallado.

    Args:
        VEC_r (numpy.ndarray): Position vector in ECI frame (km)
        VEC_v (numpy.ndarray): Velocity vector in ECI frame (km/s)

    Returns:
        numpy.ndarray: A numpy array containing the keplerian elements (a, e, i, RAAN, wp, theta)
    '''


    #useful constants:
    mu = 3.986004418e5              #Earth gravitational constant


    r = np.linalg.norm(VEC_r)       #Norm of position
    v = np.linalg.norm(VEC_v)       #Norm of velocity
    VEC_K = np.array([0, 0, 1])     #Unit vector along z direction


    #find angular momentum
    VEC_h = np.cross(VEC_r, VEC_v)
    h = np.linalg.norm(VEC_h)



    #find the vector pointing to the asc node
    VEC_n = np.cross(VEC_K, VEC_h)
    n = np.linalg.norm(VEC_n)



    #find eccentricity vector:
    E1 = v**2 - (mu/r)
    E2 = np.dot(VEC_r, VEC_v)
    # print('E2', E2)
    VEC_e = (E1 * VEC_r - np.dot(E2, VEC_v)) / mu
    e = np.linalg.norm(VEC_e)



    #find specific mechanical energy
    zeta = ( (v**2) / (2) ) - (mu/r)

    if e != 1.0:
        a = -(mu)/(2 * zeta)
        p = a * (1 - e**2)
    else:
        p = (h**2) / mu
        a = np.inf
    


    #find inclination
    inc = np.arccos(VEC_h[2] / h)
    warnings.warn('The inclination has been found to be sometimes wrong. Quadrant disambiguation?')



    #find ascending node
    RAAN = np.arccos(VEC_n[0] / n)
    if VEC_n[1] < 0 :
        RAAN = 2 * np.pi - RAAN
    


    #find argument of perigee
    W1 = np.dot(VEC_n, VEC_e)
    W2 = n * e
    wp = np.arccos(W1 / W2)
    if VEC_e[2] < 0:
        wp = np.pi * 2 - wp
    


    #find true anomaly
    T1 = np.dot(VEC_e, VEC_r)
    T2 = e * r
    theta = np.arccos(T1 / T2)
    if np.dot(VEC_r, VEC_v) < 0:
        theta = 2 * np.pi - theta



    #if elliptic equatorial:
    if e > 0 and inc < 1e-10:

        theta = np.arccos(VEC_e[0] / e)
        if VEC_e[1] < 0 :
            theta = 2*np.pi - theta
    
    #if circular inclined:
    if e < 1e-10 and inc > 0:

        dots = np.dot(VEC_n, VEC_r) / (n * r)
        theta = np.arccos(dots)
        if VEC_r[2] < 0 :
            theta = 2*np.pi - theta
    
    #if circular equatorial
    if e < 1e-10 and inc < 1e-10:

        theta = np.arccos(VEC_r[0] / r)
        if VEC_r[1] <0:
            theta = 2*np.pi -theta


    #return keplerian elements
    return np.array([a, e, inc, RAAN, wp, theta])


def perifocal_vel(theta, h, e):
    '''
    This function calculates the perifocal velocity vector from keplerian elements .

    Args:
        theta (float): True anomaly (rad)
        h (float): Angular momentum (km^2/s)
        e (float): Eccentricity

    Returns:
        numpy.ndarray: Perifocal velocity vector (km/s)
    '''
    #function returns the perifocal velocity vector at position TA


    #useful constants:
    mu = 3.986004418e5

    const  = mu / h
    velocity = const * np.array([-np.sin(theta), (e + np.cos(theta)), 0])
    return velocity


def perifocal_pos(theta, a, e):
    '''
    This function calculates the perifocal position vector from keplerian elements .

    Args:
        theta (float): True anomaly (rad)
        a (float): semi-major axis (km)
        e (float): Eccentricity

    Returns:
        numpy.ndarray: Perifocal position vector (km)
    '''
    #useful constants:
    mu = 3.986004418e5

    #function to return the position vector at TA in perifocal frame
    r = (a*(1-e**2)) / (1+e*np.cos(theta)) 
    pos =  np.array([r *np.cos(theta), r *np.sin(theta), 0])
    return pos

def peri_cartesian_rotation(angles):
    '''
    This function calculates the rotation matrix to transform vectors from the perifocal frame to the cartesian frame.

    Args:
        angles (list): List of angles (in radians): RAAN, argument of perigee (wp), inclination (INC), true anomaly (theta)

    Returns:
        numpy.ndarray: Rotation matrix (3x3)
    '''

    RAAN, wp, INC, theta = angles


    #function returns the rotation matrix to transform  r & v in perifocal into cartesian coordinates
    R11 = (np.cos(RAAN) * np.cos(wp) - np.cos(INC) * np.sin(RAAN) * np.sin(wp))
    R21 = (np.cos(wp) * np.sin(RAAN) + np.cos(RAAN) * np.cos(INC) * np.sin(wp))
    R31 = (np.sin(INC) * np.sin(wp))

    R12 = (-np.cos(RAAN) * np.sin(wp) - np.cos(INC) * np.cos(wp) * np.sin(RAAN))
    R22 = (np.cos(RAAN) * np.cos(INC) * np.cos(wp) - np.sin(RAAN) * np.sin(wp))
    R32 = (np.cos(wp) * np.sin(INC))

    R13 = (np.sin(RAAN) * np.sin(INC))
    R23 = (-np.cos(RAAN) * np.sin(INC))
    R33 = (np.cos(INC))
    


    R = np.array([[R11, R12, R13], [R21, R22, R23], [R31, R32, R33]])  #full rotation matrix


    return R

def cartesian_rsw_rotation(angles):
    '''
    This function calculates the rotation matrix to transform vectors from the cartesian frame to the RSW frame (Radial, s-plane, normal to orbit plane).
    Needed for Gaussian Osculating element computation. From Curtis
    Args:
        angles (list): List of angles (in radians): RAAN, argument of perigee (wp), inclination (i), true anomaly (theta)

    Returns:
        numpy.ndarray: Rotation matrix (3x3)
    '''
    RAAN, wp, i, theta = angles
    u = wp + theta

    R = np.zeros((3,3))


    R[0, 0] = -np.sin(RAAN) * np.cos(i) * np.sin(u) + np.cos(RAAN) * np.cos(u)
    R[0, 1] = np.cos(RAAN) * np.cos(i) * np.sin(u) + np.sin(RAAN) * np.cos(u)
    R[0, 2] =  np.sin(i) * np.sin(u)
    
    R[1, 0] = -np.sin(RAAN) * np.cos(i) * np.cos(u) - np.cos(RAAN) * np.sin(u) 
    R[1, 1] = np.cos(RAAN) * np.cos(i) * np.cos(u) - np.sin(RAAN) * np.sin(u)
    R[1, 2] = np.sin(i) * np.cos(u)

    R[2, 0] = np.sin(RAAN) * np.sin(i)
    R[2, 1] = -np.cos(RAAN) * np.sin(i)
    R[2, 2] = np.cos(i)
    return R

def gregorian_to_J0(day, month, year):
    '''
        Function to transform gregorian date ( day, month, year) into Julian Date J0. Julian date starts in antiquity at noon. 
        This function does not include the time along the day. (i.e. no difference between 2 pm and 6 pm in this function)
        From Vallado and Curtis
    '''

    A = (month + 9) / 12
    J0 = 367 * year - int( ( 7 * (year + int(A)) ) / 4 ) + int((275 * month) / 9) + day + 1721013.5

    return J0

def julianDate(UT, day, month, year):
    '''
    Function to calculate the julian date from gregorian calendar date (Day , month, year) and UT in fractional hours. 
    Universal Time is measured from the time of the Sun's
    passage over the Greenwich meridian. The day changes at noon in JD!!
    Fractional hours mean using the following operation: hours + (min / 60) + (seconds / 3600)
    '''

    J0 = gregorian_to_J0(day, month, year)
    JD = J0 + UT/24

    return JD #JD is days!!

def ModJulianDate(UT, day, month, year):
    '''
    Function to calculate the modifed Julian date from gregorian calendar date and UT in fractional hours. The day changes at midnight.
    Modified Julian date starts on the 23rd of May 1968
    
    For the GMAT modified Julian date, GMAT uses:MJD = JD - 2430000 // 5 Jan 1941
    '''
    
    JD = julianDate(UT, day, month, year)
    MJD = JD - 2400000.5   


    return MJD #in days, civilian 


def kep_to_eq(a, e, i, RAAN, wp, theta):
    """
    This function transforms keplerian elements to equinoctial elements.

    Args:
        a (float): Semi-major axis (km)
        e (float): Eccentricity
        i (float): Inclination (rad)
        RAAN (float): Right ascension of the ascending node (rad)
        wp (float): Argument of perigee (rad)
        theta (float): True anomaly (rad)

    Returns:
        numpy.ndarray: A numpy array containing the equinoctial elements (p, f, g, h, k, L)
    """

    p = a * (1 - e**2)      # semi-parameter

    f = e * np.cos(wp + RAAN)
    g = e * np.sin(wp + RAAN)
    h = np.tan(i/2) * np.cos(RAAN)
    k = np.tan(i/2) * np.sin(RAAN)
    L = RAAN + wp + theta

    return np.array([p, f, g, h, k, L])

def eq_to_kep(p, f, g, h, k, L):
    """
    This function transforms from equinoctial elements to keplerian elements.

    Args:
        p (float): Semi-parameter
        f (float): f element
        g (float): g element
        h (float): h element
        k (float): k element
        L (float): Mean longitude

    Returns:
    numpy.ndarray: A numpy array containing the keplerian elements (a, e, i, RAAN, wp, theta)
    """

    a = p / ((1 - (f**2) - (g**2)))
    e = np.sqrt((f**2) + (g**2))
    i = np.arctan2(2 * np.sqrt((h**2) + (k**2)), 1-(h**2)-(k**2))
    wp = np.arctan2(g*h - f*k, f*h + g*k)
    RAAN = np.arctan2(k,h)
    theta = L - (RAAN + wp)

    return np.array([a, e, i, RAAN, wp, theta])




def equi_to_r_cart(p, f, g, h, k, L):
    """
    This function transforms equinoctial elements to ECI position vector (r).

    Args:
        p (float): Semi-parameter
        f (float): f element
        g (float): g element
        h (float): h element
        k (float): k element
        L (float): Mean longitude

    Returns:
        numpy.ndarray: ECI position vector (km)
    """

    #compute some intermediate values
    w = 1 + f * np.cos(L) + g * np.sin(L)
    s = 1 + (h**2) + (k**2)
    alpha = (h**2) - (k**2)
    r_e = p / w
    
    #get ECI state vector
    r_ECI = np.zeros(3)
    r_ECI[0] =  (np.cos(L) + (alpha) * np.cos(L) + 2*h*k*np.sin(L)) 
    r_ECI[1] =  (np.sin(L) - (alpha) * np.sin(L) + 2*h*k*np.cos(L))
    r_ECI[2] = 2  * (h*np.sin(L) - k * np.cos(L))
    r_ECI = (r_ECI * (r_e / (s)))
    
    return r_ECI

def equi_to_v_cart(modEquinoctialElements):
    """
    This function transforms equinoctial elements to ECI velocity vector (v).

    Args:
        modEquinoctialElements (list): List containing equinoctial elements (p, f, g, h, k, L)

    Returns:
        numpy.ndarray: ECI velocity vector (km/s)
    """



    p, f, g, h, k, L = modEquinoctialElements           #Unpack the list

    #useful constants:
    mu = 3.986004418e5              #Earth gravitational parameter

    s = 1 + (h**2) + (k**2)
    alpha = (h**2) - (k**2)

    v_ECI = np.zeros(3)
    v_ECI[0] = -(np.sin(L) + (alpha) * np.sin(L) - 2 * h * k * np.cos(L) + g - 2 * f* h * k + (alpha) * g)
    v_ECI[1] = -(-np.cos(L) + (alpha) * np.cos(L) + 2 * h * k * np.sin(L) - f + 2 * g * h * k + (alpha) * f)
    v_ECI[2] = 2 * (h*np.cos(L) + k*np.sin(L) + f*h + g*k)

    v_ECI = v_ECI * np.sqrt(mu/p)/(s)

    return v_ECI



def rotation(r, angle):
    """
    NOTE: DEPRECATED FUNCTION
    This function creates a rotation matrix to transform a vector from the inertial frame to the body-fixed frame.
    This function looks at a pure rotation around the z-axis
    Args:
        r (numpy.ndarray): Input vector in the inertial frame (km)
        angle (float): Rotation angle of the frame in radians

    Returns:
        numpy.ndarray: Rotated vector in the body-fixed frame (km)

    Notes:
        - Rotation is performed around the z-axis.
        - The rotation matrix transforms the reference frame, not the vector itself.
    """
    
    R = np.zeros((3,3))
    R[0,0] = np.cos(angle)
    R[0,1] = np.sin(angle)
    R[1,0] = -np.sin(angle)
    R[1,1] = np.cos(angle)
    R[2,2] = 1

    r_rot = np.matmul(R, r)

    return(r_rot) 

def v_rotation(v_inertial, r_rot):
    """
    This function calculates the velocity vector in the rotating body-fixed frame from the inertial velocity vector.

    Args:
        v_inertial (numpy.ndarray): Inertial velocity vector (km/s)
        r_rot (numpy.ndarray): Rotated position vector in the body-fixed frame (km)

    Returns:
        numpy.ndarray: Velocity vector in the rotating body-fixed frame (km/s)
    """

    #Constant:
    VEC_EARTH_ROTATION = np.array([0, 0, 0.0000729211585530])


    v_rot = (v_inertial - np.cross(VEC_EARTH_ROTATION, r_rot))
    
    return v_rot 

def a_rotation(a_itrf, r_rot, v_rot):
    """
    This function calculates the acceleration vector in the inertial frame from the body-fixed acceleration vector.

    Args:
        a_itrf (numpy.ndarray): Acceleration vector in the ITRF frame (Earth-Fixed) (km/s^2)
        r_rot (numpy.ndarray): Rotated position vector in the body-fixed frame (km)
        v_rot (numpy.ndarray): Velocity vector in the rotating body-fixed frame (km/s)

    Returns:
        numpy.ndarray: Acceleration vector in the inertial frame (km/s^2)
    """
    
    VEC_EARTH_ROTATION = np.array([0, 0, 0.0000729211585530])
    A = np.cross(VEC_EARTH_ROTATION, r_rot)
    
    B = np.cross(VEC_EARTH_ROTATION, A) #  centrifugal acceleration
    C = 2*np.cross(VEC_EARTH_ROTATION, v_rot) #coriolis acceleration
    
    a_eci = a_itrf + B + C
    return a_eci #output in km


def LatLong(r_rot):
    """
    This function transforms a position vector in the Earth-fixed frame 
    to geocentric latitude, longitude, and radius.
    Algorithm from Vallado.

    Args:
        r_rot (numpy.ndarray): Input position vector in Earth-fixed frame (km)

    Returns:
        tuple: A tuple containing the altitude (h), geocentric latitude (phi_gc), and longitude (lambda) in km and radians.
    """

    ##some constant -- For simplicity, we assume the Earth is a perfect sphere.
    R_E = 6378.1363 # km
    R_E_polar = R_E                 #real polar radius: 6356.751 #km, polar radius
    ecc_R_E = 0                     # Earth shape ecc: 0.081819221456
    tol=1e-8
    r_ecef = r_rot

    r_dsat = np.sqrt( (r_ecef[0]**2) + (r_ecef[1]**2) )


    lmbda = np.arctan2(r_ecef[1], r_ecef[0])
    delta = np.arctan(r_ecef[2]/r_dsat)


    guess = 0
    new = delta
    while abs(new - guess) > tol:
        guess = new
        C = R_E / (np.sqrt(1 - (ecc_R_E**2) * np.sin(guess)*np.sin(guess)))
        new = np.arctan((r_ecef[2] + C * (ecc_R_E**2)*np.sin(guess)) / r_dsat)

    h = np.linalg.norm(r_ecef)
    phi_gd = new

    phi_gc = np.arctan((1-ecc_R_E**2) * np.tan(phi_gd))

    return h, phi_gc, lmbda



def spherical_to_r_cart(radius, phi, lmbda):
    """
    This function transforms spherical coordinates (radius, latitude, longitude) to a position vector in the Cartesian frame.

    Args:
        radius (float): Radius (km)
        phi (float): Latitude (radians)
        lmbda (float): Longitude (radians)

    Returns:
        numpy.ndarray: Position vector in Cartesian frame (km)
    """
    position = radius * np.array([np.cos(phi)*np.cos(lmbda), np.cos(phi)*np.sin(lmbda), np.sin(phi)])
    
    return position





def kpl2rv(a, e, inc, raan, wp, ta, mu):
    """
    This function transforms Keplerian elements to position (r) and velocity (v) vectors in the ECI frame.
    Algorithm from Vallado
    Args:
        a (float): Semi-major axis (km)
        e (float): Eccentricity
        inc (float): Inclination (rad)
        raan (float): Right ascension of the ascending node (rad)
        wp (float): Argument of perigee (rad)
        ta (float): True anomaly (rad)
        mu (float): Gravitational parameter


    Returns:
        tuple: A tuple containing the position vector (numpy.ndarray) and velocity vector (numpy.ndarray) in ECI frame (km, km/s).
    """



    u = ta + wp
    w_true = raan + wp

    if e== 0 and inc == 0:
        wp = 0
        raan = 0
        

    elif e ==0  and inc !=0:
        wp = 0
        ta = u

    elif e!=0 and inc ==0:

        raan = 0
        wp = w_true

        
    p = a * (1 - e**2)
    r = p / (1 + e* np.cos(ta))
    const = np.sqrt(mu / p)


    r_peri = np.array([r * np.cos(ta), r * np.sin(ta), 0])

    v_peri = np.array([-const * np.sin(ta), const * (e + np.cos(ta)), 0])

    angles = np.array([raan, wp, inc, ta])
    rot_matrix = peri_cartesian_rotation(angles)

    vec_r_cart = np.matmul(rot_matrix, r_peri)
    
    vec_v_cart = np.matmul(rot_matrix, v_peri)

    return vec_r_cart, vec_v_cart





