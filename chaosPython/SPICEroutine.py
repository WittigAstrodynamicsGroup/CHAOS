'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


SPICEroutine.py



This script provides functions that interact with the NASA Ames Information System (NAIF) SPICE software library. 
SPICE is used for planetary ephemeris calculations and spacecraft mission planning. 

The script defines functions for:

* Converting dates to SPICE time format
* Retrieving positions of celestial bodies (Sun, Moon, etc.) using SPICE
* Applying and removing Earth's rotation to vectors using SPICE
* Loading and unloading SPICE kernel files

**NOTE:** This script assumes the SPICE kernels are loaded before using these functions. 
'''

import numpy as np
import os
import spiceypy




def strMonth(month):
    """
    Converts a month number to its corresponding abbreviation (e.g., 1 becomes 'JAN').
    This is used for compatibility with SPICE date and name format
    Args:
        month (int): Month number (1-12)

    Returns:
        str: Abbreviation of the month name (e.g., 'JAN', 'FEB', etc.)
    """

    monthDict = {1:'JAN', 2:'FEB', 3:'MAR', 4:'APR', 5:'MAY', 6:'JUN', 7:'JUL', 8:'AUG', 9:'SEP', 10:'OCT', 11:'NOV', 12:'DEC'}

    return monthDict[month]


def spiceSunPos(t, epoch, str_refObj):
    """
    Retrieves the Sun's position in kilometers from a reference object using SPICE.

    This function takes the time since the simulation started (t), the start epoch date (epoch), 
    and the reference object name (str_refObj) as inputs. It uses SPICE to calculate 
    the Sun's position relative to the reference object in the J2000 inertial frame at the 
    specified time. The position is returned in kilometers.

    Args:
        t (float): Time since simulation epoch in seconds
        epoch (tuple): Epoch date as a tuple of (day, month, year)
        str_refObj (str): Reference object name (e.g., 'Earth')

    Returns:
        numpy.ndarray: Sun's position vector relative to the reference object (km)
    """

    #convert the epoch to TDB [Temps Dynamique Barycentrique]
    day, month, year = epoch
    str_month = strMonth(month)
    str_time = '{} {} {} 12:00:00'.format(year, str_month, day)


    et = spiceypy.str2et(str_time)                  #number of TDB seconds since J2000 until epoch
    et+= t                                          #add time that has past since simulation


    vec_r_sun, corr = spiceypy.spkpos('sun', et,
                                       'J2000',
                                        'none',
                                        str_refObj) #position of sun wrt an
                                                                                #object in the J2000Eq frame.


    return vec_r_sun #km


def spiceMoonPos(t, epoch, str_refObj):
    """
    Retrieves the Moon's position in kilometers from a reference object using SPICE.

    This function is similar to `spiceSunPos` but retrieves the Moon's position instead.

    Args:
        t (float): Time since simulation epoch in seconds
        epoch (tuple): Epoch date as a tuple of (day, month, year)
        str_refObj (str): Reference object name (e.g., 'Earth')

    Returns:
        numpy.ndarray: Moon's position vector relative to the reference object (km)
    """
    
    #convert the epoch to TDB [Temps Dynamique Barycentrique]
    day, month, year = epoch
    str_month = strMonth(month)
    str_time = '{} {} {} 12:00:00'.format(year, str_month, day)

    
    et = spiceypy.str2et(str_time) #number of TDB seconds since J2000 until epoch
    et+= t #add time that has past since simulation



    vec_r_moon, corr = spiceypy.spkpos('moon',
                                        et,
                                        'J2000',
                                        'none',
                                        str_refObj)         #position of sun wrt an object in the J2000Eq frame.



    return vec_r_moon #km

def spiceTargPos(t, epoch, str_targObj, str_refObj):
    """
    Retrieves the position of a target celestial body in kilometers from a reference object using SPICE.

    This function allows retrieving the position of any celestial body specified by its name 
    (str_targObj) relative to a reference object (str_refObj) at a specific time.

    Args:
        t (float): Time since simulation epoch in seconds
        epoch (tuple): Epoch date as a tuple of (day, month, year)
        str_targObj (str): Target celestial body name (e.g., 'Mars', 'ISS')
        str_refObj (str): Reference object name (e.g., 'Earth')

    Returns:
        numpy.ndarray: Target body's position vector relative to the reference object (km)
    """
    
    #convert the epoch to TDB [Temps Dynamique Barycentrique]
    day, month, year = epoch
    str_month = strMonth(month)
    str_time = '{} {} {} 00:00:00'.format(year, str_month, day)

    
    et = spiceypy.str2et(str_time)                  #number of TDB seconds since J2000 until epoch
    et+= t                                          #add time that has past since simulation
    vec_r_targ, corr = spiceypy.spkpos(str_targObj,
                                        et,
                                        'J2000',
                                        'none',
                                        str_refObj) #position of sun wrt an object in the J2000Eq frame.



    return vec_r_targ #km

def applyEarthRotation(t, epoch, r_vec):
    """
    Applies Earth's rotation to a vector using SPICE.

    This function takes a vector (r_vec) and applies Earth's rotation at a specific time 
    (defined by epoch and t) to the vector. The rotation is modeled by transforming the 
    vector from the J2000 inertial frame to the ITRF93 reference frame, which accounts for 
    Earth's rotation.

    Args:
        t (float): Time since simulation epoch in seconds
        epoch (tuple): Epoch date as a tuple of (day, month, year)
        r_vec (numpy.ndarray): Input vector (in J2000 frame)

    Returns:
        tuple: A tuple containing the rotated vector (r_rot) and the rotation matrix
    """

    #convert time to a format usable by SPICE
    day, month, year = epoch
    str_month = strMonth(month)
    str_time = '{} {} {} 12:00:00'.format(year, str_month, day)

    et = spiceypy.str2et(str_time)                  #number of TDB seconds since J2000 until epoch
    et+= t                                          #add time that has past since simulation

    #rotation matrix between inertial and Earth-fixed frame (itrf) at time "et"
    rot_matrix = spiceypy.pxform('J2000', 'itrf93', et)

    #rotate the vector
    r_rot = np.matmul(rot_matrix, r_vec)


    return r_rot, rot_matrix

def undoEarthRotation(r_rot, rot_matrix):
    """
    Transforms a vector in the Earth-fixed frame to the inertial frame.

    Args:
        r_rot (numpy.ndarray): Vector that is in the ITRF frame
        rot_matrix (numpy.ndarray): The rotation matrix that was used to apply the rotation

    Returns:
        numpy.ndarray: The original vector with Earth's rotation removed
    """

    matrix = np.transpose(rot_matrix)

    r_vec = np.matmul(matrix, r_rot)

    return r_vec




def loadKernels(kernels_to_load):
    """
    Loads specified SPICE kernels from the local "Data/kernels" directory.

    Args:
        kernels_to_load (list): List of kernel filenames to load
    """

    ##Looks for the kernels in the local "Data/kernels" folder
    pathToKernels = os.path.dirname(__file__)
    for kernel in kernels_to_load:
        spiceypy.furnsh(pathToKernels + '/Data/kernels/' + kernel)
    return None


def unloadKernels(kernels_to_load):
    """
    Unloads specified SPICE kernels from the local "Data/kernels" directory.

    Args:
        kernels_to_load (list): List of kernel filenames to unload
    """

    ##Looks for the kernels in the local "Data/kernels" folder
    pathToKernels = os.path.dirname(__file__)
    for kernel in kernels_to_load:
        spiceypy.unload(pathToKernels + '/Data/kernels/' + kernel)
    return None
