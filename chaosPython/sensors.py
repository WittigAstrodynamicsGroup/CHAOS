# Sensor trial script: Testing and validating various sensors

"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


sensors.py


This Python script tests and validates the functionalities of different sensors 
used in the spacecraft environment. The script includes functions for simulating 
sensor measurements and their associated errors.

The script defines functions for:

- Faraday cup current: Calculates the current measured by a Faraday cup sensor 
  based on ion density, spacecraft velocity, collection area, and angle of attack.
- Faraday cup error: Estimates the measurement error in the Faraday cup current 
  due to a specified error fraction.
- Faraday cup angle: Simulates the measured angle of the incoming flow based on 
  the actual angle and the error characteristics of the Faraday cup.
- Effective cup angle: Calculates the effective angle of a sensor considering 
  the signal-to-noise ratio.
- Signal-to-noise ratio: Calculates the signal-to-noise ratio based on ion density, 
  collection area, spacecraft velocity, and error factor.
- MEMS gyroscope: Simulates the measurement of a MEMS gyroscope, including the 
  effects of bias drift and random noise.
    
    The script uses functions from `numpy.random` for generating random numbers.
"""
import numpy as np


def faradayCupCurrent(alpha, v, A, rho_ions, k):
    """
    Calculates the current measured by a Faraday cup sensor.

    Args:
        alpha (float): Angle of attack (radians).
        v (float): Spacecraft velocity (m/s).
        a (float): Collection area of the Faraday cup (m^2).
        rho_ions (float): Ion density (particles/m^3).
        k (float): Constant factor (implementation specific).

    Returns:
        float: Measured current (A).
    """
    q_e = 1.6e-19 #Coulomb
    I = q_e * rho_ions * A * v * np.cos(alpha)
    if I < 0:    
        return 0
    else:
        return I

def faradayCupErrorCst(alpha, k, epsilon, v, A, rho_ions):
    """
    Estimates the measurement error in the Faraday cup current.

    Args:
        alpha (float): Angle of attack (radians).
        k (float): Constant factor (implementation specific).
        epsilon (float): Error fraction.
        v (float): Spacecraft velocity (m/s).
        a (float): Collection area of the Faraday cup (m^2).
        rho_ions (float): Ion density (particles/m^3).

    Returns:
        float: Measurement error (radians).
    """

    q_e = 1.6e-19 #Coulomb

    err = ((1/k) * np.arccos( np.cos(alpha)  + (epsilon / (q_e * A * v * rho_ions)))) - alpha
    return err

def faradayCupAngle(alpha, epsilon, v, A, rho_ions, rng):
    """
    Simulates the measured angle based on the error of the Faraday cup.

    Args:
        alpha (float): Actual angle of attack (radians).
        epsilon (float): Error fraction.
        v (float): Spacecraft velocity (m/s).
        a (float): Collection area of the Faraday cup (m^2).
        rho_ions (float): Ion density (particles/m^3).
        rng (numpy.random.Generator): Random number generator.

    Returns:
        float: Measured angle (radians).
    """
    q_e = 1.6e-19 #Coulomb

    #random value of current error:
    rand = rng.normal(0, epsilon)

    #peak current:
    I_max = (q_e * A * v * rho_ions)
    
    internal = np.cos(alpha)  + ( rand / I_max )
    
    # if current greater than peak, limit to peak
    if abs(internal) > 1:
        internal = 1 * np.sign(internal)

    angle =  np.arccos( internal ) 
    return angle


def effectiveCupAngle(SN, epsilon, rho_ions_total, aperture, v):

    """
    Calculates the effective angle of a sensor considering the signal-to-noise ratio.

    Args:
        signal_to_noise (float): Signal-to-noise ratio.
        epsilon (float): Error fraction.
        rho_ions_total (float): Total ion density (particles/m^3).
        aperture (float): Collection area of the sensor (m^2).
        velocity (float): Spacecraft velocity (m/s).

    Returns:
        float: Effective angle of the sensor (radians).
    """
    q_e = 1.6e-19 #Coulomb

    ratio = q_e / epsilon

    angle = np.arccos(SN /( (ratio) * rho_ions_total * aperture * v) ) 

    return angle #rad


def signalToNoise(rho_ions , A , v, epsilon):
    """
    Calculates the signal-to-noise ratio based on sensor properties.

    Args:
        ion_density (float): Ion density (particles/m^3).
        collection_area (float): Collection area of the sensor (m^2).
        velocity (float): Spacecraft velocity (m/s).
        error_fraction (float): Error fraction.

    Returns:
        float: Signal-to-noise ratio.
    """
    q_e = 1.6e-19

    SN = rho_ions * A * v *q_e / epsilon

    return SN

def MEMSgyro(omega, delta_t, biais_n, sigma_n, sigma_b, rng, rng2):
    """
    Simulates the measurement of a MEMS gyroscope, including bias drift and random noise.

    Args:
        angular_velocity (numpy.ndarray): True angular velocity (rad/s).
        time_delta (float): Sampling time interval (s).
        bias_noise (numpy.ndarray): Initial bias instability (rad/s).
        noise_sigma (float): Standard deviation of random noise (rad/s).
        bias_drift_sigma (float): Standard deviation of bias drift (rad/s).
        random_number_generator1 (numpy.random.Generator): Random number generator for bias drift.
        random_number_generator2 (numpy.random.Generator): Random number generator for noise.

    Returns:
        tuple: (numpy.ndarray, numpy.ndarray) - Updated bias instability and measured angular velocity (rad/s).
    """
    rand_b = rng.normal(0, 1, size=3) 
    rand_n = rng2.normal(0, 1, size=3) 
    biais_new = biais_n + sigma_b * np.sqrt(delta_t) * rand_b
    measurement = omega + 0.5 * (biais_new + biais_n) + np.sqrt(((sigma_n**2)/delta_t) + ((delta_t*(sigma_b**2))/12)) * rand_n
    # print('biais',biais_new, delta_t)
    return biais_new, measurement



def scalar_MEMSgyro(omega, delta_t, biais_n, sigma_n, sigma_b, rng, rng2):
    """
    Simulates the measurement of a scalar MEMS gyroscope, including bias drift and random noise.

    This function assumes a single-axis gyroscope.

    Args:
        angular_velocity (float): True angular velocity (rad/s).
        time_delta (float): Sampling time interval (s).
        bias_noise (float): Initial bias error (rad/s).
        noise_sigma (float): Standard deviation of random noise (rad/s).
        bias_drift_sigma (float): Standard deviation of bias drift (rad/s).
        random_number_generator1 (numpy.random.Generator): Random number generator for bias drift.
        random_number_generator2 (numpy.random.Generator): Random number generator for noise.

    Returns:
        tuple: (float, float) - Updated bias error and measured angular
    """
    
    rand_b = rng.normal(0, 1) 
    rand_n = rng2.normal(0, 1) 
    biais_new = biais_n + sigma_b * np.sqrt(delta_t) * rand_b
    measurement = omega + 0.5 * (biais_new + biais_n) + np.sqrt(((sigma_n**2)/delta_t) + ((delta_t*(sigma_b**2))/12)) * rand_n
    # print(measurement)
    return biais_new, measurement


