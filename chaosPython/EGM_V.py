#python code to implement the non-Keplerian gravity model --- Vallado method

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


EGM.py



Python code to generate the non-keplerian gravity model, based on the EGM08 data and the Vallado approach. Multiple functions are defined:
-   matricise():    which extract the values of C and S, de-normalises them and store them in respective matrices
-   ALP():          which generates non-normalised Associated Legendre Polynomials through recursion as stated Vallado
-   dU_d[item]:     which calculates the derivatives of the potential with respect to [item]
-   static_egm08(): which calculates the acceleration due to the non-keplerian gravity in the body-fixed frame (Earth rotation not included)
-   egm08():        which calculates the acceleration due to the non-keplerian gravity in the inertial frame (Earth rotation included)
'''



import pandas as pd
import numpy as np
import spiceypy
from scipy.special import factorial
from scipy.special import lpmn
from .transformations import julianDate, spherical_to_r_cart
from .SPICEroutine import *



def matricise(location, filename = 'EGM96_n36', max_degree=6):

    '''
    Function to extract the C and S coefficient from the EGM file and store them in respective matrices. 
    The coefficients are de-normalised.
    loc = 2: C coeff 
    loc = 3: S coeff
    '''
    
    dataframe = pd.read_fwf(filename)
    
    
    # +1 added as buffer to allow loop to run with an offset (Because we start the matrix at coefficients C of degree 2)
    C = np.empty((max_degree +1 , max_degree +1 )) 
    C[:] = 0  #set all the indexes as 0 in the matrix
    
    '''
    Define the loop for the pyramid/ matrix
    '''

    loc = 0  # set variable to iterate through textfile

    for n in range(2, max_degree + 1): # -1):   #this loop points towards a row at each iteration
        m = 0  #set variable to iterate through matrix row
        
        while m <= (n):# + start_degree) :   #create offset to have 3 data point in the first row // this loop points at each index within a row of the matrix
            
            if m == 0 :
                k = 1
            else:
                k = 2

            nominator =  factorial(n + m)
            denominator = (factorial(n-m)) * k * (2 * n + 1) 
            scalor = np.sqrt(nominator / denominator)  #de-scale the coefficients for ease s
            
            # print(n, m, scalor)
            C[n][m] =(float(dataframe.iloc[loc][location]))#.replace('D', 'e')))
            C[n][m] = C[n][m] / scalor
            m += 1  
            loc += 1
    
    C[0:2] = 0   #buffer so that the matrix row / column corresponds to the degree/order of the coeff
    
    if location == 2:
        C[0][0] = 1
    return C



def ALP(phi, degree):
    '''
    Function to create a matrix of unitless, non-normalised Associated Legendre Polynomials based on the recursion from Vallado
    NOTE: phi in rad
    '''
    degree += 1
    order = degree
    grid = np.empty((degree , order))
    grid[:] = 0

    #initial values 
    grid[0][0] = 1
    grid[1][0] = np.sin(phi)
    grid[1][1] = np.cos(phi)

    for n in range(degree - 2 ):
        n = n + 2 #degree between 2 and max degree
        m = 0   # order, restarted at 0 for each degree

        while m < n:
            if m == 0 :
                k = 1
            else:
                k = 2


            
            
            grid[n][m] = grid[n-2][m] + ( 2* n - 1) * np.cos(phi) * grid[n-1][m-1] #mid-values of each row
            m += 1 
        


        grid[n][0] =  ( (2 * n - 1) * np.sin(phi) * grid[n-1][0] - (n - 1) * grid[n-2][0] ) / (n) # first value of each row

        grid[n][n] =  (2 * n - 1) * np.cos(phi) * grid[n-1][n-1]  # last value of each row

    return grid

def ALPS(phi, degree, order):
    '''
    Function to create a matrix of unitless, non-normalised Associated Legendre Polynomials based on the recursion from Vallado
    NOTE: phi in rad
    NOTE: this is equivalent to ALP, but faster
    '''
    X = lpmn(order, degree, (np.sin(phi)))[0]

    #I'm proud of how unreadable that is. I feel smart.
    #This is to format the matrix into the same shape as the results of ALP
    ALPS = np.matrix.transpose(np.array([X[i, :]*(-1)**i for i in range(len(X))]))
    return ALPS



#debugging: Verify ALP or ALPS are correct
if __name__ =='__main__' and 1==1:

    it = 4
    phi = 79.3*np.pi/180
    grid = ALP(phi, degree =6)
    
    #Matrix P is provided analytically by Vallado:
    P = np.empty((5, 5))
    P[:] = 0
    P[0][0] = 1
    P[1][0] = np.sin(phi)
    P[1][1] = np.cos(phi)
    P[2][0] = 0.5 * (3 * (np.sin(phi)**(2))  - 1)
    P[2][1] = 3 * np.sin(phi) * np.cos(phi)
    P[2][2] =  3 * (np.cos(phi)**(2))
    P[3][0] = 0.5 * ( 5 * (np.sin(phi))**(3) - 3 * np.sin(phi))
    P[3][1] = 0.5 * np.cos(phi) * (15 * (np.sin(phi))**(2) - 3)
    P[3][2] = 15 * (np.cos(phi))**(2) * np.sin(phi)
    P[3][3] = 15 * (np.cos(phi))**(3)
    P[4][0] = (1 / 8) * (35 * (np.sin(phi))**(4) - 30 * (np.sin(phi))**(2) + 3)
    P[4][1] = (5 / 2) * np.cos(phi) * (7 * (np.sin(phi))**(3) - 3 * np.sin(phi)) 
    P[4][2] = (15 / 2) * (np.cos(phi))**(2) * (7 * (np.sin(phi))**(2) - 1)
    P[4][3] = 105 * (np.cos(phi)**(3)) * np.sin(phi)
    P[4][4] = 105 * (np.cos(phi)**(4))



    print('func', ALP(phi, 2))
    print( '\nvallad/o',P )
    print('ALPS', ALPS(phi, 2, 2))
    # print("vallado recursion", grid[it])
    # print("Vallado analytical", P[it])






       

def dU_dr(r, phi, lmbda, C, S, degree, order):
    """
    This function calculates the derivative of the non-Keplerian potential 
    with respect to the radial distance (r) in the geocentric spherical 
    coordinate system.

    Args:
        r (float): Radial distance from the Earth's center in kilometers (km).
        phi (float): Geocentric latitude in radians (rad).
        lmbda (float): Geocentric longitude in radians (rad).
        C (numpy.ndarray): Array containing the C coefficients of the EGM08 model.
        S (numpy.ndarray): Array containing the S coefficients of the EGM08 model.
        degree (int): Maximum degree of the spherical harmonics used.
        order (int): Maximum order of the spherical harmonics used 
                     (must be less than or equal to degree).

    Returns:
        float: The derivative of the non-Keplerian potential with respect to r (km^2/s^2).
    """
    mu = 398600.4418   #km-5-vallado
    RE = 6378.1363# -- km --EGM08
    degree +=1 

    derivative = 0
    

    P =  ALPS(phi, degree, degree)
    for n in range(2, int(degree) ):
        m = 0   # order (derivative order)
        
        
        while m <= order:
            derivative +=  -(mu / (r**2)) * ((RE / r)**(n) * (n + 1) * P[n][m] * (C[n][m] * np.cos(m * lmbda) + S[n][m] * np.sin(m * lmbda)))
            m += 1

    return derivative
            


def dU_dphi(r, phi, lmbda, C, S,  degree, order):
    """
    This function calculates the derivative of the non-Keplerian potential 
    with respect to the lattitude (phi) in the geocentric spherical 
    coordinate system.

    Args:
        r (float): Radial distance from the Earth's center in kilometers (km).
        phi (float): Geocentric latitude in radians (rad).
        lmbda (float): Geocentric longitude in radians (rad).
        C (numpy.ndarray): Array containing the C coefficients of the EGM08 model.
        S (numpy.ndarray): Array containing the S coefficients of the EGM08 model.
        degree (int): Maximum degree of the spherical harmonics used.
        order (int): Maximum order of the spherical harmonics used 
                     (must be less than or equal to degree).

    Returns:
        float: The derivative of the non-Keplerian potential with respect to phi (km^2/s^2).
    """
    mu = 398600.4418   #km-5-vallado
    RE = 6378.1363# -- km --EGM08


    degree +=1 
    derivative = 0

    P = ALPS(phi,degree, degree)
    for n in range(2, int(degree) ):
        m = 0   # order (derivative order)

        
    
        while m <= order:
            derivative +=  (mu / (r)) * ((RE / r)**(n) * ( P[n][m+1] - m * np.tan(phi) * P[n][m] ) * (C[n][m] * np.cos(m * lmbda) + S[n][m] * np.sin(m * lmbda)))
            m += 1
    return derivative




def dU_dlmbda(r, phi, lmbda,C, S,  degree, order):
    """
    This function calculates the derivative of the non-Keplerian potential 
    with respect to the longitude (lmbda) in the geocentric spherical 
    coordinate system.

    Args:
        r (float): Radial distance from the Earth's center in kilometers (km).
        phi (float): Geocentric latitude in radians (rad).
        lmbda (float): Geocentric longitude in radians (rad).
        C (numpy.ndarray): Array containing the C coefficients of the EGM08 model.
        S (numpy.ndarray): Array containing the S coefficients of the EGM08 model.
        degree (int): Maximum degree of the spherical harmonics used.
        order (int): Maximum order of the spherical harmonics used 
                     (must be less than or equal to degree).

    Returns:
        float: The derivative of the non-Keplerian potential with respect to lmbda (km^2/s^2).
    """

    mu = 398600.4418   #km-5-vallado
    RE = 6378.1363# -- km --EGM08

    degree +=1 
    derivative = 0

    P = ALPS(phi, degree, degree)
    for n in range(2, int(degree) ):  
        m = 0   # order (derivative order)


        while m <= order:
        
            derivative +=  (mu / (r)) * ((RE / r)**(n) * m * P[n][m] * (S[n][m] * np.cos(m * lmbda) - C[n][m] * np.sin(m * lmbda)))
            m += 1
    return derivative




def static_egm08(r, phi, lmbda, C, S, deg, order):
    """
    This function computes the non-Keplerian gravity acceleration at a point 
    specified in geocentric spherical coordinates (r, phi, lambda).

    Args:
        r (float): Radial distance from the Earth's center in kilometers (km).
        phi (float): Geocentric latitude in radians (rad).
        lmbda (float): Geocentric longitude in radians (rad).
        C (numpy.ndarray): Array containing the C coefficients of the EGM08 model.
        S (numpy.ndarray): Array containing the S coefficients of the EGM08 model.
        deg (int): Maximum degree of the spherical harmonics used.
        order (int): Maximum order of the spherical harmonics used 
                     (must be less than or equal to degree).

    Returns:
        numpy.ndarray: The non-Keplerian gravity acceleration in the geocentric 
                       Earth-fixed (ECEF) frame in kilometers per square second squared (km^2/s^2).
    """

    #get cartesian coordinates
    rx, ry, rz = spherical_to_r_cart(r, phi, lmbda)

   #get potential derivatives
    U_r = dU_dr(r, phi, lmbda, C, S, deg, order)
    U_phi = dU_dphi(r, phi, lmbda, C, S, deg, order)
    U_lmbda = dU_dlmbda(r, phi, lmbda, C, S, deg, order) 

    #Compute acceleration in ECEF
    a_i = ((1 / r) * U_r - ( rz / ( (r**2) * np.sqrt((rx**2) + (ry**2)) ) ) * U_phi ) * rx - ( ( 1 / ((rx**2) + (ry**2)) ) * U_lmbda ) * ry #- ((mu * rx) / (r**3)) 

    a_j = ((1 / r) * U_r - ( rz / ( (r**2) * np.sqrt((rx**2) + (ry**2)) ) ) * U_phi ) * ry + ( ( 1 / ((rx**2) + (ry**2)) ) * U_lmbda ) * rx #- ((mu * ry) / (r**3)) 

    a_k = (1 / r) * U_r * rz + ( ( np.sqrt((rx**2) + (ry**2)) ) / (r**2) ) * U_phi #- ((mu * rz) / (r**3)) 
    
    acceleration = np.array([a_i, a_j, a_k])
    
    return acceleration




    



    

def egm08(t, r_cartesian, v_inertial, C, S, epoch, deg, order):
    """
    This function computes the non-Keplerian gravity acceleration acting on a 
    point mass in the inertial frame at a given time. It utilizes the EGM08 
    data and the Vallado model for Earth's gravity field.

    Args:
        t (float): Time in seconds since a reference epoch.
        r_cartesian (numpy.ndarray): Position vector of the point mass in the 
                                      inertial frame (km).
        v_inertial (numpy.ndarray): Velocity vector of the point mass in the 
                                      inertial frame (km/s).
        C (numpy.ndarray): Array containing the C coefficients of the EGM08 model.
        S (numpy.ndarray): Array containing the S coefficients of the EGM08 model.
        epoch (tuple): A tuple containing the year, month, and day of the 
                       reference epoch.
        deg (int): Maximum degree of the spherical harmonics used.
        order (int): Maximum order of the spherical harmonics used 
                     (must be less than or equal to degree).

    Returns:
        numpy.ndarray: The non-Keplerian gravity acceleration acting on the 
                       point mass in the inertial frame (km/s^2).
    """

    #rotate input vector : from inertial to body-fixed rotation(r_cartesian, angle)
    r_rot, rot_matrix = applyEarthRotation(t, epoch, r_cartesian) 
    r, lmbda, phi= spiceypy.reclat(r_rot)
    a_itrf = static_egm08(r, phi, lmbda, C, S, deg, order) #compute body-fixed acceleration (itrf)

    a_eci =  undoEarthRotation(a_itrf, rot_matrix) #transform acceleration into inertial frame

    return a_eci

