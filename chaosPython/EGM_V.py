#python code to implement the non-Keplerian gravity model --- Vallado method

'''
Python code to generate the non-keplerian gravity model, based on the EGM08 data and the Vallado approach. Multiple functions are defined:
-   matricise():    which extract the values of C and S, de-normalises them and store them in respective matrices
-   ALP():          which generates non-normalised Associated Legendre Polynomials through recursion as stated Vallado
-   dU_d[item]:     which calculates the derivatives of the potential with respect to [item]
-   static_egm08(): which calculate the acceleration due to the non-keplerian gravity in the body-fixed frame (Earth rotation not included)

NOTE: debugging code is provided in if statements. internal debugging output is also present
NOTE: For faster run at low degree and order, use the EGM08_n6 file. (valid up to degree n = 6) 
'''



import pandas as pd
import numpy as np
import spiceypy
from scipy.special import factorial
from scipy.special import lpmn
from .transformations import julianDate, spherical_to_r_cart
from .SPICEroutine import *
#extract data into a dataframe 

#set some constants
loc_C = 2
loc_S = 3
start_degree = 2
max_degree = 6
RE = 6378.1363# -- km --EGM08
VEC_EARTH_ROTATION = np.array([0, 0, 0.0000729211585530]) #from Vallado// rad/ s
mu = 398600.4418   #km-5-vallado

'''
Store the value of C and S coefficients in respective matrices
'''
def matricise(location, filename = 'EGM96_n36', max_degree=6):

    '''
    Function to extract the C and S coefficient from the EGM file and store them in respective matrices. 
    The coefficients are de-normalised.
    loc = 2: C coeff , 
    loc = 3: S coeff
    '''
    
    dataframe = pd.read_fwf(filename)
    
    C = np.empty((max_degree +1 , max_degree +1 )) # +1 added as buffer to allow loop to run with an offset (Because we start the matrix at coefficients C of degree 2)
    C[:] = 0  #set all the indexes as NaN in the matrix
    
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
    
    C[0:2] = 0 #np.zeros((2, max_degree + 1))  #buffer so that the matrix row / column corresponds to the degree/order of the coeff
    
    # C = np.concatenate((buffer, C))
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

            # nominator =  factorial(n + m)
            # denominator = ( factorial(n-m)) * (2 * n + 1) * k
            # scalor = np.sqrt(nominator / denominator)
            
            
            grid[n][m] = grid[n-2][m] + ( 2* n - 1) * np.cos(phi) * grid[n-1][m-1] #mid-values of each row
            # grid[n][m] = grid[n][m] / scalor
            m += 1 
        
        # nominator_init =  factorial(n)
        # denominator_init = ( factorial(n)) * (2 * n + 1)
        # nominator_end =  factorial(2 * n)
        # denominator_end = (1) * 2 * ( 2 * n + 1)
        # scalor_init = np.sqrt(nominator_init / denominator_init)
        # scalor_end = np.sqrt(nominator_end / denominator_end)

        grid[n][0] =  ( (2 * n - 1) * np.sin(phi) * grid[n-1][0] - (n - 1) * grid[n-2][0] ) / (n) # first value of each row
        # grid[n][0] = grid[n][0] / scalor_init

        grid[n][n] =  (2 * n - 1) * np.cos(phi) * grid[n-1][n-1]  # last value of each row
        # grid[n][n] = grid[n][n] / scalor_end

    return grid

def ALPS(phi, degree, order):
    # X = lpmn(degree, order, (np.sin(phi)))[0]
    X = lpmn(order, degree, (np.sin(phi)))[0]

    ALPS = np.matrix.transpose(np.array([X[i, :]*(-1)**i for i in range(len(X))]))
    return ALPS



#debugging
if __name__ =='__main__' and 1==1:

    it = 4
    phi = 79.3*np.pi/180
    grid = ALP(phi, degree =6)
    
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

    from scipy.special import lpmn
    n, m = 3, 1


    # X = lpmn(4, 4, (np.sin(phi)))[0]
    # ALPS = np.matrix.transpose(np.array([X[i, :]*(-1)**i for i in range(len(X))]))

    print('func', ALP(phi, 2))
    print( '\nvallad/o',P )
    print('ALPS', ALPS(phi, 2, 2))
  
    # print("vallado recursion", grid[it])
    # print("Vallado analytical", P[it])
    # print(grid)
    # exit()





       

def dU_dr(r, phi, lmbda, C, S, degree, order):
    '''
    Inputs in km and rad where relevant
    '''
    degree +=1 

    derivative = 0
    
    # C = matricise(loc_C)
    
    # S = matricise(loc_S)
    
    # print(degree)
    P =  ALPS(phi, degree, degree)
    for n in range(2, int(degree) ): # - 2):
        m = 0   # order (derivative order)
        
        
        while m <= order:
            derivative +=  -(mu / (r**2)) * ((RE / r)**(n) * (n + 1) * P[n][m] * (C[n][m] * np.cos(m * lmbda) + S[n][m] * np.sin(m * lmbda)))
            # print("r//r={}, mu={}, n={}, m={}, C={}, S={}, P={}, deriv={} ".format(r, mu, n, m, C[n][m], S[n][m], P[n][m], derivative))
            m += 1
    return derivative
            


def dU_dphi(r, phi, lmbda, C, S,  degree, order):
    '''
    Inputs in km and rad
    '''
    degree +=1 
    derivative = 0
    # C = matricise(loc_C)
    # S = matricise(loc_S)
    # print(C)
    P = ALPS(phi,degree, degree)
    # print(P)
    for n in range(2, int(degree) ):#- 2):
        m = 0   # order (derivative order)

        
    
        while m <= order:
            derivative +=  (mu / (r)) * ((RE / r)**(n) * ( P[n][m+1] - m * np.tan(phi) * P[n][m] ) * (C[n][m] * np.cos(m * lmbda) + S[n][m] * np.sin(m * lmbda)))
            # print("phi//r={}, mu={}, n={}, m={}, C={}, S={}, P={}, deriv={} ".format(r, mu, n, m, C[n][m], S[n][m], P[n][m+1], derivative))
            m += 1
    return derivative




def dU_dlmbda(r, phi, lmbda,C, S,  degree, order):
    '''
    Inputs in km and rad
    '''
    degree +=1 
    derivative = 0
    # C = matricise(loc_C)
    # S = matricise(loc_S)
    # print("C:",C)
    # print("S:", S)
    P = ALPS(phi, degree, degree)
    # print("ALP:", P)
    for n in range(2, int(degree) ):  
        m = 0   # order (derivative order)


        while m <= order:
        
            derivative +=  (mu / (r)) * ((RE / r)**(n) * m * P[n][m] * (S[n][m] * np.cos(m * lmbda) - C[n][m] * np.sin(m * lmbda)))
            # print("lmbda//r={}, mu={}, n={}, m={}, C={}, S={}, P={}, deriv={} ".format(r, mu, n, m, C[n][m], S[n][m], P[n][m], derivative))
            m += 1
    return derivative


if __name__ == '__main__' and 1==0:
    #debugging
    phi = 0* np.pi/180
    lmbda = 0 * np.pi /180

    C = matricise(2, max_degree=2)
    S = matricise(3, max_degree=2)
    r = 6771
    deg = 2
    P = ALP(phi, deg)
    U = dU_dr(r, phi, lmbda,C, S, degree=deg)
    U_p = dU_dphi(r, phi, lmbda,C, S, degree = deg)
    U_l = dU_dlmbda(r, phi, lmbda,C, S, degree = deg)
    # print('P', P)
    # print('C', C, 'S', S)
    print( 'U_r', U, 'U_p', U_p, 'U_l', U_l)

    exit()
    # corrective values of derivatives for n = 2
    # corr_dr = (-mu / (r**2)) * 3 * ((RE/r)**2) * ((-0.5) * C[2][0] + 3 * C[2][2])
    # print('corr dr:',corr_dr)
    # corr_dphi = (mu / r) * ( RE / (r))**2 * (P[2][1] * C[2][0] + P[2][2] * C[2][1])
    # print('corr dphi:', corr_dphi)
    # corr_lmb = (mu / r) * ((RE / r)**2) * (P[2][1] * S[2][1] + 2 * P[2][2] * S[2][2])
    # print('corr lmb:', corr_lmb)

    
    

# exit()

def static_egm08(r, phi, lmbda, C, S, deg, order):
    '''
    Compute non-keplerian gravity acceleration based on geocentric Earth-fixed polar input. r in km, phi / lmbda in rad
    '''
    rx, ry, rz = spherical_to_r_cart(r, phi, lmbda)#astropy.coordinates.spherical_to_cartesian(r, phi, lmbda)#  spherical_to_r_cart(r, phi, lmbda)
    # r1, r2, r3 = spherical_to_r_cart(r, phi, lmbda)
    # print('r in static', np.allclose(np.array([rx, ry, rz]), np.array([r1, r2, r3])))
    # print('rotation\n', 'phi', phi*180/np.pi, 'lmbda', lmbda*180/np.pi, 'r', r)
   
    U_r = dU_dr(r, phi, lmbda, C, S, deg, order)
    U_phi = dU_dphi(r, phi, lmbda, C, S, deg, order)
    U_lmbda = dU_dlmbda(r, phi, lmbda, C, S, deg, order) 
    
    a_i = ((1 / r) * U_r - ( rz / ( (r**2) * np.sqrt((rx**2) + (ry**2)) ) ) * U_phi ) * rx - ( ( 1 / ((rx**2) + (ry**2)) ) * U_lmbda ) * ry #- ((mu * rx) / (r**3)) 

    a_j = ((1 / r) * U_r - ( rz / ( (r**2) * np.sqrt((rx**2) + (ry**2)) ) ) * U_phi ) * ry + ( ( 1 / ((rx**2) + (ry**2)) ) * U_lmbda ) * rx #- ((mu * ry) / (r**3)) 

    a_k = (1 / r) * U_r * rz + ( ( np.sqrt((rx**2) + (ry**2)) ) / (r**2) ) * U_phi #- ((mu * rz) / (r**3)) 
    
    acceleration = np.array([a_i, a_j, a_k]) #/ 1000 #transform to km/s^2
    
    return acceleration




    



    

def egm08(t, r_cartesian, v_inertial, C, S, epoch, deg, order):
    '''
    Compute the non-keplerian gravity acceleration in the inertial frame, according to the EGM08 data and Vallado model.
    NOTE: output in km/s^2
    '''
    UT = t / 3600
    UT += 12# time offset--in hours

    #get julian date
    day, month, year = epoch

    JD = julianDate(UT, day, month, year) 


    #Initial Epoch, at the start of the simulation
    EARTH_OFFSET = 2 * np.pi * (0.7790572732640 + 1.00273781191135448 * (JD - 2451545))        

    # angle = EARTH_OFFSET #+ Earth_rate * t # in rad
    r_rot, rot_matrix = applyEarthRotation(t, epoch, r_cartesian) #rotate input vector : from inertial to body-fixed rotation(r_cartesian, angle)
    r, lmbda, phi= spiceypy.reclat(r_rot)#astropy.coordinates.cartesian_to_spherical(r_rot[0], r_rot[1], r_rot[2])## #
    #r, phi, lmbda = LatLong(r_rot)
    a_itrf = static_egm08(r, phi, lmbda, C, S, deg, order) #compute body-fixed acceleration (itrf)
    # print('itrf', a_itrf)
    a_eci =  undoEarthRotation(a_itrf, rot_matrix) #transform acceleration into inertial frame rotation(a_itrf, -angle)

    return a_eci

if __name__ =='__main__' and 0 == 1:
    import math
    #testing of rotations and latlong functions
    r = 8378+406#km
    
    # r = np.array([rx, ry, rz])


    r = np.array([2384.46, -5729.01, 5050.46])
    r_norm = np.linalg.norm(r)
    phi =0* math.pi / 180
    lmbda = 0* math.pi / 180
    
    rx = r_norm * np.cos(lmbda)* np.cos(phi)
    ry = r_norm * np.cos(phi) * np.sin(lmbda)
    rz = r_norm * np.sin(phi)
    C = matricise(2)
    S = matricise(3)
    
    
    v = np.array([-7.36138,  -2.98997, 1.64354])
    time = 0
    # angle = np.linalg.norm(VEC_EARTH_ROTATION) * time
    # v_rot = v_rotation(v, rotation(r, angle))
    # print('v', v)
    # print('v_rot', v_rot)
    # print('omgea earth', VEC_EARTH_ROTATION[2])
    # print(v_rot[1] * VEC_EARTH_ROTATION[2] * 2)

    
    acc = egm08(time, r, v, C, S, np.array([0, 0, 0]), 2)
    # tot_acc = np.linalg.norm(acc)
    # print('input r_norm', r_norm)
    print('acc vec', acc)
    # print('acc norm', tot_acc)
    # print(v_rot[0] * VEC_EARTH_ROTATION[2] * 2)
    # print('theoretical value', r *  mu / np.linalg.norm(r**3))


    #comparing to analytical values of Curtis
    # J2 = 0.001082626173852227
    # P1 = (3/2) * (J2 * mu * RE**2) / (r_norm**4)
    # P2 = (rz**2) / (r_norm**2)

    # P_curtis_analytical = np.array([P1 * ( (rx/r_norm) * (5 * P2 - 1) ), P1 * ( (ry/r_norm) * (5 * P2 - 1) ), P1 * ( (rz/r_norm) * (5 * P2 - 3) )])
    # print('Curtis analytical:', P_curtis_analytical)
    # print('Norm Curtis:', np.linalg.norm(P_curtis_analytical))

    