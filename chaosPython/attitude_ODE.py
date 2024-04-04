#python code to define attitude ODE.


'''
This code defines the function(s) representing the ODE for the attitude dynamics of spacecrafts.
This is based of the Euler attitude equations and the quaternion kinematics and allows for calculating the velocity of a spacecraft when thrusting with given torque and initial
rotational velocity.
'''
#import relevant modules and files
import numpy as np




# Developed form
def analytical(t, omega_state, torque_thruster, torque_perturb, satellite,  grid):
    #get omega components and quaternion --- quaternion is first part of array
    quaternion = omega_state[:4]
    omega_rot = omega_state[4:7]# omega rot in body-fixed frame.
    B_hyst = omega_state[7:10]
    # grid.burnTime[pixel] =  omega_state[10]
    #create derivative array
    omega_state_dot = np.zeros_like(omega_state)
    
    # print('attitude ODE: state', omega_state)
    
    #update values for satellite class 
    satellite.quaternion = quaternion
    satellite.angularVel = omega_rot
    satellite.hysteresisInduction = B_hyst
    wx, wy, wz = omega_rot 
    

    #extract other constants:
    inertia_inv = satellite.InverseinertiaMatrix

    I = satellite.inertiaMatrix
    #get orbital position and velocity
    r_eci = satellite.eci_position() 
    v_eci = satellite.eci_velocity()

    #get magnetic dipole data
    #deprecated 
    # B_inertial, B_inertial_dot = dipole_model(t, r_eci, v_eci)
 
 
    #call perturbing torques and sum them to thruster torques
    #Magnetic function deprecated
    # torque_mag, B_dot = satellite.perturbingAccMag(B_inertial, B_inertial_dot, omega_state, B_hyst, t)
    # # orbital_perturb, torque_perturb = satellite.envPerturbations(t)
    B_dot = np.zeros_like(B_hyst)
    #sum torques, SI
    torque = torque_thruster + torque_perturb #sum all torques// +  torque_mag 
    # print('Att ODE, thruster Torques', torque_thruster)
    tx, ty, tz = torque

    #compute angular acceleration
    # matrix = np.array([tx + (I[1][1] - I[2][2]) * wy * wz, ty + (I[2][2] - I[0][0]) * wx * wz, tz + (I[0][0] - I[1][1]) * wy * wx])
    matrix = torque - np.cross(omega_rot, np.matmul(I, omega_rot))

    omega_dot = np.matmul(inertia_inv, matrix)
    # print('ode', omega_dot)

    #compute quatenion rate of change - numerical
    lmbda_matrix = np.array([[0, wz, -wy, wx], [-wz, 0, wx, wy], [wy, -wx, 0, wz], [-wx, -wy, -wz, 0]])
    q_dot = 0.5 * np.matmul(lmbda_matrix, quaternion)
    
    #fuel ODE:
    fuel_consumption = -1 * grid.state #do not comsume when grid is off
    #debug statement
    # print('quaternion', quaternion, 'omega rot', omega_rot, 'rot matrix', lmbda_matrix, 'q_dot', q_dot, 'omega_dot', omega_dot, 'inertia matrix', I)
 
    # Merge both array for integration [q_dot, omega_dot]

    #populate the derivative array
    omega_state_dot[:4] = q_dot #quaternion
    omega_state_dot[4:7] = omega_dot # angular vel
    omega_state_dot[7:10] = B_dot #magnetic ---deprecated
    omega_state_dot[10] = fuel_consumption #fuel consumption
    omega_state_dot[11] = 0 #just track the state, do not alter it
    omega_state_dot[-1] = 0 #just track the state, do not alter it
        #progress:
    # print('progress', t/t_stop ) 
    # print('acc:', omega_dot[:4])
    return omega_state_dot





