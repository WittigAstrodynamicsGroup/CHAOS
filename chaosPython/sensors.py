# sensor trial script. THis python script tests and 
# validates the uses of various sensors 

import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt

from .EnvironmentModelling import rhoIonsFunc

def faradayCupCurrent(alpha, v, A, rho_ions, k):
    q_e = 1.6e-19 #Coulomb
    I = q_e * rho_ions * A * v * np.cos(alpha)
    if I < 0:    
        return 0
    else:
        return I

def faradayCupErrorCst(alpha, k, epsilon, v, A, rho_ions):
    #function to return the measurement error based on the error fraction in the current

    q_e = 1.6e-19 #Coulomb

    # err = ((1/k) * np.arccos( np.cos(alpha)  * (1 + epsilon))) - alpha
    err = ((1/k) * np.arccos( np.cos(alpha)  + (epsilon / (q_e * A * v * rho_ions)))) - alpha
    return err

def faradayCupAngle(alpha, epsilon, v, A, rho_ions, rng):
    #function to return the measured angle based on the error of the faraday cup
    q_e = 1.6e-19 #Coulomb
    #random value of current error:
    rand = rng.normal(0, epsilon)
    #peak current:
    I_max = (q_e * A * v * rho_ions)
    
    internal = np.cos(alpha)  + ( rand / I_max )
    # if current greater than peak, limit to peak
    if abs(internal) > 1:
        internal = 1 * np.sign(internal)
    
    #set faraday cup detection limit
    # if alpha >= limit:
    #     internal -= 0.5

    angle =  np.arccos( internal ) 
    return angle


def effectiveCupAngle(SN, epsilon, rho_ions_total, aperture, v):
    q_e = 1.6e-19 #Coulomb

    ratio = q_e / epsilon

    angle = np.arccos(SN /( (ratio) * rho_ions_total * aperture * v) ) 

    return angle #rad


def signalToNoise(rho_ions , A , v, epsilon):
    q_e = 1.6e-19

    SN = rho_ions * A * v *q_e / epsilon

    return SN

def MEMSgyro(omega, delta_t, biais_n, sigma_n, sigma_b, rng, rng2):
    
    rand_b = rng.normal(0, 1, size=3) 
    rand_n = rng2.normal(0, 1, size=3) 
    biais_new = biais_n + sigma_b * np.sqrt(delta_t) * rand_b
    measurement = omega + 0.5 * (biais_new + biais_n) + np.sqrt(((sigma_n**2)/delta_t) + ((delta_t*(sigma_b**2))/12)) * rand_n
    # print('biais',biais_new, delta_t)
    return biais_new, measurement



def scalar_MEMSgyro(omega, delta_t, biais_n, sigma_n, sigma_b, rng, rng2):
    
    rand_b = rng.normal(0, 1) 
    rand_n = rng2.normal(0, 1) 
    biais_new = biais_n + sigma_b * np.sqrt(delta_t) * rand_b
    measurement = omega + 0.5 * (biais_new + biais_n) + np.sqrt(((sigma_n**2)/delta_t) + ((delta_t*(sigma_b**2))/12)) * rand_n
    # print(measurement)
    return biais_new, measurement


if __name__ == '__main__':


    # #MEMS:
    # def trendline(sigma, t):
    #     return np.sqrt(t) * sigma


    # time_range = np.linspace(0, 86400*1000, 10000)
    # itera = 250
    # omega = 0
    # biais_n = 0
    # sigma_n = ((0.15)/60) * np.pi/180 #rad
    # sigma_b = ((0.3)/3600) * np.pi/180 #rad
    # prev = 0
    # biais_n1 = biais_n
    # biais_n2 = biais_n
    # biais_n3 = biais_n
    # biais_n4 = biais_n
    # dict = {}
    # randx = []
    # randy = []
    # randz = []
    # chi = []
    # rng = default_rng()

    # for j in range(itera):
    #     prev = 0
    #     biais_n1 = biais_n
    #     dict[j] = []

    #     rng1 = default_rng()
    #     rng2  =default_rng()
    #     rng3 = default_rng(seed=5)
    #     for i in range(len(time_range)):
    #         delta_t = time_range[i] - prev
    #         # omega = np.sin(0.1*i) 
    #         val1 = scalar_MEMSgyro(omega, delta_t, biais_n1, sigma_n, sigma_b, rng1, rng2)
    #         dict[j].append(val1[1] *180/np.pi)
    #         prev = time_range[i]
    #         biais_n1 = val1[0]





    #     plt.plot(time_range/86400, dict[j], c='green', alpha=0.05, zorder=1)
    
 
    # std = []
    # mn = []
    # for i in range(len(time_range)):
    #     sigma = []
    #     for j in range(itera):
    #         sigma.append(dict[j][i]) #already in deg
    #     std.append(np.std(sigma))
    #     mn.append(np.mean(sigma))


    # plt.plot(time_range/86400, trendline(sigma_b, time_range)*180/np.pi, linestyle='dashed', c='darkorange',zorder=3, label='predicted std')
    # plt.plot(time_range/86400, trendline(-sigma_b, time_range)*180/np.pi, linestyle='dashed', c='red', zorder=2,label='predicted std')
    # plt.scatter(time_range/86400, std, marker='x', c='black', zorder=2, label='computed std')
    # plt.scatter(time_range/86400, mn, marker='o', c='red', zorder=2, label='computed mean')
    # plt.xlabel('days')
    # plt.ylabel('angular velocity [deg/s]')
    # plt.legend()
    # # plt.savefig('/home/kash/Desktop/PhD/Work/report-writing/Plots/CHAOS/applications/ESA/MEMS_gyro/VG1703SPE_msrmt_{}.pdf'.format(itera))
    # plt.savefig('localPlot/sensor.png')

    #faraday cup



    EGG_current = np.array([0.6, 0.58, 0.61, 0.5, 0.2, 0, 0])

    EGG_angle = np.array([0, 15, 30, 45, 60, 75, 90])

    mu_E = 3.986e14
    aperture =  3.85e-4#np.pi * ((30/2)*1e-3)**2 #np.pi * (7e-3)**2
    test_rho = 3.75e12 / (3.85e-5 * 1e3)
    rho_ion = 1e11

    f_rhoO, f_rhoH = rhoIonsFunc()
    # print(test_rho)
    # exit()
    # v_sc=7.7e3
    angle_range = np.linspace(0, np.pi/2, 300)
   
    alt_range = np.linspace(150, 2000, 300)

    SN_mat = np.zeros((len(alt_range), len(angle_range)))
    current = []
    total_rho = []
    rho_H_plot = []
    rho_O_plot = []
    angle_one = []
    eff_thrust = []
    eff_thrust_2 = []
    false_firing = []
    #altitudes
    for i in range(len(alt_range)):
        print(i)
        # compute  velocity
        a = 6378 + alt_range[i] 
        v_sc = np.sqrt(mu_E /( a*1e3)) #m/s
        print('a', a, 'v', v_sc)
        #compute rho ions

        rho_O = f_rhoO(alt_range[i])
        rho_H = f_rhoH(alt_range[i])

        angle_one.append(np.arccos(1 / (faradayCupCurrent(0, v_sc, aperture, (rho_O + rho_H)*1e6, 1) / 1e-9)) *180/np.pi)
        eff_thrust.append((1  - np.cos(angle_one[-1]*np.pi/180)**2) / 4)
        eff_thrust_2.append(eff_thrust[-1]*2)
        false_firing.append((1 - (1 - np.cos(angle_one[-1]*np.pi/180))) / 2)

        current.append(faradayCupCurrent(0, v_sc, aperture, (rho_O + rho_H)*1e6, 1))
        rho_H_plot.append(rho_H)
        rho_O_plot.append(rho_O)
        total_rho.append(rho_H + rho_O)


        # for j in range(len(angle_range)):
        #     print(j)
        #     I = faradayCupCurrent(angle_range[j], v_sc, aperture, (rho_O + rho_H)*1e6, 1)
        #     # print('I', I)
        #     SN = I / 1e-9
        #     if SN < 0.8:
        #         SN = 0
        #     SN_mat[i][j] = SN

            # coordinates[i][j] = (i, j)



    # print(SN_mat)
    plt.plot(alt_range, current)
    plt.xlabel('Altitude [km]')
    plt.ylabel('Current [A]')
    plt.show()
    plt.plot(total_rho, alt_range, label='total ion density')
    plt.plot(rho_H_plot, alt_range, label='H+ density')
    plt.plot(rho_O_plot, alt_range, label='O+ density')
    plt.ylabel('Altitude [km]')
    plt.xlabel('Ion densities [per m^3]')
    plt.legend()
    plt.show()
    plt.plot(alt_range, angle_one)
    plt.xlabel('Altitude [km]')
    plt.ylabel('Cone angle for S-to-N=1 [deg]')
    plt.show()    
    plt.plot(alt_range, eff_thrust, label='single-sided')
    plt.plot(alt_range, eff_thrust_2, label='double-sided')
    plt.xlabel('Altitude [km]')
    plt.ylabel('Effective thrust coefficient [-]')
    plt.legend()
    plt.show()
    plt.plot(alt_range, false_firing)
    plt.xlabel('Altitude [km]')
    plt.ylabel('Poportion of acceleration thrust [-]')
    plt.show()
    # plt.imshow(SN_mat)
    # plt.contourf(angle_range*180/np.pi, alt_range, SN_mat, cmap='tab20c')
    # plt.colorbar()
    # plt.show()