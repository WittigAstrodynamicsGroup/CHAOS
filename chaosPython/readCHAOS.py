# PYTHON FILE TO READ THE DATAFILES CREATED BY CHAOS.PY




'''
This python code is procedural code to extract and analyse the data generated when running the C.H.A.O.S. code.
'''

import os
import numpy as np

from .class_extractor import read_data, read_grid_data
from .quaternions import body_to_inertial
import matplotlib.pyplot as plt
from .transformations import eq_to_kep, equi_to_r_cart, equi_to_v_cart
import time
#change path to where the data is

def plotData(it, datapath="./output", figpath='./figs/'):

    print('Looking for data in ', datapath)
    print('Figures will be stored in', figpath)


    #include the run in figpath
    figpath = figpath + 'run{}'.format(it)
    if os.path.exists(figpath) != True:
        os.mkdir(figpath)
    # os.chdir('Total combined')
    # #extract velocity and quaternion data
    start = time.time()


    t, wx, wy, wz, w_norm = read_data(datapath + '/run{}/timetrack/timetrack{}.txt'.format(it, it))

    t1, qx, qy, qz, qr = read_data(datapath + '/run{}/attitude/quaternion{}.txt'.format(it, it))

    t_grid, gridname, grid_state, grid_fuel, total_fuel = read_grid_data(datapath +'/run{}/objects/grid1_data{}.txt'.format(it, it))
    t_grid2, gridname2, grid_state2, grid_fuel2, total_fuel2 = read_grid_data(datapath +'/run{}/objects/grid2_data{}.txt'.format(it, it))
    #extract orbital data
    ts_orb = []
    p = []
    f = []
    g = []
    h = []
    k = []
    L = []
    testItem1 = []
    testItem2 = []
    testItem3 = []
    testItem4 = []

    datafile = open(datapath + '/run{}/orbital/equinoctial{}.txt'.format(it, it), 'r')
    for line in datafile.readlines(): #split each line
        item = line.split('\t') #split each element
        p.append(float(item[1]))
        f.append(float(item[2]))
        g.append(float(item[3]))
        h.append(float(item[4]))
        k.append(float(item[5]))
        L.append(float(item[6]))



    print('Reading files took (s)', time.time()-start, 'file length', len(t))
    start = time.time()

    t_o = np.append([], t)
    t_o_s = np.append([], ts_orb)
    t = t_o/3600
    t_h_s = t_o_s /3600




    # #compute quatenrion norm for error analysis

    q_norm = np.sqrt(np.square(qx) + np.square(qy) + np.square(qz) + np.square(qr)) - 1
    print('computing quaternion (s)', time.time()-start)
    start = time.time()
    ## transform to angles/s for comparison:

    wx = np.append([], wx )
    wy = np.append([], wy )
    wz = np.append([], wz )
    w_norm = np.append([], w_norm)

    # #rad to deg
    wx = wx * 180/np.pi
    wy = wy * 180/np.pi
    wz = wz * 180/np.pi
    w_norm = w_norm * 180/np.pi


    I = np.array([[0.0067,0,0], [0,0.0419,0], [0,0,0.0419]])
    m = 0.5
    # #compute pointing error between z-axis and ----hard-coded--- magnetic field vector 
    z = np.array([0, 0, 1])
    y = np.array([0, 1, 0])
    x = np.array([1, 0, 0])
    H_cste = np.array([0,  0, 22464.7e-9 ])
    pointingz_list = []
    pointingy_list = []
    pointingx_list = []
    magnitude_list = []
    SMA = []
    ECC = []
    INC = []
    RAAN = []
    WP = []
    TA =[]
    ra = []
    rp = []
    offset = []
    rel_v_dir_x = []
    rel_v_dir_y = []
    rel_v_dir_z = []

    thruster_switch_on = []
    thruster_switch_off = []

    print(len(q_norm))

    a0, ecc0, inc0, raan0, wp0, ta0 = eq_to_kep(p[0], f[0], g[0], h[0], k[0], L[0])
    for i in range(len(q_norm)):
        print(i)
        q = np.array([qx[i], qy[i], qz[i], qr[i]])
        z_inertial = body_to_inertial(z, q)
        y_inertial = body_to_inertial(y, q)
        x_inertial = body_to_inertial(x, q)



        #orbital data conversion
        a, ecc, inc, raan, wp, ta = eq_to_kep(p[i], f[i], g[i], h[i], k[i], L[i])

        SMA.append(a)
        ECC.append(ecc)
        INC.append((inc)*180/np.pi)
        RAAN.append((raan)*180/np.pi)
        WP.append((wp)*180/np.pi)
        TA.append(ta*180/np.pi)
        ra = np.append(ra, a * (1 + ecc))
        rp = np.append(rp, a* (1 - ecc))

        v_cart = equi_to_v_cart(np.array([p[i], f[i], g[i], h[i], k[i], L[i]]))
        r_cart = equi_to_r_cart(p[i], f[i], g[i], h[i], k[i], L[i])
        v_dir = (v_cart / np.linalg.norm(v_cart))
        # angle between axis and velocity direction--velocity tracking
        vcossy = np.dot(y_inertial, v_dir)
        vcossz = np.dot(z_inertial, v_dir)
        vcossx = np.dot(x_inertial, v_dir)
        vpointingz = np.arccos(vcossz / (np.linalg.norm(z_inertial) * np.linalg.norm(v_dir)))
        vpointingy = np.arccos(vcossy / (np.linalg.norm(y_inertial) * np.linalg.norm(v_dir)))
        vpointingx = np.arccos(vcossx / (np.linalg.norm(x_inertial) * np.linalg.norm(v_dir)))
        rel_v_dir_z.append(vpointingz*180/np.pi)
        rel_v_dir_y.append(vpointingy*180/np.pi)
        rel_v_dir_x.append(vpointingx*180/np.pi)

        




        #Thruster on/off code:
        if rel_v_dir_y[-1] <= 90:
            thruster_switch_on.append(1)

        else:thruster_switch_on.append(0)







    print('parsing took (s)', time.time()-start)
    #############################Test variables##################

    plt.scatter(t_grid, total_fuel, label=r'fuel evolution')
    plt.xlabel('Time [hrs]')
    plt.ylabel(r'remaining firing time [s]')
    # plt.title('event function')
    plt.legend(loc=0)
    plt.savefig(figpath + "/Total_fuel_time_evolution.pdf")
    plt.clf()

    plt.scatter(t_grid, grid_fuel, label='grid1 fuel')
    plt.scatter(t_grid, grid_fuel2, label='grid2 fuel')
    plt.xlabel('Time')
    plt.ylabel('Fuel remaining [time]')
    plt.savefig(figpath + "/grid_fuel_time_evolution.pdf")
    plt.clf()

    plt.scatter(t_grid, grid_state, label='grid1 state', c='firebrick')
    plt.scatter(t_grid, grid_state2, label='grid2 state2', c='forestgreen')
    plt.xlabel('Time')
    plt.ylabel('Grid state [on/off]')
    plt.savefig(figpath + "/grid_fuel_time_evolution.pdf")
    plt.clf()


    #############################Orbital variables##################
    #plot orbital data:
    #SMA
    plt.plot(t, SMA, label='sma')
    plt.xlabel('Time [hrs]')
    plt.ylabel('kilometers')
    plt.legend()
    plt.savefig(figpath + "/SMA_vs_time.pdf")
    plt.clf()




    #ra/rp:
    plt.plot(t, ra, label='apogee')
    plt.plot(t, rp, label='perigee')
    plt.xlabel('Time [hrs]')
    plt.ylabel('kilometers')
    plt.legend()
    plt.savefig(figpath + "/perigee_apogee_vs_time.pdf")
    plt.clf()

    #ecc
    plt.plot(t, ECC, label='ecc')
    plt.xlabel('Time [hrs]')
    plt.ylabel('[-]')
    plt.legend()
    plt.savefig(figpath + "/ecc_vs_time.pdf")
    plt.clf()

    plt.plot(t, INC, label='inc')
    plt.xlabel('Time [hrs]')
    plt.ylabel('degrees')
    plt.legend()
    plt.savefig(figpath + "/inc_vs_time.pdf")
    plt.clf()

    plt.plot(t, RAAN, label='raan')
    plt.xlabel('Time [hrs]')
    plt.ylabel('degrees')
    plt.legend()
    plt.savefig(figpath + "/raan_vs_time.pdf")
    plt.clf()

    plt.plot(t, WP, label='wp')
    plt.xlabel('Time [hrs]')
    plt.ylabel('degrees')
    plt.legend()
    plt.savefig(figpath + "/AOP_vs_time.pdf")
    plt.clf()

    plt.plot(t, TA, label='ta')
    plt.xlabel('Time [hrs]')
    plt.ylabel('degrees')
    plt.legend()
    plt.savefig(figpath + "/TA_vs_time.pdf")
    plt.clf()


    #############################Attitude variables##################
    # plot results-qauternion norm 
    plt.plot(t, q_norm,label='Quaternion norm')
    plt.xlabel('time [hrs]')
    plt.ylabel('Quaternion norm [-]')
    plt.legend()
    plt.savefig(figpath + "/quaternion_norm.pdf")
    plt.clf()

    # angular velocity norm 
    plt.scatter(t, w_norm,label='velocity norm', s=0.1, c='darkorange')
    plt.xlabel('time [hrs]')
    plt.ylabel('angular velocity norm [deg/s]')
    plt.legend()
    plt.savefig(figpath + "/angular_Velocity_norm.pdf")
    plt.clf()


    # angular velocity components
    plt.plot(t, wx,label='wx' )
    plt.plot(t, wy, label='wy')
    plt.plot(t, wz,label='wz' )
    plt.xlabel('time [hrs]')
    plt.ylabel('angular velocity [deg/s]')
    plt.legend()
    plt.savefig(figpath + "/angVel_components.pdf")
    plt.clf()



