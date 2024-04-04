#python code storing the class for the satellite object


'''
Satellite class definition. The satellite object holds the attributes of the spacecraft, such as 

    -mass: Mass of the object [kg]
    -Cd: Coefficient of Drag [-]
    area: Can be a constant, although it will be implemented as a function at some point [m^2]
    quaternion: Currently None, but will hold the quaternion/ pointing vector of the satellite and its angular velocity.
    velocity: Currently None, but will be used when implemented with the propagator (orbital state, x, x_dot) [km, [km/s]]
    inertiaMatrix: inertia matrix of the satellite, which holds infomration about the shape and mass distribution of the object [SI]



NOTE: Default values for a 1U CubeSat. Do not confuse inertiaMatrix with InverseinertiaMatrix!!!
'''

import numpy as np
from stl import mesh
from .transformations import kep_to_eq, rv2kpl, kpl2rv, equi_to_r_cart, equi_to_v_cart
from .quaternions import body_to_inertial, inertial_to_body
from .support_fun import gridsBurntime
import warnings


class Satellite:
    #keyword input argument--default values for standard 1U CubeSat.
    def __init__(self, grids, environment,  **kwargs):
    
        #get environmental info
        self.environment = environment
        self.grids = grids
        #satellite state information
        self.quaternion = kwargs.get('quaternion', None)
        self.angularVel = kwargs.get('angularVel', None)
        self.velocity = kwargs.get('velocity', None) #orbital
        self.position = kwargs.get('position', None)
        self.keplerianElements = kwargs.get('keplerianElements', np.array([6781.36, 0.0021168, 51.6*np.pi/180, 90*np.pi/180, 0, 0]))
        self.modEquinoctialElements = kwargs.get('equinoctial', None)


        #satellite data
        self.mass = kwargs.get('mass', 1.33)
        self.Cd = kwargs.get('Cd', 2.2)
        self.Cr = kwargs.get('Cr', 0.8) #reflection, for SRP
        self.area = kwargs.get('area', 0.015)# m^2 #lf MODEL
        self.Tw = kwargs.get('Tw', 300) #k
        self.alpha_coeff = kwargs.get('alpha_coeff', 0.9)


        #Geometry definition
        self.panelList = kwargs.get('panelList', [])
        self.STLfile = kwargs.get('STLfile', None)
        self.scale_factor = kwargs.get('scale_factor', 1e6) #how many units of STL file in 1 m2?
        self.minFacet = kwargs.get('minFacet', 4.6e-3)
        self.inertiaMatrix = kwargs.get('inertiaMatrix', np.diag([0.00216667]*3))
        self.InverseinertiaMatrix = np.linalg.inv(self.inertiaMatrix)
        self.centreOfMass = kwargs.get('centreOfMass', np.array([0, 0, 0]))

        #magnetic perturbation data
        self.hysteresisOrientation = [np.array([0, 0,1])] #orientation of hsyterisis material, in body-fixed frame
        self.hysteresisVolume = np.array([0])
        self.hysteresisInduction = [np.array([0])]
        self.hysteresisMoment = np.array([0])
        self.hysteresisSaturation = np.array([0.7])
        self.hysteresisRemanence = np.array([0.4])
        self.hysteresisCoercive = np.array([1.59])
        self.B0 = np.array([0])  
        self.MagnetMoment = [np.array([0, 0, 0])]
        self.t_permMagnet = 0
        self.t_hysteresis = 0
        self.totalPerturbation = 0
        self.t_previous = 0
        

        #Constants
        self.mu_fs = (4e-7) * np.pi #free space permeability
        self.mu = environment.mu #398600.4418   #km-5-vallado
        self.mu_E = self.mu * 1e9 # Earth gravity parameters
        self.epsilon = 1.1  #tolerance for value of B_hyst
        self.it = None #iteration number for writing


        #initialisation of orbital state
        #Based on which variable is given, transform to modifiedEquinoctial
        if self.position is not None and self.velocity is not None:
            #get keplerian elements -r,v in km - km/s
            self.keplerianElements = rv2kpl(self.position, self.velocity)

        else:
            # get eci r, v in km-km/s
            self.position, self.velocity = kpl2rv(self.keplerianElements[0],        #SMA
                                                    self.keplerianElements[1],      #ECC
                                                    self.keplerianElements[2],      #INC 
                                                    self.keplerianElements[3],      #RAAN 
                                                    self.keplerianElements[4],      #AOP 
                                                    self.keplerianElements[5],      #THETA
                                                    self.mu)                        #MU--in km^3/s^2
        
            #transform into modEq
        a, e, inc, raan, wp, theta = self.keplerianElements
        self.modEquinoctialElements = kep_to_eq(a, e, inc, raan, wp, theta)



        #initialisation of geometry
        #using custom panel class
        if len(self.panelList) != 0:
            self.areas = np.array([i.area for i in self.panelList])
            self.normals = [i.surfNormal for i in self.panelList]
            self.positions = [i.surfPos for i in self.panelList]
            self.Tws = [i.Tw for i in self.panelList]
            self.Crs = [i.Cr for i in self.panelList]
            self.alpha_coeffs = [i.alpha for i in self.panelList]

        if self.STLfile is not None:
            #read stl file
            self.minFacet = self.minFacet * self.scale_factor #minFacet given in m2
            HFstl = mesh.Mesh.from_file(self.STLfile)

            self.areas = np.array([HFstl.areas[i][0] / self.scale_factor for i in range(len(HFstl.areas))])
            #[HFstl.areas[i][0] / self.scale_factor for i in range(len(HFstl.areas)) if HFstl.areas[i][0] >= self.minFacet]

            self.normals = HFstl.get_unit_normals()
            self.positions = [np.array([np.sum(HFstl.x[i])/3, np.sum(HFstl.y[i])/3, np.sum(HFstl.z[i])/3])/np.sqrt(self.scale_factor) for i in range(len(HFstl.x))]
            self.Tws = [self.Tw for i in range(len(HFstl.areas))]
            self.Crs = [self.Cr for i in range(len(HFstl.areas))]
            self.alpha_coeffs = [self.alpha_coeff for i in range(len(HFstl.areas))]

            if len(HFstl.areas) > 30:
                warnings.warn("WARNING: {} facets detected. Is this normal? Slow run time can be expected.".format(len(HFstl.areas)))


            
    
    
    #define any methods here


    def envPerturbations(self, t):

        #fetch the spacecraft state
        vec_r_eci = self.eci_position()
        vec_v_eci  =self.eci_velocity()
        quaternion = self.quaternion

        disturb_force = np.array([0., 0., 0.])
        disturb_torque = np.array([0., 0., 0.])

        #compute pure orbital forces: EGM  + LS
        disturb_force += self.environment.orbitalPerturbation(t, vec_r_eci, vec_v_eci, self.area, self.mass, self.Cd, self.Cr)
        
        # print('pure', disturb_force)
        #check that a HF model is implemented
        # if len(self.panelList) != 0:
            #HF computation of panel-dependent forces


        #compute the torque and force for each panel/facet
        # print('sat t', t)


        #for each grid, compute parameters
        totalFuelMass = 0
        totalBurntime = 0
        burntimeVal = gridsBurntime(self.grids)
        for grid in self.grids:
            totalFuelMass += grid.fuelMass
            totalBurntime += grid.totalburnTime

        #print('fuel mass', totalFuelMass)
        #compute  fuel consumption:
        consumed_fuel = totalFuelMass*(1 - (burntimeVal / totalBurntime))
        remaining_mass = self.mass - consumed_fuel

        # print('remaining mass', remaining_mass)


        #compute "effective CoM"
        COMpos = self.positions - self.centreOfMass
        force, torque = self.environment.panelPerturbation(t, vec_r_eci, vec_v_eci, remaining_mass, quaternion, self.areas, self.normals, self.Tws, self.alpha_coeffs, self.Crs, COMpos)
            # sums up
        disturb_force += force
        disturb_torque += torque
        # # # else:

        #call low fidelity drag and srp forces
        # disturb_force += self.environment.LFpanelPerturbation(t, vec_r_eci, vec_v_eci, self.area, self.mass, self.Cd, self.Cr)
        


        return disturb_force, disturb_torque

    def geometrySurfacePosition(self):
        #returns a matrix with all 6 vectors defining the surface positions. defined
        # as a  method as it needs to be updated post initialisation
        #NOTE: TO DEL--DEPRECATED

        # position of X, -X, Y, -Y, Z, -Z
        posMatrix = np.array([self.geometryPosX, -self.geometryPosX, self.geometryPosY, -self.geometryPosY, self.geometryPosZ, -self.geometryPosZ])
        normalMatrix = np.array([self.geometryX, -self.geometryX,self.geometryY, -self.geometryY, self.geometryZ, -self.geometryZ])
        surfMatrix = np.array([self.geometrySurfX, self.geometrySurfX, self.geometrySurfY, self.geometrySurfY, self.geometrySurfZ, self.geometrySurfZ])
        
        return surfMatrix, posMatrix, normalMatrix



    def eci_position(self):
        #compute eci position from modified equinoctial elements
       
        p, f, g, h, k, L = self.modEquinoctialElements          ##Get the equinoctial elements
        
        r_ECI = equi_to_r_cart(p, f, g, h, k, L)                ##Function in "transformations.py"
        
        self.position = r_ECI                                   ##Update the class attribute

        return r_ECI






    def eci_velocity(self):
        #compute eci velocity from modified equinoctial elements
        p, f, g, h, k, L = self.modEquinoctialElements
        # print('eci_vel, mod equi:', self.modEquinoctialElements)
        w = 1 + f * np.cos(L) + g * np.sin(L)
        s = 1 + (h**2) + (k**2)
        alpha = (h**2) - (k**2)

        v_ECI = np.zeros(3)
        v_ECI[0] = -(np.sin(L) + (alpha) * np.sin(L) - 2 * h * k * np.cos(L) + g - 2 * f* h * k + (alpha) * g)
        v_ECI[1] = -(-np.cos(L) + (alpha) * np.cos(L) + 2 * h * k * np.sin(L) - f + 2 * g * h * k + (alpha) * f)
        v_ECI[2] = 2 * (h*np.cos(L) + k*np.sin(L) + f*h + g*k)

        v_ECI = v_ECI * np.sqrt(self.mu/p)/(s)
        self.velocity = v_ECI
        # print('v_Eci', v_ECI)
        return v_ECI






    #magnetic methods here

    # def initial_induction(self):
    #     Bs = self.hysteresisSaturation[-1]
    #     Br = self.hysteresisRemanence[-1]
    #     Hc = self.hysteresisCoercive[-1]
    #     u = self.hysteresisOrientation[1:]
    #     a, e, i, RAAN, wp , theta = self.keplerianElements
    #     r_eci = self.position
    #     v_eci = self.velocity

    #     #call magnetic field:
    #     B_env, B_env_dot = dipole_model(0, r_eci, v_eci)
    #     warnings.warn('initial Epoch not specified in Dipole model. This msg is in class_satellite. The attitude is randomised here')
    #     angle_min = 0 * np.pi / 180
    #     angle_max = 360 * np.pi / 180
    #     angle = angle_max * 2
    #     if self.quaternion is None:

    #         while angle < angle_min or angle > angle_max:

    #         #intialise random quaternion
    #             v = np.random.rand(1,3) 
    #             angle = np.random.rand(1) * 2 * np.pi
    #             v = v / np.linalg.norm(v)
    #             q = np.append(v * np.sin(angle/2), np.cos(angle/2))
    #             q = q / np.linalg.norm(q)

    #             #compute z angle
    #             z = np.array([1, 0, 0])
    #             z_in = body_to_inertial(z, q)
    #             coss = np.dot(z_in, B_env) / np.linalg.norm(B_env)
    #             angle = np.arccos(coss) 
    #     else:
    #         q = self.quaternion 
    #     #compute H_parallel
    #     H_env = B_env / (4*np.pi * 1e-7)
    #     p =  (1/Hc) * np.tan((np.pi * Br) / (2 * Bs))
    #     H_parallel = []
    #     for i in range(len(u)):
    #         H_parallel.append(np.dot(H_env, body_to_inertial(u[i], q)))
    #         #compute associated random value of B:
    #         #Upper value
    #         B_lim1 = (2/np.pi) * Bs * np.arctan(p * (H_parallel[i] - Hc))
    #         #lower value
    #         B_lim2 = (2/np.pi) * Bs * np.arctan(p * (H_parallel[i] + Hc))

    #         it = [B_lim1, B_lim2]
    #         max_value = max(it)
    #         min_value = min(it)

    #         #random--- assign the induction to the satelltie attribute
    #         self.hysteresisInduction.append((max_value - min_value) * np.random.default_rng().random() + min_value)

    #     #assign randomised orientation
    #     self.quaternion = q

    #     return q   


    def computeMoment(self, B, V, u):
        # compute the magnetic dipole moment based on given data --body-fixed frame // Am^2
        m = ( (B * V) / self.mu_fs ) * u
        # print('Inputs: \t', 'B', B, 'V', V,'u', u)
        # print('COmputeMoment:', m)
        return m

    def permanentMagnet(self, B, V, u):
        #appends the permanent magnetic dipole value to the list
        self.MagnetMoment.append(self.computeMoment(B, V, u))
        

    def hysteresisMaterial(self, V, Bs, Br, H0, u):
        #appends the hysteresis factors values to list ---body-fixed frame
        self.hysteresisCoercive = np.append(self.hysteresisCoercive, H0)
        self.hysteresisRemanence = np.append(self.hysteresisRemanence, Br)
        self.hysteresisSaturation = np.append(self.hysteresisSaturation, Bs)
        self.hysteresisVolume = np.append(self.hysteresisVolume, V)
        self.hysteresisOrientation.append(u) #body-fixed frame
        self.B0 = [0] * len(self.hysteresisVolume)
        self.B0 = np.append([], self.B0)

    def dH_dt(self, r, v, m):
        #rate of change of the magnetic field for a genenric orbit in an invariant dipole model //H
        pass
        

    def dH_parallel_dt(self, omega, quaternion, H_in, H_dot_in):
        #compute the rate of change of the parallel component of magnetic flux -- rate of change of scalar, so frame doesn't matter. 
        S = []
        rate = []
        omega_in = body_to_inertial(omega, quaternion)
        # print('omega dh', omega_in, 'quat dh', quaternion, 'H in', H_in )

        for i in range(len(self.hysteresisOrientation)):
            hyst_orientation_in = body_to_inertial(self.hysteresisOrientation[i], quaternion)
            S.append( np.cross(omega_in, hyst_orientation_in) )
            # print('hyst orientation rot', self.hysteresisOrientation[i])
            # print('hyst orientation in', body_to_inertial(self.hysteresisOrientation[i], quaternion))
            # print('cross product', S[i])
            rate.append(np.dot(H_dot_in, hyst_orientation_in) + np.dot(H_in, S[i])) 
        # print('rate', rate)
        return rate #returns a list

    def __magneticInduction(self, H_parallel, H_dot_parallel):
        #calculate the magnetic induction based on H_body, dH_||/dt etc.. for each hysteresis material provided.
        A = []
        B = []
        C = []
        C1 = []
        C2 = []
        induction = []
        increasing_loop = []
        decreasing_loop = []
        for i in range(len(self.hysteresisOrientation)):
            A.append((2 * self.hysteresisSaturation[i]) / (np.pi))
            B.append((np.pi * self.hysteresisRemanence[i]) / (2 * self.hysteresisSaturation[i]))
            C.append((H_parallel[i] -  np.sign(H_dot_parallel[i]) * self.hysteresisCoercive[i]))
            C1.append((H_parallel[i] -  self.hysteresisCoercive[i]))
            C2.append((H_parallel[i] +  self.hysteresisCoercive[i]))
            increasing_loop.append(A[i] * np.arctan( (1 / self.hysteresisCoercive[i]) * np.tan(B[i]) * C1[i] ))
            decreasing_loop.append(A[i] * np.arctan( (1 / self.hysteresisCoercive[i]) * np.tan(B[i]) * C2[i] ))
            induction.append(A[i] * np.arctan( (1 / self.hysteresisCoercive[i]) * np.tan(B[i]) * C[i] ))
        
        return induction, decreasing_loop, increasing_loop # returns a list 


    def magnetic_induction_ODE(self, B0, H_parallel, H_dot_parallel):
        
        #ODE FUNCTION FOR DYNAMIC HYSTEReSIS RESPONSE.
        
        A = []
        B = []
        # B0 = self.B0
        B_dot = []
        hysteresisRemanence = np.array([1, 1.696, 1.696])
        sign  = np.sign(H_dot_parallel)
        B = np.append(B, (np.pi * B0 ) / (2 * self.hysteresisSaturation))
        A = np.append(A, (2 * self.hysteresisSaturation) / (hysteresisRemanence * np.pi))
        B_dot = np.append(B_dot, H_dot_parallel * A * ((( sign * H_parallel +  self.hysteresisCoercive)  * np.cos(B) - sign * hysteresisRemanence * np.sin(B) ) / (2 * self.hysteresisCoercive) )**2)

        return B_dot


    def FlatleyAndHenretty(self, B0, H_parallel, H_dot_parallel):
        sign = np.sign(H_dot_parallel)
  
        Bs = self.hysteresisSaturation
        Br = self.hysteresisRemanence
        Hc  =self.hysteresisCoercive

        H_lim = []
        # B0 = self.B0
        dB_dH = []
        B_dot = []

        k = 0.6
        p = 4.75
        q = 0.0885

        H_lim = (1/k) * np.tan((np.pi * B0) /(2 * Bs) ) - Hc
        BP = 2 * k * (Bs / np.pi) * np.cos((np.pi*B0) / (2*Bs))*2
        F = (H_parallel - H_lim) / (2 * Hc)

        for i in range(len(sign)):
            if sign[i] < 0:
                F[i] = 1 - F[i]
       
        Q = q + (1-q) * F * p
        B_dot = Q * BP * H_dot_parallel
        B_dot[0] = 0
        return B_dot

    def Flatley(self, B0, H_parallel, H_dot_parallel):
        A1 = []
        A = []
        B = []
        k = []
        H_lim = []
        sign = np.sign(H_dot_parallel)
        dB_dH = []
        B_dot = []
        Bm = self.hysteresisSaturation
        Br = self.hysteresisRemanence
        Hc = self.hysteresisCoercive
       

        A1 = (np.pi / 2) * (B0 / Bm)
        A = A1.astype(float)
        k =  (1/Hc) * np.tan((np.pi*Br)/(2*Bm))
        # print('A',A , 'type element A', type(A[0]) , type(A[1]), type(A[2]), 'type(A)', type(A))
        H_lim =  (1/k) * np.tan(A) - sign * Hc  #not - sign?
        B =  ((H_parallel - H_lim) / (2 * Hc))**2
        dB_dH =  (2/np.pi) * k * Bm * np.cos(A) * np.cos(A) * B
       
        B_dot = dB_dH * H_dot_parallel
        # print(B_dot)
        return B_dot

    def magneticInduction(self, H_parallel, H_dot_parallel, t, t0):
        
        #TO BE DEVELOPPED FOR DYNAMIC ODE MODELLING OF THE HYSTEReSIS RESPONSE.
        
        induction = []

        # #euler
        B_dots = self.magnetic_induction_ODE(H_parallel, H_dot_parallel)
        #time step
        dt = t - t0 
        if dt >= 0:
            #euler steps
            induction = B_dots * dt
            #update the previous time position
            
            self.B0 = self.B0 + induction
            # _previous = t
            # print('B', sel.f.B0)
        else:
            warnings.warn('Negative time step in magnetic induction: t: {}, t0:{}'.format(t, t0))
            
        return self.B0
        



    def perturbingAccMag(self, *args):
        self.t_permMagnet = 0
        self.t_hysteresis = 0
        H_parallel = []
        self.hysteresisMoment = []
        # split input
        B_inertial = args[0] #permanent magnet induction 
        B_inertial_dot = args[1]
        q = args[2][:4]
        omega = args[2][4:7] #in rot frame
        B_hyst = args[3]
        t = args[4]
        # t0 = args[4]
        self.totalPerturbation = 0
        # print('class_sat\n', 'q', q, 'omega', omega)
        # q = q / np.linalg.norm(q) #unitise quaternion.
        #rotate magnetic field vector into body frame
        B_env_body = inertial_to_body(B_inertial, q)
        H_body_str = B_env_body / self.mu_fs
        H_inertial_str = B_inertial / self.mu_fs
        H_inertial_dot = B_inertial_dot / self.mu_fs

        #hsyterisis #compute alignment for each rod
        for i in range(len(self.hysteresisOrientation)):
            H_parallel.append(np.dot(H_body_str, self.hysteresisOrientation[i])) # list, for each element


        #compute the rate of change of alignment and the magnetic induction for each rod
        H_dot_parallel = self.dH_parallel_dt(omega, q,  H_inertial_str, H_inertial_dot) #list, for each element
        # major, top_curve, low_curve = self.__magneticInduction(H_parallel, H_dot_parallel)
        # to delete when above line is uncommented
        major = [0, 0, 0] 
        
        B_dot = self.Flatley(B_hyst, H_parallel, H_dot_parallel)

        #ensure the computed B_dot is within the major loops:
        # for i in range(len(B_hyst[1:])):
        #     it = [top_curve[i+1], low_curve[i+1]]
        #     max_val = max(it)
        #     min_val = min(it)
        #     max_limit = max_val * self.epsilon**np.sign(max_val)
        #     min_limit = min_val * self.epsilon**(-np.sign(min_val))
        #     if B_hyst[i+1] > max_limit or B_hyst[i+1] < min_limit:
        #         pass
                # warnings.warn(' The computed B dot is :{}-- outside of the bounds:{}-{} at t:{} , q:{}'.format(B_hyst[i], min_limit, max_limit, t, q))

       
        
        #Compute the moment for each rod
        for i in range(len(self.hysteresisOrientation)):
            self.hysteresisMoment.append(self.computeMoment(B_hyst[i], self.hysteresisVolume[i], self.hysteresisOrientation[i]))
        #compute torque 
        #permanent magnet
      
 
        for i  in range(len(self.MagnetMoment)):
            self.t_permMagnet += np.cross(self.MagnetMoment[i], B_env_body)
            # print('quaternion', q, '\nBinertial', B_inertial, '\nH_body', B_env_body, '\nmoment body-fixed:', self.MagnetMoment[i], '\ntorque', self.t_permMagnet)
       
        # print('H_body', H_body)
        # print('perm magnet torque :', self.t_permMagnet)
        for i in range(len(self.hysteresisOrientation)):
            # print('rod #', i)
            # print('\nquaternion', q, 'H_inertial', H_inertial, 'H_body', H_body)
            # # print('omega', omega, 'omega in', body_to_inertial(omega, q))
            # print('H_parallel', H_parallel[i], 'H_parallel_dot', H_dot_parallel[i])
            # print('magnetic induction', B_hyst, 'Moment', self.hysteresisMoment[i] )
            # print('t', t)
        
            self.t_hysteresis += np.cross(self.hysteresisMoment[i], B_env_body)
            # print('torque produced',np.cross(self.hysteresisMoment[i], B_env_body))

        #compute gravity gradient torque:

        t_gradient =0# self.gravity_gradient_torque(q)
        #return torque
        self.totalPerturbation = self.t_permMagnet + self.t_hysteresis + t_gradient
        
        
        #write data
        if len(self.hysteresisCoercive) != 1:
            # print('writing')
            f = open('run{}/hyst_loop'.format(self.it) + '.txt', 'a')
            f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(t, H_parallel[1], B_hyst[1],H_parallel[2], B_hyst[2], major[1], major[2], np.linalg.norm(B_inertial)))
            f.close()

        # f = open('run{}/torques'.format(self.it) + '.txt', 'a')
        # f.write('{}\t{}\t{}\t{}\n'.format(t,  np.linalg.norm(self.t_permMagnet), np.linalg.norm(self.t_hysteresis), np.linalg.norm(t_gradient)))
        # f.close()
        # print('in class_sat:\n', 'total hyst', self.t_hysteresis, '\ntotal mT', self.t_permMagnet, '\ngradient', t_gradient, '\ntotal', self.totalPerturbation )
        return self.totalPerturbation, B_dot

'''
test code
'''
if __name__ == '__main__':
    r_vec = np.array([0, 6767.01, 0])
    v_vec = np.array([4.7681, 0, 6.0244])
    # kep = np.array([6378, 0.1, 58, 90, 58, 20])
    cubesat = Satellite()

    print(cubesat.position, cubesat.velocity)
    print(cubesat.eci_velocity())