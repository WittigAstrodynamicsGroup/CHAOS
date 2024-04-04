#Script to hold the environment models such as solar activity, atmospehric temperature, and and Sentman models
from .transformations import julianDate
from .quaternions import *
import numpy as np 
import scipy
from scipy.special import erf
from scipy.interpolate import interp1d
from .support_fun import floatPower


def readSpaceWeatherFile(filename):
    """
    Reads a space weather file and extracts F10 solar radio flux, geomagnetic indices, and dates.

    This function assumes the space weather file format is as specified by GMAT
    NOTE: The file is modified to contain no space or non numerical characters.

    """


    yy=[]
    mm=[]
    dd=[]
    B=[]
    n=[]
    K=[]
    k2=[]
    k3=[]
    k4=[]
    k5=[]
    k6=[]
    k7=[]
    k8=[]
    tot_k=[]
    a=[]
    a2=[]
    a3=[]
    a4=[]
    a5=[]
    a6=[]
    a7=[]
    a8=[]
    avg_a=[]
    cp=[]
    isn=[]
    f10=[]
    q=[]
    ctr81=[]
    lst81=[]
    f10_obs=[]
    ctr812=[]
    lst812=[]

    date = {}
    i = 0
    datafile = open(filename, 'r')
    for line in datafile.readlines(): #split each line
        item = line.split(' ') #split each element
        
        item = [i for i in item if i != '']

        f10.append(float(item[26]))
        ctr81.append(float(item[28]))
        tot_k.append(float(item[13])/80)
        date[i] = [float(item[2]), float(item[1]), float(item[0])]
        i+=1

    return  f10, ctr81, tot_k, date



def interpF10Index_ext(filename):
    """
    Reads and interpolates F10 and geomagnetic indices from a space weather file, 
    and retruns interpolating function.

    Returns:
        tuple: A tuple containing two functions:
            - funcF10: Function to interpolate F10 index at a given time.
            - funcCtr81: Function to interpolate geomagnetic index at a given time.
    """


    #read the data file
    f10, ctr81, Kp, date = readSpaceWeatherFile(filename)
    JDtime = []


    #make the dte readable for the interpolation
    for i in range(len(date.values())):
        JDtime.append(julianDate(0, date[i][0], date[i][1], date[i][2]))


    #interpolate between the data points
    f10InterpCubic = scipy.interpolate.interp1d(JDtime, f10, kind='cubic')
    ctr81InterpCubic = scipy.interpolate.interp1d(JDtime, ctr81, kind='cubic')


    #create the interpolation functions
    def funcF10(t):

        #for times provided by the datafile
        if JDtime[0] <= t <= JDtime[-1]:
            f = f10InterpCubic(t)
        
        #if time is beyond what is provided in datafile, switch to scaled sinusoidal motion
        else:
            f = 145 + 75 * np.cos(0.001696*t + 0.35 * np.sin(0.00001695*t))
        
        return f

    def funcCtr81(t):
        if JDtime[0] <= t <= JDtime[-1]:
            f = ctr81InterpCubic(t)
        else:
            f =  145 + 75 * np.cos(0.001696*t + 0.35 * np.sin(0.00001695*t))
        return f

    return funcF10, funcCtr81

def atmTemperature(f10, f10_avg, vec_r_sun, vec_r_sat, EarthFlattening):
    """
    Calculates the uncorrected exospheric temperature of the atmosphere.

    This function considers solar activity and the Sun's Local Hour Angle (LHA)
    to estimate the uncorrected exospheric temperature.

    Args:
        f10 (float): 10.7 cm solar radio flux index.
        f10_avg (float): Average 10.7 cm solar radio flux index.
        vec_r_sun (np.ndarray): Sun's position vector in satellite reference frame (x, y, z).
        vec_r_sat (np.ndarray): Satellite's position vector in Earth-centered inertial frame (ECI) (x, y, z).
        EarthFlattening (float): Earth's flattening factor.

    Returns:
        float: Uncorrected exospheric temperature in Kelvin.
    """

    #nighttime global exospheric temperature
    Tc = 379 + 3.24 * f10_avg + 1.3 * (f10 - f10_avg) #Kelvin

    # Uncorrected exospheric temperature

    #Positions
    rs_x, rs_y, rs_z = vec_r_sun
    ri, rj, rk = vec_r_sat

    decl = np.arccos( (np.sqrt(rs_x**2 + rs_y**2)) / np.linalg.norm(vec_r_sun))
    phi = np.arctan( (1 /(1 - EarthFlattening)**2) * (rk / (np.sqrt(ri**2 + rj**2))) )
   


    inside = ((rs_x * ri) + (rs_y * rj)) / (np.sqrt(rs_x**2 + rs_y**2) * np.sqrt(ri**2 + rj**2))
    
    #clamping to avoid floating point issues.
    if inside <=-1:
        inside = -1

    elif inside >=1:
        inside = 1
        
    elif -1 < inside < 1:
        pass #inside = inside

    #Local Hour Angle
    LHA = ((rs_x * rj - rs_y * ri) / abs(rs_x * rj - rs_y * ri)) * np.arccos(inside)
    


    LHA_deg = LHA * 180/ np.pi

    eta = abs(phi - decl) / 2
    theta = abs(phi + decl) / 2

    tau_deg = LHA_deg - 37 + 6 * np.sin(LHA + (43 * np.pi/180))
    tau = tau_deg * np.pi / 180
    
    #compute real solution of exponent 2.2 
    A1 = np.sin(theta)
    B1 = np.cos(eta)
    C = np.cos(tau/2)

    A = floatPower(A1, 2.2)
    B = floatPower(B1, 2.2)

    #compute uncorrected exospheric temperature
    T_unc = Tc * ( 1 + 0.3 * (A + (B - A) * C**3 )  )

    return T_unc

    

def MoeSentman(arr_angle, arr_A, Aref, Ta, arr_Tw, Vorb, arr_alpha):
    """
    Calculates drag (Cd) and lift (Cl) coefficients for a spacecraft surface
    using the Moe-Sentman model.

    This function implements the Moe-Sentman model to estimate the drag (Cd)
    and lift (Cl) coefficients experienced by a spacecraft surface due to
    reflected molecules. It takes various parameters as input:

    Args:
        arr_angle (np.ndarray): Array of incidence angles (radians).
        arr_A (np.ndarray): Array of accommodation coefficients (0 to 1).
        Aref (float): Reference area of the spacecraft surface (m^2).
        Ta (float): Ambient temperature (K).
        arr_Tw (np.ndarray): Array of wall temperatures of the surface (K).
        Vorb (float): Orbital velocity of the spacecraft (m/s).
        arr_alpha (np.ndarray): Array of thermal accommodation coefficients.

    Returns:
        tuple: A tuple containing two NumPy arrays:
            - arr_Cd (np.ndarray): Array of drag coefficients.
            - arr_Cl (np.ndarray): Array of lift coefficients.
    """

    #set universal constants:
    R = 8.3145              # gas constant
    Ma = 28.964
    
    
    #define parameters
    s = Vorb / (2 * R * Ta / Ma)**0.5
    Ti = Ma * (Vorb**2) / (3 * R)
    arr_gamma = np.cos(arr_angle)
    arr_l = np.sin(arr_angle)

    
    
    #compute helpers
    arr_P = np.exp(-(arr_gamma**2) * (s**2)) / s
    G = ( 1 / (2 * (s**2)) )
    Q = 1 + ( 1 / (2 * (s**2)) )
    arr_Z = 1 + erf(arr_gamma*s)


    #compute Vr, velocity of reflected particles
    Vr = (2/3)**(0.5) * Vorb * (1 + arr_alpha * ((arr_Tw/Ti) - 1))**0.5


    #Cd
    arr_Cd = (arr_A / Aref) * ( (arr_P/np.sqrt(np.pi)) + (arr_gamma*Q*arr_Z) + ( (arr_gamma * Vr)/(2 * Vorb) )  * (arr_gamma*np.sqrt(np.pi)*arr_Z + arr_P) )
    #Cl
    arr_Cl = (arr_A / Aref) * ( (arr_l*G*arr_Z) + ( (arr_l * Vr)/(2*Vorb) ) * (arr_gamma*np.sqrt(np.pi)*arr_Z + arr_P) )

    return arr_Cd, arr_Cl





def rhoIonsFunc():
    """
    Provides functions for interpolating ion densities based on altitude.

    This function defines two interpolation functions, `funcAO` and `funcH`,
    which estimate ion densities (O+ and H+) for a given altitude (h).

    It uses pre-defined altitude (h) and ion density (O+, H+) data points stored in lists.
    The function employs linear interpolation to estimate ion densities for altitudes
    that fall within the range of the provided data points.

    Returns:
        tuple: A tuple containing two functions:
            - funcAO: Function to estimate O+ ion density for a given altitude.
            - funcH: Function to estimate H+ ion density for a given altitude.
    """




    ##AO read graph vals
    O_ions = [192, 192.39499180074446, 250.35970056089553, 320.6295377461545, 373.124101961492, 462.83830908411693, 565.0329026967139, 759.1148882722711, 1036.2697849656568, 1472.2008006608266, 2025.807400828576, 2721.6478038741084, 3928.779850797452, 5537.147982918537, 7866.478220341349, 11912.441462435127, 17333.741935660863, 26671.308122962924, 38500.82946028165, 65194.567312008454, 92620.1801769758, 111280.35596055829, 129499.55632787442, 142514.02686584866, 155589.81527961403, 167175.70507846994, 153126.23697157818, 133699.99495840858, 115810.355343675, 101928.3402998044, 91153.65059494431, 77089.05765053419, 62146.50957443905, 41699.21333456813, 27979.437696964487, 19382.652643700334, 12397.372985258204, 8318.418917498131, 5031.4921476797235, 2947.7435677687095, 1768.8049065540579, 972.1798483448839, 574.1234598078929, 358.5291035548677, 222.1149340009125, 139.81786267358868, 104.07079878925109,  70, 45, 15, 5, 0]
    O_ions_h = [133.25581395348837, 148.25581395348837, 151.1627906976744, 154.06976744186045, 156.97674418604652, 159.88372093023256, 162.7906976744186, 168.60465116279067, 171.51162790697674, 174.41860465116278, 180.2325581395349, 186.04651162790697, 188.953488372093, 191.86046511627904, 197.67441860465115, 203.48837209302326, 209.30232558139534, 212.20930232558138, 220.93023255813952, 223.83720930232556, 235.46511627906975, 241.27906976744185, 247.09302325581393, 250, 258.72093023255815, 279.06976744186045, 299.41860465116275, 316.86046511627904, 331.39534883720927, 340.1162790697674, 351.74418604651163, 366.27906976744185, 383.7209302325581, 409.8837209302325, 438.953488372093, 459.30232558139534, 488.3720930232558, 505.8139534883721, 529.0697674418604, 552.3255813953488, 572.6744186046511, 595.9302325581396, 616.2790697674418, 636.6279069767442, 656.9767441860465, 680.2325581395348, 691.860465116279, 786.8852459016393, 1034.2771982116244, 1305.5141579731744, 1737.704918032787, 2000]
    
    
    ## H read graph lists
    H_ions = [102.42296238083073, 102.42296238083073, 126.0397549865452, 157.59749964349658, 206.72155404294074, 279.9531186032767, 376.1136375169116, 530.0874440450517, 753.0810710352223, 1472.2008006608266, 2042.0385239894013, 3296.177471017853, 4645.57015870548, 6495.33690490405, 8452.250172073447, 10401.178759232058, 11446.478359946557, 12007.886026633996, 12104.095309372206, 12104.095309372206, 11538.18954628246, 10911.31832481692, 10484.51482236907, 9836.0814701533, 9451.336575278814, 9081.64123388567, 8796.324269742678, 8385.06754455816, 7993.038406809362, 7741.921962790426, 7498.694817589695, 7148.106483272828, 6868.503524939906, 6599.837423030635, 6495.33690490405,6445.33690490405]
    H_ions_h = [133 ,247.09302325581393, 252.90697674418604, 261.62790697674416, 270.3488372093023, 284.8837209302325, 302.3255813953488, 319.7674418604651, 343.0232558139535, 395.3488372093023, 415.6976744186046, 441.86046511627904, 459.30232558139534, 482.5581395348837, 508.7209302325581, 540.6976744186046, 563.953488372093, 598.8372093023255, 633.7209302325581, 677.3255813953488, 747.0930232558139, 828.4883720930233, 895.3488372093022, 970.9302325581394, 1049.4186046511627, 1119.1860465116279, 1183.139534883721, 1273.2558139534883, 1398.2558139534883, 1494.1860465116279, 1593.0232558139535, 1703.4883720930231, 1802.3255813953488, 1918.6046511627906, 1982.5581395348836, 2100]
    
    
    f_AO2 = interp1d(O_ions_h, O_ions, kind='linear')
    f_H2 = interp1d(H_ions_h, H_ions, kind='linear')


    #interpolate for O
    def funcAO(h):
        if O_ions_h[0] <= h <= O_ions_h[-1]:
            return f_AO2(h)
        else:
            return f_AO2(O_ions_h[0])
    
    #interpolate for H
    def funcH(h):
        if H_ions_h[0] <= h <= H_ions_h[-1]:
            return f_H2(h)
        else:
            return f_H2(H_ions_h[0])

    return funcAO, funcH

