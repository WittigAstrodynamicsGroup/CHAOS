#python code holding the shape functions


'''
Python code to define the shape functions. Ideally, they have as input only an information digit (grid size, line size, total amount of pixels, etc...)
and the shape function automatically generates the position dict. 
some shape functions are:

-uniformGrid

-assymetricGrid

-???
NOTE: what if the different shapes were subclasses of the main grid, which would hold the margin data. each subclass would then have the position 
generated differently, with everything else in common

'''

import numpy as np
from .quaternions import inertial_to_body




def uniformGrid(size, y_spacing=None, z_spacing=None):
    margin_y = 7.5e-3 # m // separation between edge and frst pixel
    margin_z = 7.5e-3 # m 3mm + 4.5 to centre
    sideLength =  10e-2 #m length of the side of a CubeSat

    rowLength = size
    #creates a dict holding the posiiton of pixel for a uniform grid
    if y_spacing == None:
        y_spacing = (sideLength - 2 * margin_y) / (rowLength - 1) #horizontal space between the pixels
    if z_spacing == None:
        z_spacing = (sideLength - 2 * margin_z) / (rowLength - 1) #vertical space between pixels

    r_dict = {}
    for i in range(rowLength):
        for k in range(rowLength):
            rx = 5e-2
            ry = ((margin_y + k * y_spacing) - sideLength/2)
            rz = ((margin_z + i * z_spacing) - sideLength/2)
            r = [rx, ry, rz] #r position vector
            r_dict["'{} {}'".format(k - 5,i - 5)] = r
    
    return r_dict



def isolateQuadrant(grid_dict):
    
    
    # grid_dict = uniformGrid(side)
    
    r_dict = {}
    for i in range(len(grid_dict.values())):
        pos = list(grid_dict.values())[i]
        if  pos[1] > 0 and pos[1] * pos[2] > 0:
            
            r_dict[list(grid_dict.keys())[i]] = pos

    return r_dict
            

def AA_FC_grid():

    r_dict = uniformGrid(10)
    FC_radius = 15e-3 + 4.5e-3 # 3cm +pixel radius
    to_del = []
    for pix in r_dict.keys():
        if np.linalg.norm(r_dict[pix][1:]) <= FC_radius:
           
            to_del.append(pix)
        
    for item in to_del:
        del r_dict[item]

    return r_dict

def changeFace(r_dict, angle, vector):
    angle = angle*np.pi/180
    q_rot = np.append(vector * np.sin(angle/2), np.cos(angle/2))

    new_dict = {k:inertial_to_body(v, q_rot) for k, v in r_dict.items()}
    return new_dict

def removePix(grid_dict, pix_list):

    for i in pix_list:
        del grid_dict[i]

    return grid_dict



# def applySymmetry(grid_dict):


def AA_2FC_grid():
    #call normal FC grid
    grid_dict = AA_FC_grid()

    #pixels to be removed
    Q1 = ["'2 0'",  "'4 0'", "'3 1'", "'1 1'", "'1 2'", "'3 2'", "'0 3'", "'2 3'", "'4 3'", "'3 4'", "'1 4'"]
    Q2 = ["'-3 0'",  "'-5 0'", "'-4 1'", "'-2 1'", "'-2 2'", "'-4 2'", "'-1 3'", "'-3 3'", "'-5 3'", "'-4 4'", "'-2 4'"]
    Q4 = ["'-3 -1'",  "'-5 -1'", "'-4 -2'", "'-2 -2'", "'-2 -3'", "'-4 -3'", "'-1 -4'", "'-3 -4'", "'-5 -4'", "'-4 -5'", "'-2 -5'"]  
    Q3 = ["'2 -1'",  "'4 -1'", "'3 -2'", "'1 -2'", "'1 -3'", "'3 -3'", "'0 -4'", "'2 -4'", "'4 -4'", "'3 -5'", "'1 -5'"]  
    r_dict = removePix(grid_dict.copy(), Q1)
    r_dict = removePix(r_dict, Q2)
    r_dict = removePix(r_dict, Q3)
    r_dict = removePix(r_dict, Q4)
    return r_dict

def findCenter(r_dict):
    x_m = [k[0] for k in r_dict.values()]
    y_m = [k[1] for k in r_dict.values()]
    z_m = [k[2] for k in r_dict.values()]
    mean = np.array([np.mean(x_m), np.mean(y_m), np.mean(z_m)])

    return mean


def oppositeFace(grid_dict):
    new_dict = {}

    for pix in grid_dict.keys():
        # print(grid_dict[pix], type(grid_dict[pix]))
        pos = - np.append([], grid_dict[pix])
        new_dict[pix] =  list(pos)

    return new_dict



def assignQuadrant(grid_dict_A, grid_dict_B, thrust_axis):
    ##deprecated
    new_dict = {}

    for pix in grid_dict_A.keys():

        org_torque = np.cross(grid_dict_A[pix], thrust_axis) #compute torque for given pixel on face A 

        for pixB in grid_dict_B.keys():
            
            assess_torque = np.cross(grid_dict_B[pixB], -thrust_axis) #compute torque for given pixel on face B
        
            if np.allclose(org_torque, assess_torque): #if the torques are equal, 

                new_dict[pix] = grid_dict2_B[pixB]  #then assign the pixel the same name as on face A

    return new_dict





if __name__ =='__main__':
    import matplotlib.pyplot as plt
    import scipy.optimize 
    from scipy.optimize import NonlinearConstraint

    
    spacing = None #9e-3 + 0.35e-3
    # grid_dict = uniformGrid(10, y_spacing=spacing, z_spacing=spacing)
    grid_dict1 = AA_FC_grid()
    grid_dict2 = AA_2FC_grid()
    grid_dict2_B = changeFace(grid_dict2, 180, np.array([0, 1, 0]))
    r_dict = grid_dict2
    # print('OLD', r_dict.values())




    grid_dict_B = assignQuadrant(grid_dict2, grid_dict2_B, np.array([1, 0, 0]))
    grid_dictC = oppositeFace(grid_dict2)

    print(grid_dict_B,'\n', grid_dictC   )
    # exit()
    X = []
    Y = []
    X_grid = []
    Y_grid = []
    for i in r_dict.values():

        X.append(i[1])
        Y.append(i[2])

    for i in grid_dict2_B.values():
        X_grid.append(i[1])
        Y_grid.append(i[2])

    grid_dict2_B = isolateQuadrant(grid_dict2_B)
    x_m = [k[0] for k in grid_dict2_B.values()]
    y_m = [k[1] for k in grid_dict2_B.values()]
    z_m = [k[2] for k in grid_dict2_B.values()]
    mean = np.array([np.mean(x_m), np.mean(y_m), np.mean(z_m)])
    print('mean', mean)


    # pixPos = np.zeros((len(r_dict.values()), 3))

    # for i in range(len(r_dict.values())):
    #     vec = list(r_dict.values())[i]
    #     for j in range(len(vec)):
    #         pixPos[i][j] = vec[j]
    

    # D2pos = pixPos[:, 1:]
    # print(constraints_border(D2pos, 0.0075))
    # exit()

    


    borderSpace = 3
    pixelSpace = 9e-3
    quadrantBorder = 0
    minX = quadrantBorder
    maxX = 5e-2 - borderSpace

    # 02 04, 
    



    # print(r_dict.keys())
    plt.plot([-5e-2 + 3e-3]*6, np.linspace(-5e-2 + 3e-3, 5e-2 -3e-3, 6), c='black')
    plt.plot([5e-2 - 3e-3]*6, np.linspace(-5e-2 +3e-3, 5e-2 - 3e-3, 6), c='black')
    plt.plot(np.linspace(-5e-2 + 3e-3, 5e-2 - 3e-3, 6), [5e-2 - 3e-3]*6, c='black')
    plt.plot(np.linspace(-5e-2 + 3e-3, 5e-2 - 3e-3, 6), [-5e-2 + 3e-3]*6, c='black')

    plt.plot([quadrantBorder]*6, np.linspace(quadrantBorder, 5e-2-7.5e-3, 6), c='blue', linestyle='dashed')
    plt.plot([5e-2-7.5e-3]*6, np.linspace(quadrantBorder, 5e-2-7.5e-3, 6), c='blue', linestyle='dashed')
    plt.plot(np.linspace(quadrantBorder, 5e-2-7.5e-3, 6), [5e-2-7.5e-3]*6, c='blue', linestyle='dashed')
    plt.plot(np.linspace(quadrantBorder, 5e-2-7.5e-3, 6), [quadrantBorder]*6, c='blue', linestyle='dashed')

    plt.scatter(X, Y)
    # # plt.scatter(sol.x[0], sol.x[1], marker='x', c='green')
    # # plt.scatter(otherPixel[0:, 0], otherPixel[0:, 1], marker='x', c='green')
    # plt.scatter(X_grid, Y_grid, marker='x')


    # plt.scatter(mean[1], mean[2])
    plt.ylim((-5e-2, 5e-2))
    plt.xlim((-5e-2, 5e-2))
    plt.show()

    


