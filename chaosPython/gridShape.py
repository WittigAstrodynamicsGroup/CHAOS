#python code holding the shape functions


'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


gridShape.py


This Python file defines functions for generating various configurations of grids 
on a CubeSat surface. The grids are represented as dictionaries, where keys are 
string representations of pixel positions (row, column), and values are lists 
containing the (x, y, z) coordinates of the pixel (in meters).

**Assumptions:**

* The CubeSat is a cubic shape with a known side length.
* The sensor pixels have a fixed radius size.

**Functions:**

* `uniformGrid(size, y_spacing=None, z_spacing=None)`: Generates a uniform grid.
* `AA_FC_grid()`: Generates a grid with a central Full Coverage (FC) region.
* `removePix(grid_dict, pix_list)`: Removes pixels from a grid based on a list of keys.
* `AA_2FC_grid()`: Generates a grid with an FC region and removes specific quadrants.
* `findCenter(r_dict)`: Calculates the center position of all pixels in a grid.
* `oppositeFace(grid_dict)`: Creates a new dictionary with negated positions (opposite face).

**Dependencies:**

* `numpy` for numerical operations
'''

import numpy as np





def uniformGrid(size, y_spacing=None, z_spacing=None):

    """
    Generates a dictionary representing a uniform grid of pixels on a CubeSat face.

    Args:
        size (int): The size (number of pixels) of the grid along one edge.
        y_spacing (float, optional): The spacing between pixels in the y-edge of the face. 
            If None, it's calculated based on sideLength and size.
        z_spacing (float, optional): The spacing between pixels in the z-edge of the face. 
            If None, it's calculated based on sideLength and size.

    Returns:
        dict: A dictionary where keys are string representations of pixel positions (row, column) 
              and values are lists containing the (x, y, z) coordinates of the pixel (meters).
    """


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




            

def AA_FC_grid():
    """
    Generates a grid with a central Faraday Cup (FC) and removes pixels overlapping with the faraday cup.

    This function represents a specific thruster configuration on a CubeSat face.
    It represents the configuration of CUbe-de-ALPS when on only one face.
    Returns:
        dict: A dictionary representing the pixel positions after removing unwanted pixels.
    """
    r_dict = uniformGrid(10)
    FC_radius = 15e-3 + 4.5e-3 # 3cm +pixel radius
    to_del = []
    for pix in r_dict.keys():
        if np.linalg.norm(r_dict[pix][1:]) <= FC_radius:
           
            to_del.append(pix)
        
    for item in to_del:
        del r_dict[item]

    return r_dict






def removePix(grid_dict, pix_list):

    """
    Removes pixels from a grid dictionary based on a list of pixel names (keys).

    Args:
        grid_dict (dict): A dictionary representing a grid of pixels, where keys are string 
                          representations of pixel positions and values are lists containing 
                          the (x, y, z) coordinates of the pixel.
        pix_list (list): A list of pixel names (keys) to be removed from the grid.

    Returns:
        dict: A new dictionary with the specified pixels removed.
    """

    for i in pix_list:
        del grid_dict[i]

    return grid_dict






def AA_2FC_grid():
    """
    Generates a grid with a central Faraday Cup (FC) region and removes specific quadrant pixels.

    This function represents a sensor configuration with an FC region on a CubeSat face 
    after removing pixels from designated quadrants (Q1, Q2, Q3, Q4).
    It represents the split of Cube-de-ALPS over 2 faces.

    Returns:
        dict: A dictionary representing the pixel positions after removing unwanted pixels 
              from specific quadrants.
    """
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
    """
    Calculates the average (center) position of all pixels in a grid.

    Args:
        r_dict (dict): A dictionary representing a grid of pixels, where keys are strings 
                          and values are lists containing the (x, y, z) coordinates of the pixel.

    Returns:
        numpy.array: A NumPy array containing the calculated center coordinates (x, y, z) in meters.
    """
    x_m = [k[0] for k in r_dict.values()]
    y_m = [k[1] for k in r_dict.values()]
    z_m = [k[2] for k in r_dict.values()]
    mean = np.array([np.mean(x_m), np.mean(y_m), np.mean(z_m)])

    return mean


def oppositeFace(grid_dict):
    """
    Creates a new dictionary with pixel positions negated, representing the opposite face of a CubeSat.
        i.e:     x -> (-x)


    Args:
        grid_dict (dict): A dictionary representing a grid of pixels on a CubeSat face.

    Returns:
        dict: A new dictionary with the same pixel names (keys) but negated positions (values).
    """
    new_dict = {}

    for pix in grid_dict.keys():
        # print(grid_dict[pix], type(grid_dict[pix]))
        pos = - np.append([], grid_dict[pix])
        new_dict[pix] =  list(pos)

    return new_dict



