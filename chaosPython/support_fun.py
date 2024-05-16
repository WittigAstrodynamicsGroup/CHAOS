"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


support_fun.py


This file contains a collection of support functions, including those for:

- Mathematical operations
- Geometric calculations
- Grid-related operations
- Position partitioning

"""


import numpy as np 




def floatPower(number, exp):
    """
    Calculates a number raised to an exponent, preserving sign for negative exponents.
    This form avoids complex numebr when the exponent is fractional.
    Args:
        number: The number to be exponentiated.
        exp: The exponent.

    Returns:
        float: The result of the exponentiation, maintaining sign for negative exponents.
    """

    sol = 0 #default 
    if number > 0:
        sol = number ** exp

    elif number < 0:
        sol = - (abs(number)** exp)         # Preserve sign for negative bases

    elif number == 0:
        sol = 0

    return sol

    
def findCenter(r_dict):

    """
    Calculates the mean (center) position from a dictionary of positions.

    Args:
        r_dict: A dictionary containing positions as 3-element tuples (x, y, z).

    Returns:
        numpy.ndarray: The mean position as a 3-element array.
    """

    x_m = [k[0] for k in r_dict.values()]          #store all x values
    y_m = [k[1] for k in r_dict.values()]          #store all y values
    z_m = [k[2] for k in r_dict.values()]          #store all z values

    #compute the mean of each element
    mean = np.array([np.mean(x_m), np.mean(y_m), np.mean(z_m)])

    return mean


def gridsBurntime(grids):
    """
    Calculates the total burn time of all grids in a list of grids.

    Args:
        grids: A list of grid objects, each having a 'burnTime' attribute.
        burnTime must be a dictionary!
    Returns:
        float: The total burn time of all grids.
    """
    #gets the total fuel in all of the grids
    totalBurntime = 0
    for grid in grids:
        totalBurntime += sum(grid.burnTime.values())

    return totalBurntime


def selectGrid(grids):
    """
    Selects a grid to ignite from a list based on fuel levels and firing state.
    If no grids should be firing, returns a default grid
    Args:
        grids: A list of grid objects.

    Returns:
        grid: The selected grid object, based on criteria.
    """


    #check if there are any grids with pixels left with fuel
    #Criterion is the total fuel remaining on the grid 
    #averages to greater than 1 second of firing per pixel
    fuelled = [g for g in grids if sum(g.burnTime.values()) > g.numHeads]

    #For grids that have fuel remaining, check if the grid is on
    firing = [c for c in fuelled if c.state == 1]

    #if none have fuel
    if len(fuelled) == 0:
        #it doesn't matter which one, because they don't fire.
        selected_grid = grids[0]


    #if none  are firing
    elif len(firing) == 0:
        #it doesn't matter which one, because they don't fire.
        selected_grid = fuelled[-1]

    #one is firing
    elif len(firing) == 1:
        selected_grid = firing[0]
    
    #if you selected more than 1 grid (not supported currently)
    elif len(firing) > 1:
        print('LEN OF SELECTED GRID > 1!!')
        raise ValueError

    return selected_grid


def circlesIntersection(r_a, r_b, h, theta):
    """
    Calculates the area of intersection between two circles.
    This function is explicitly designed for the Faraday cup illumination area. It computes
    the intersection area of the faraday cup sensor based on the height of the cup and the 
    flow angle.
    Args:
        r_a: Radius of the first circle (faraday cup sensor).
        r_b: Radius of the second circle (faraday cup aperture).
        h: Height of the cup.
        theta: Orientation of the Faraday cup relative to the flow.


        NOTE: theta is always positive!!
        NOTE: For CHAOS, r_a = r_b
    Returns:
        float: The area of intersection between the circles.
    """



    A_ill = 0
    #centre to centre distance between the circles
    d = h * np.tan(theta)

    #If centre-to-centre distance is greater than 
    #the sum of the radii
    if d > r_a  + r_b or d < 0:
        A_ill = 0

    #If both centres are over-imposed
    elif d == 0:
        A_ill = np.pi * r_a**2
        #A_ill = np.pi * min(r_a, r_b)** 2      #If both circles are not the same size, use this     
    
    #If there is a partial overlap of the circles
    elif 0 < d <= r_a + r_b:

        #x position of intersection point
        x = ((r_a**2) - (r_b**2) + (d**2)) / (2 * d)

        #y positions of intersection points
        y1 = np.sqrt((r_a**2) - (x**2))

        #angles of sectors formed
        gamma_a = 2 * np.arcsin(y1 / r_a)
        gamma_b = 2 * np.arcsin(y1 / r_b)

        #surface of sector formed
        S_a = (r_a**2) * gamma_a / 2
        S_b = (r_b**2) * gamma_b / 2

        #area of triangle formed
        A_a = y1 * x
        A_b = y1 * (d - x)

        #total area
        A_ill = S_a - A_a + S_b - A_b

    return A_ill 


def basicQuadrant(r_dict):
    """
    Separates positions in a dictionary into four quadrants based on x-y coordinates.
    NOTE: Works only for even numbers of positions.

    Args:
        r_dict: A dictionary containing positions as 3-element numpy array (x, y, z).

    Returns:
        tuple: Four dictionaries containing positions in each quadrant (I, II, III, IV).
    """
    I = {}
    II = {}
    III = {}
    IV = {}

    for i in r_dict.keys():
        if r_dict[i][1] > 0 and r_dict[i][2] > 0: #top left
            I[i] = r_dict[i]

        elif r_dict[i][1] > 0  and r_dict[i][2] < 0: #bottom left
            IV[i] = r_dict[i]

        elif r_dict[i][1] < 0 and r_dict[i][2] > 0: #top right
            II[i] = r_dict[i]

        elif r_dict[i][1] < 0 and r_dict[i][2] < 0: #bottom right
            III[i] = r_dict[i]

    return I, II, III, IV



def updateFuelMass(grids, satelliteClass):
    """
    Calculates the mass of fuel consumption and the remaining mass.

    This function iterates through a list of grids 
    provided in `grids` and calculates the  mass of fuel 
    consumption based on a burning rate model.

   
    **Arguments:**
        - grids: List of objects representing fuel tank grids (content 
                 assumed).
        - satelliteClass: Reference to a class holding information about the 
                          spacecraft, including its initial mass (`mass`).

    **Returns:**
        float: The updated mass of the spacecraft (kg) after accounting for 
              fuel consumption.
    """
    
    #for each grid, compute parameters
    totalFuelMass = 0
    totalBurntime = 0
    burntimeVal = gridsBurntime(grids)
    for grid in grids:
        totalFuelMass += grid.fuelMass
        totalBurntime += grid.totalburnTime

    #compute  fuel consumption:
    consumed_fuel = totalFuelMass*(1 - (burntimeVal / totalBurntime))

    #compute remaining fuel
    #NOTE: Use sat.mass0 as mass0 is not updated!
    remaining_mass = satelliteClass.mass0 - consumed_fuel

    return remaining_mass