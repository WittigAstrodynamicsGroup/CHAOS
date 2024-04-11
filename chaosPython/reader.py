'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


reader.py


This file defines two functions for reading data from tab-separated files:

- `read_data(filename)`: Reads numerical data from a file. It expects each line to contain five tab-separated floating-point values (time, x, y, z, norm).

- `read_grid_data(filename)`: Reads data from a file, likely grid-related information

'''




def read_data(filename):
    """
    Reads data from a file containing numerical values separated by tabs.

    Args:
        filename (str): Path to the data file.

    Returns:
        tuple: A tuple containing separate lists of data points.
            - t (list): List of time values (first column).
            - x (list): List of x values (second column).
            - y (list): List of y values (third column).
            - z (list): List of z values (fourth column).
            - norm (list): List of norm values (fifth column).
    """



    t=[]
    x=[]
    y=[]
    z=[]
    norm=[]
    datafile = open(filename, 'r')
    for line in datafile.readlines(): #split each line
        item = line.split('\t') #split each element
        t.append(float(item[0]))
        x.append(float(item[1]))
        y.append(float(item[2]))
        z.append(float(item[3]))
        norm.append(float(item[4]))

    return t, x, y, z, norm







def read_grid_data(filename):

    """
    Reads data from a file containing grid data (similar format to read_data).


    Args:
        filename (str): Path to the grid data file.

    Returns:
        tuple: A tuple containing separate lists of data points, likely related to grid information.
            - t (list): List of time values (first column).
            - x (list): List of x values (second column)
            - y (list): List of y values (third column)
            - z (list): List of z values (fourth column)
            - norm (list): List of norm values (fifth column)
    """



    t=[]
    x=[]
    y=[]
    z=[]
    norm=[]
    datafile = open(filename, 'r')
    for line in datafile.readlines(): #split each line
        item = line.split('\t') #split each element
        t.append(float(item[0]))
        x.append(item[1])
        y.append(float(item[2]))
        z.append(float(item[3]))
        norm.append(float(item[4]))

    return t, x, y, z, norm

