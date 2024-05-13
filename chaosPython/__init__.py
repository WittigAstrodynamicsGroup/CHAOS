"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


__init__.py




This file initializes the chaosPython Python module.

It imports all the relevant classes and functions from various files within the 
module structure, making them accessible when using the module.

"""






#init file to import all the relevant classes to CHAOS
#this way all classes and functions are accessible through the module.
from .CHAOS import CHAOS
from .class_sensor import Sensor
from .class_satellite import Satellite
from .class_control_system import control_system
from .gridShape import AA_2FC_grid, AA_FC_grid, oppositeFace
from .class_perturbationManager import Environment
from .class_grid import Grid
from .class_dataset import Dataset
from .support_fun import basicQuadrant
from .transformations import eq_to_kep, equi_to_r_cart
from .smart_fun import pixel_fuel, r_stop, sensorMeasurement, assessGridState, assessGridState2
from .SPICEroutine import loadKernels, unloadKernels
from .readCHAOS import plotData
import os


#if any checks must br run, run them here?

moduleDirName = os.path.dirname(__file__)