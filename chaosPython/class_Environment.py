"""


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_Environment.py




This Python module defines the `ionDensity` class, which calculates the total ion density 
(H+ and O+) at the spacecraft's position based on an altitude-dependent model.

**Key Functionalities:**

-
**Assumptions and Considerations:**

- The ion density model assumes an altitude-dependent relationship for H+ and O+ ions.
- The specific implementation of the interpolation functions (`self.interpHions` 
  and `self.interpOions`) is defined in the `EnvironmentModelling` module.

"""


import numpy as np
from .EnvironmentModelling import rhoIonsFunc






class Environment:
    """
    Calculates the ion density at the spacecraft's position.

    This function calculates the total ion density (H+ and O+) at the 
    spacecraft's position (`vec_r`) based on an altitude-dependent model.
    """


    ####################
    #INSTANTIATION
    ####################
    def __init__(self):

        #create interpolation function
        self.interpOions, self.interpHions = rhoIonsFunc()
        
        #Earth radius
        self.R_E = 6378.1363




    ####################
    # METHODS
    ###################


    def ionDensity(self, vec_r):

        #compute altitude
        h = np.linalg.norm(vec_r) - self.R_E

        #compute ion density
        H_density = self.interpHions(h)*1e6             #from cm3 to m3
        O_density = self.interpOions(h)*1e6

        #return total ion density
        return H_density + O_density
