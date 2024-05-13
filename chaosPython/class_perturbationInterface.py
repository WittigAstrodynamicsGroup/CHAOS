
"""

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_perturbationsInterface.py





This Python module defines the `perturbationInterface` class, 
which serves as the base class for all spacecraft perturbation models.

The `perturbationInterface` class enforces a common interface for 
calculating environmental forces and torques acting on a spacecraft. 

Key functionalities include:

- **Abstract Method (`computePerturb`):** 
  This method serves as a placeholder, requiring implementation in 
  derived classes. It defines the expected behavior for calculating 
  the total force and torque vectors acting on a spacecraft due to 
  environmental effects.
  
- **Enforced Implementation:** 
  The `__init_subclass__` method ensures that `perturbationInterface` 
  cannot be directly instantiated. Derived classes inheriting from 
  `perturbationInterface` must provide their own implementation for 
  `computePerturb`.
"""

class perturbationInterface:


  ####################
  # INSTANTIATION
  ####################


  # This method is automatically called when a subclass of MyInterface is created
  def __init_subclass__(cls, **kwargs):

      # Check if the subclass being created is MyInterface itself
      if cls.__name__ == 'MyInterface':

          # If so, raise a TypeError indicating that MyInterface cannot be instantiated
          raise TypeError("Can't instantiate abstract class perturbationInterface. Implement the force model in a derived class of perturbationInterface")
      
      # Call the __init_subclass__ method of the superclass (if any)
      super().__init_subclass__(**kwargs)



  ####################
  # METHODS
  ####################

  # This method is marked as abstract and must be implemented by subclasses
  def computePerturb(self, time, satelliteClass):
      """
      This function is designed to compute the environmental perturbations acting on a spacecraft.

      **Abstract Method:** This function serves as an abstract method (or a placeholder) 
      and requires implementation in derived classes. It takes three arguments:

      - `self`: Reference to the current object instance.
      - `time`: Simulation time 
      - `satelliteClass`: Reference to a class holding information about the spacecraft 
                              (e.g., mass, area, coefficients).


      NOTE: This function raises a `NotImplementedError` to enforce implementation 
      in child classes, as it's an abstract method.

      **Returns:**
          tuple(numpy.ndarray): A tuple containing two NumPy ndarrays representing:
              - forceVector (km/s^2): Total force vector acting on the spacecraft in the simulation.
              - torqueVector (Nm): Total torque vector acting on the spacecraft in the simulation.
      """
      ##OVERRIDE WITH ACTUAL IMPLEMENTATION
      raise NotImplementedError("The function [computePerturb] must be implemented, with inputs (self, time, satelliteClass) and return (forceVector, torqueVector)!")
      