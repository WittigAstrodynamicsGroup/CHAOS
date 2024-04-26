#python code to hold dataset class


'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


class_dataset.py


Storage class that handles storing data and writing it to a textfile.
The methods available are the following
        -write_value:   Write n different values (1 number) to a textfile. Each column is a value

'''




class Dataset:
    """
    **Dataset Class Initialization**

    This constructor initializes the `Dataset` class with no input arguments.

    **Attributes:**

    * `breaker` (int): This attribute tracks if the integration has moved in time.
    * `t_start` (float): This attribute stores the start time of the simulation (initially 0).
    * `y0` (None): This attribute is intended to hold the initial condition from a checkpointing file (initially None).
    """

    #No input argument when instantiated
    def __init__(self):
        
        self.breaker = 0                    #Breaker variable
        self.t_start = 0                    #Start time of the simulation
        self.y0 = None                      #initial condition from checkpointing file
    
    # define methods here


    
    def write_value( self, filename, *args):
        """
        **Write Values to Text File**

        This method writes n numerical values (passed as arguments) to a text file.

        **Arguments:**

        * `filename` (str): The name of the text file to write to (assumed to be already opened for writing).
        * `*args` (float): A variable number of numerical arguments representing the values to be written.

        **Functionality:**

        The method iterates through the arguments, converting them to strings separated by tabs (`\t`), and writes them to the 
        specified file followed by a newline character (`\n`).
        """       
        #NOTE: file must be open beforehand
        f = filename
        f.write(''.join('%s\t' % item for item in args))
        f.write('\n')




    

