import numpy as np

from ..storage_tanks.tank import Tank

from scipy.optimize import Bounds, minimize

# Visualisation
import matplotlib.pyplot as plt

class Opti:
    """
    Class Opti
    ------------------
    The Opti class is designed to perform optimization on a cryogenic storage tank system.
    It takes a Tank object and various simulation parameters to set up the optimization
    problem, focusing on the tank's aspect ratio and other physical properties.
    
    Parameters
    ----------
    Tank_obj : Tank
        An instance of the Tank class representing the cryogenic tank to be optimized.
        Must have all required heat transfer properties and a defined cryogen.
    time : float, optional
        Duration of the simulation in hours. Default is 720 hours.
    dz : float, optional
        Vertical spacing in meters for the simulation grid. Default is 0.1 m.
    thickness : float, optional
        Percentage of the difference between the internal and external diameter of the tank,
        representing wall thickness. Default is 0.02.
    x0 : float, optional
        Initial value for the aspect ratio to be optimized. Default is 1.
    bounds : list of float, optional
        Minimum and maximum practical values for the aspect ratio. Default is [0.1, 5].
    
    Attributes
    ----------
    tank : Tank
        The provided Tank object to be optimized.
    time : float
        Simulation time in seconds.
    dz : float
        Vertical spacing in meters.
    thick : float
        Wall thickness as a percentage.
    a_opt : float
        Initial aspect ratio value.
    bounds : Bounds
        Bounds object specifying the minimum and maximum aspect ratio.
        
    Methods
    -------
    aspect(verbose=1, tol=1e-8)
        Optimizes the aspect ratio of the tank by minimizing the boil-off rate (BOR).
        Returns the optimal aspect ratio found.
    BOR_array(a_array)
        Calculates and returns an array of boil-off rates (BOR) for a given array of aspect ratios.
    rs_plot()
        Plots the response surface of the boil-off rate (BOR) as a function of the aspect ratio.
    """

    def __init__(self, Tank_obj, time = 720, dz = 0.1, thickness = 0.02, x0 = 1, bounds = [0.1, 5]):

        self.tank   = Tank_obj                      # Tank class to optimize
        self.time   = time * 3600                   # [h] Time of the simulation
        self.dz     = dz                            # [m] Vertical spacing
        self.thick  = thickness                     # % of the difference between internal diameter and external diameter.
        self.a_opt  = x0                            # initial value of the aspect ratio.
        self.bounds = Bounds(bounds[0], bounds[1])  # Minimum and maximum practical ranges of the aspect ratio

        # Input validation
        if not isinstance(Tank_obj, Tank):
            raise TypeError("The provided object is not an instance of the Tank class. Please provide a valid Tank object.")

        if Tank_obj.cryogen.name == 'EmptyCryogen':
            raise ValueError("No cryogen defined in the Tank object. Please assign a valid cryogen to the Tank.")

        required_attrs = ['U_L', 'U_V', 'Q_roof', 'Q_b_fixed', 'T_air', 'eta_w']
        for attr in required_attrs:
            if not hasattr(Tank_obj, attr):
                raise AttributeError(f"Tank object missing required attribute '{attr}'. Please ensure all necessary heat transfer properties are defined.")

    def __grid(self):
        """ Update the z_grid in the class tank """
        # Calculate number of nodes
        n_z = 1 + int(np.round(self.tank.l_V/self.dz, 0))

        # Define dimensionless computational grid
        self.tank.z_grid = np.linspace(0, 1, n_z)
        return
        
    def __ofunction_BOR(self, a):
        """ Objective function to minimize the Boil-off rate value
            with respect of the aspect ratio """
        # Update aspect ratio
        # Internal diameter of the tank defined as using the aspect ratio
        self.tank.d_i = ((4 * self.tank.V)/(np.pi * a))**(1/3)
        
        # External diameter of the tank defined as a porcentage of the internal diameter
        self.tank.d_o = self.tank.d_i * (1 + self.thick)
        
        # Update the zgrid
        self.__grid()
        
        # Reinitialize cryogen
        self.tank.cryogen.set_coolprops(self.tank.cryogen.P)
        
        # Run the simulation
        self.tank.evaporate(self.time)
        
        # return the BOR value
        return self.tank.BOR()
    
    def aspect(self, verbose = 1, tol = 1e-8):
        """
        Optimize the aspect ratio of the tank by minimizing the boil-off rate (BOR).
        
        Inputs
        ----------
        verbose : int, optional. Controls the verbosity of the optimization output.
            - 1: No output (default)
            - 2: Display optimization process information
        tol : float, optional. Tolerance for the optimization process, default is 1e-8
        
        Return
        -------
        The optimal aspect ratio found (a_opt), float.
        """
        
        self.res     = minimize(self.__ofunction_BOR, self.a_opt, method='trust-constr',tol = tol , options={'verbose': verbose}, bounds=self.bounds)
        self.BOR_opt = self.tank.BOR()
        self.a_opt   = self.res.x[0]
        return self.a_opt
    
    def BOR_array(self, a_array):
        """
        Returns an array of the BOR given an array of Aspect ratios.
        
        Input
        ----------
        a_array : numpy array. Array of aspect ratios for which to compute the boil-off rates.

        Return
        -------
        BOR : list. List of boil-off rates corresponding to each aspect ratio in a_array.
        """

        self.BOR = []
        for a in a_array:
            self.BOR.append(self.__ofunction_BOR(a))
        return self.BOR
    
    def __BOR_plot(self, a_array):
        """ Private BOR array function to plot """
        self.__BOR = []
        for a in a_array:
            self.__BOR.append(self.__ofunction_BOR(a))
        return self.__BOR
    
    def rs_plot(self, range = [0.2, 2], n_points = 60):
        """
        Plots the response surface of the boil-off rate (BOR) as a function of the aspect ratio.

        Inputs
        ----------
        range : list of float, optional. The lower and upper bounds for the aspect ratio to plot. 
            - Default is [0.2, 2].
        n_points : int, optional. Number of points to evaluate between the range. 
            - Default is 60.

        Return
        -------
        None
            - Displays a matplotlib plot of BOR versus aspect ratio.
        """

        # Create a plot
        self.__aspectarray = np.linspace(range[0], range[1], n_points)
        self.__BOR_plot(self.__aspectarray)
        plt.plot(self.__aspectarray, self.__BOR,"-o")

        # Add labels and title
        plt.xlabel('Aspect Ratio | a m/m')
        plt.ylabel('Boil-off Ratio | %/day')
        plt.title(f'Tank Volume: {self.tank.V} m^3, Compound: {self.tank.cryogen.name}')

        return plt.show()
        
