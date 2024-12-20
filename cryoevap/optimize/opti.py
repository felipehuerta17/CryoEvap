import numpy as np

from ..storage_tanks.tank import Tank

from scipy.optimize import Bounds, minimize

# Visualisation
import matplotlib.pyplot as plt

class Opti:
    """ Class to be used to optimize values of Tank class"""
    def __init__(self, Tank, time = 720, dz = 0.1, thickness = 0.02):
        self.tank = Tank                # Tank class to optimize
        self.time = time * 3600         # [h] Time of the simulation
        self.dz = dz                    # [m] Vertical spacing
        self.thick = thickness          # % of the difference between internal diameter and external diameter.
        self.a_opt = 1                      # initial value of the aspect ratio.
        self.bounds = Bounds([0.1], [5])  # Minimum and maximum practical ranges
                                        # of the aspect ratio
        #self.bor
    ## thing to happen in order to optimize
    ## have a defined tank
    ## have a defined HeatTransfers parameters
    ## have a defined cryogen
    ## have a defined parameters of simulation:
    ## vertical spacing, time. 
    ## could add a function when initialize raising an error 
    ## in case any of these conditions aren't met

    def __grid(self):
        """ update the z_grid in the class tank """
        # Calculate number of nodes
        n_z = 1 + int(np.round(self.tank.l_V/self.dz, 0))
        # Define dimensionless computational grid
        self.tank.z_grid = np.linspace(0, 1, n_z)
        return
        
    def __ofunction_BOR(self, a):
        """ Objective function to minimize the Boil-off rate value
            with respect of the aspect ratio """
        # update aspect ratio
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
    
    def aspect(self):
        """ Optimize the Aspect Ratio of the tank minimizing the BOR """
        self.res = minimize(self.__ofunction_BOR, self.a_opt, method='trust-constr',tol=1e-8 , options={'verbose': 1}, bounds=self.bounds)
        self.BOR_opt = self.tank.BOR()
        self.a_opt = self.res.x[0]
        return self.a_opt
    
    def BOR_array(self, aarray):
        """ Returns an array of the BOR given an array of Aspect ratios"""
        self.BOR = []
        for a in aarray:
            self.BOR.append(self.__ofunction_BOR(a))
        return self.BOR
    
    def __BOR_plot(self, aarray):
        """ Private BOR array function to plot """
        self.__BOR = []
        for a in aarray:
            self.__BOR.append(self.__ofunction_BOR(a))
        return self.__BOR
    def rs_plot(self):
        """ Plot the Response Surface """
        # Create a plot
        self.__aspectarray = np.linspace(0.2,2,60)
        self.__BOR_plot(self.__aspectarray)
        plt.plot(self.__aspectarray, self.__BOR,"-o")

        # Add labels and title
        plt.xlabel('Aspect Ratio | a m/m')
        plt.ylabel('Boil-off Ratio | %/day')
        plt.title(f'Tank Volume: {self.tank.V} m^3, Compound: {self.tank.cryogen.name}')

        return plt.show()
        
