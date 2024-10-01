# Open source Coolprop module
import CoolProp.CoolProp as CP

# Simpson's rule for integration with 2nd order accuracy
from scipy.integrate import simpson

# Check potential nans
import numpy as np


class Cryogen:
    """ Class which contains a cryogen thermodynamic
    and thermophysical properties """

    R = 8.314  # Ideal gas constant

    def __init__(
        self,
        name="EmptyCryogen",
        P=0,
        T_sat=0,
        T_V=0,
        rho_L=0,
        rho_V=0,
        h_L=0,
        h_V=0,
        k_V=0,
        cp_V=0,
        MW=0,
    ):
        """Constructor"""
        self.name = name
        self.P = P  # Pressure / Pa
        self.T_sat = T_sat  # Saturation temperature / K
        self.T_V = T_sat  # Saturation temperature / K
        self.Tv_avg = T_sat  # Saturation temperature / K
        self.rho_L= rho_L  # Liquid Density / kg*m^-3
        self.rho_V = rho_V  # Vapour density / kg*m^-3
        self.rho_V_sat = rho_V  # Initialize vapour density at the interface
        self.h_L = h_L  # Liquid enthalpy J/kg
        self.h_V = h_V  # Vapour enthalpy J/kg
        self.k_V = k_V  # Thermal conductivity of the vapour W/mK
        self.k_int = k_V  # Thermal conductivity at the vapour-liquid interface
        self.cp_V = cp_V  # Heat capacity at constant pressure / J/kg/K
        self.MW = MW

    def rho_ig(self, T=None, P=None):
        """Returns ideal gas density in kg/m^3"""
        if P is None:  # Don't ask for P on isobaric conditions
            P = self.P
        if T is None:
            T = self.T_sat  # If T is not provided, use saturation temperature
        # Convert R (J mol^-1 K^-1) to R_gas (J kg^-1 K^-1)
        r_gas_si = self.R / self.MW * 1e3
        # Element-wise division work for numpy arrays
        rho_ig = P / (r_gas_si * T)
        return rho_ig
    
    def set_coolprops(self, p=100000):
        """
        Set thermophysical and thermodynamic
        properties calling CoolProp Helmholtz-based
        equation of state
        """


        fluid = self.name

        self.P = p # Pressure / Pa
        self.T_sat = CP.PropsSI('T','P',p,'Q',1,fluid) # Saturation temperature / K

        # The initial vapour temperature is at equilibrium with the liquid
        self.T_V = self.T_sat                          # Initialise vapour temperature profile / K
        self.Tv_avg = self.T_sat                      # Vapour temperature / K
        
        self.rho_L = CP.PropsSI('D','P',p,'Q',0,fluid) # Liquid mass density
        self.rho_V = CP.PropsSI('D','P',p,'Q',1,fluid) # Vapour mass density
        self.rho_V_avg = self.rho_V
        self.rho_V_sat = self.rho_V

        self.h_L = CP.PropsSI('H','P',p,'Q',0,fluid)  # Liquid enthalpy J/kg
        self.h_V = CP.PropsSI('H','P',p,'Q',1,fluid)  # Vapour enthalpy J/kg
        self.k_V = CP.PropsSI('L','P',p,'Q',1,fluid)  # Thermal conductivity of the vapour W/mK
        self.k_int = self.k_V  # Thermal conductivity at the vapour-liquid interface
        self.k_V_avg = self.k_V # Initialise average thermal conductivity
        
        self.cp_V = CP.PropsSI('C','P',p,'Q',1,fluid)  # Heat capacity at constant pressure / J/kg/K
        self.cp_L = CP.PropsSI('C','P',p,'Q',0,fluid)  # Heat capacity at constant pressure / J/kg/K
        self.cp_V_avg = self.cp_V # Initialise cp_avg
        self.MW = MW = CP.PropsSI(fluid,'molemass')
    
    def update_rho_V(self, z_grid, T_V):
        '''
        Update vapour density

        z_grid (np.array): array with the computational grid
        for the vertical coordinate

        T_v (np.array): array with the vapour temperature

        l_V: vapour length / m
        '''
        # Shift temperature 1e-3 K to avoid CoolProp non convergence
        T_V_shift = np.copy(T_V)
        T_V_shift[0] = T_V_shift[0] + 1e-3

        # Compute vapour density field
        rho_V = CP.PropsSI('D','P', self.P,'T',T_V_shift, self.name)

        # Update properties if calling CoolProp was successful
        if np.any(np.isnan(rho_V)) or np.any(np.isinf(rho_V)):
            return
        
        # Compute average density on a unit-length grid
        self.rho_V_avg =  simpson(rho_V, x = z_grid)
        return
    
    def update_k_V(self, z_grid, T_V):
        '''
        Update vapour thermal conductivity

        z_grid (np.array): array with the computational grid
        for the vertical coordinate

        T_v (np.array): array with the vapour temperature

        l_V: vapour length / m
        '''

        # Shift temperature 1e-3 K to avoid CoolProp non convergence
        T_V_shift = np.copy(T_V)
        T_V_shift[0] = T_V_shift[0] + 1e-3

        k_V = CP.PropsSI('L','P',self.P,'T',T_V_shift, self.name)

        # Update properties  if calling CoolProp was successful
        if np.any(np.isnan(k_V)) or np.any(np.isinf(k_V)):
            return

        # Update average vapour density
        self.k_V_avg = simpson(k_V, x = z_grid)
    
    def update_cp_V(self, z_grid, T_V):
        '''
        Update vapour specific heat capacity

        z_grid (np.array): array with the computational grid
        for the vertical coordinate

        T_v (np.array): array with the vapour temperature

        l_V: vapour length / m
        '''
        
        # Shift temperature 1e-3 K to avoid CoolProp non convergence
        T_V_shift = np.copy(T_V)
        T_V_shift[0] = T_V_shift[0] + 1e-3

        # Compute vapour specific heat capacity field
        cp_V = CP.PropsSI('C','P',self.P,'T',T_V_shift, self.name)

        # Update average vapour density if calling CoolProp was successful
        if np.any(np.isnan(cp_V)) or np.any(np.isinf(cp_V)):
            return
        self.cp_V_avg = simpson(cp_V, x = z_grid)

