# Open source Coolprop module
import CoolProp.CoolProp as CP

class Cryogen:
    """ Class which contains a cryogen thermodynamic
    and thermophysical properties """

    R = 8.314  # Ideal gas constant

    def __init__(
        self,
        name="EmptyCryogen",
        P=0,
        T_sat=0,
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
        self.rho_L = rho_L  # Liquid Density / kg*m^-3
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
    
    def set_coolprops(self, fluid ="methane", p=100000):
        """
        Set thermophysical and thermodynamic
        properties calling CoolProp Helmholtz-based
        equation of state
        """

        self.T_sat = CP.PropsSI('T','P',p,'Q',1,fluid) # Saturation temperature / K
        self.rho_L = CP.PropsSI('D','P',p,'Q',0,fluid) # Liquid mass density
        self.rho_V = CP.PropsSI('D','P',p,'Q',1,fluid) # Vapour mass density
        self.rho_V_sat = self.rho_V

        self.h_L = CP.PropsSI('H','P',p,'Q',0,fluid)  # Liquid enthalpy J/kg
        self.h_V = CP.PropsSI('H','P',p,'Q',1,fluid)  # Vapour enthalpy J/kg
        self.k_V = CP.PropsSI('L','P',p,'Q',1,fluid)  # Thermal conductivity of the vapour W/mK
        self.k_int = self.k_V  # Thermal conductivity at the vapour-liquid interface
        self.cp_V = CP.PropsSI('C','P',p,'Q',1,fluid)  # Heat capacity at constant pressure / J/kg/K
        self.MW = MW = CP.PropsSI(fluid,'molemass')








