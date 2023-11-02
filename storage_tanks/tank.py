import numpy as np
# from ..cryogens.cryogen import Cryogen
from cryogens import Cryogen

class Tank:
    """ Class to be used as a container for the
    evaporation of pure cryogens"""

    def __init__(self, d_i, d_o, V, LF=0.97):
        """ Class constructor """
        # Compulsory parameters
        self.d_i = d_i  # [m] Tank internal diameter
        self.d_o = d_o  # [m] Tank external diameter
        self.V = V  # [m^3] Tank volume
        self.A_T = np.pi * d_i ** 2 / 4  # [m^2] Area of the surface
        # perpendicular to the vertical axis
        self.l = V / self.A_T  # [m] Tank height
        self.roof_BC = "Neumann"  # Roof Temperature boundary condition,
        # "Neumann" or "Robin"
        self.thermophysical_it = False  # Thermophysical iteration
        self.LF = LF
        self.cryogen = Cryogen()  # Empty Cryogen, see Cryogen class
        # switch for the non-eq model
        pass

    def set_HeatTransProps(self, U_L, U_V, T_air, Q_b_fixed=None, Q_roof=0):
        """Set separately tank heat transfer properties
        Usage: set_HeatTransProps(self, U_L, U_V, Q_b, Q_roof, T_air)"""
        self.U_L = U_L  # [W*m^-2*K^-1]Overall heat transfer coefficient
        # for the liquid phase stored in the tank
        self.U_V = U_V  # [W*m^-2*K^-1] Overall heat transfer coefficient
        # for the vapour phase stored in the tank
        # Is there any Q_b_fixed?
        self.Q_b_fixed = Q_b_fixed
        self.Q_roof = Q_roof  # [W] Heat ingress through the roof
        self.T_air = T_air  # [K] Temperature of the surrounding air /K
        pass

    @property
    def l_V(self):
        """Update liquid filling and vapour length"""
        return self.l * (1 - self.LF)  # [m] sets vapour length

    @property
    def A_L(self):
        """Tank wall area in contact with the liquid"""
        return np.pi * self.d_o * self.l * self.LF

    @property
    def Q_L_in(self):
        """ Liquid heat ingress through the walls
        in W """
        return self.U_L * self.A_L * (self.T_air - self.cryogen.T_sat)

    @property
    def v_z(self):
        """Update advective velocity with respect to tank liquid filling"""
        # Initial evaporation rate mol/s
        BL_0 = (self.Q_L_in + self.Q_b) / ((self.cryogen.h_V - self.cryogen.h_L))
        v_z = 4 * BL_0 / (self.cryogen.rho_V * np.pi * self.d_i ** 2)
        return v_z

    @property
    def b_l_dot(self):
        """Returns evaporation rate in kg/s
        """
        return self.v_z * self.A_T * self.cryogen.rho_V

    @property
    def Q_b(self):
        if self.Q_b_fixed is None:
            "If Q_b_fixed is not set, calculate"
            return self.U_L * self.A_T * (self.T_air - self.cryogen.T_sat)
        else:
            return self.Q_b_fixed
