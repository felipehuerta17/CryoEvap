import numpy as np
# from ..cryogens.cryogen import Cryogen
from ..cryogens.cryogen import Cryogen

# Open source Coolprop module
import CoolProp.CoolProp as CP

# Numerical integration of ODE systems
from scipy.integrate import solve_ivp

# Simpson's rule for integration with 2nd order accuracy
from scipy.integrate import simps

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

        # Initialise dimensionless grid with 100 nodes as default
        self.z_grid = np.linspace(0, 1, 100)

        # Solution object
        self.sol = None

        # Time interval in second to plot vapour temperature profile
        self.time_interval = 3600

        pass

    def set_HeatTransProps(self, U_L, U_V, T_air, Q_b_fixed=None, Q_roof=0):
        """Set separately tank heat transfer properties
        Usage: set_HeatTransProps(self, U_L, U_V, Q_b, Q_roof, T_air)"""
        self.U_L = U_L  # [W*m^-2*K^-1]Overall heat transfer coefficient
        # for the liquid phase stored in the tank
        self.U_V = U_V  # [W*m^-2*K^-1] Overall heat transfer coefficient
        
        # Roof heat transfer coefficient
        self.U_roof = U_V
        # for the vapour phase stored in the tank
        # Is there any Q_b_fixed?
        self.Q_b_fixed = Q_b_fixed
        self.Q_roof = Q_roof  # [W] Heat ingress through the roof
        self.T_air = T_air  # [K] Temperature of the surrounding air /K
        pass

    def evaporate(self, t_f):
        ''' Simulates the isobaric evaporation of the stored cryogen for the time t_f
        '''
        
        # Convert initial vapour temperature in an array
        self.cryogen.T_V = self.cryogen.T_V * np.ones(len(self.z_grid))

        # Define time span to evaluate the solution
        t_eval = np.arange(0, t_f + 1, self.time_interval)

        # Initial liquid volume
        VL_0 = self.V * self.LF

        # Initial vapour temperature
        Tv_0 = np.ones(len(self.z_grid)) * self.cryogen.T_sat

        # Interfacial temperature
        Tv_0[0] = self.cryogen.T_sat

        # Roof temperature
        dz = self.z_grid[1] - self.z_grid[0]

        # Robin BC initial condition
        Tv_0[-1] = ((2 * self.U_roof * dz * self.T_air/self.cryogen.k_V_avg +
                    4 * Tv_0[-2] - Tv_0[-3])/(3 + 2 * self.U_roof * dz * self.cryogen.k_V_avg))

        # Concatenate initial conditions in a single vector
        IC = np.append(VL_0, Tv_0)

        sol = solve_ivp(self.sys_isobaric, (0, t_f), IC, t_eval = t_eval, method='Radau', atol=1e-6, rtol=1e-6)
        self.sol = sol
    
    def sys_liq_volume(self, t, y):
        '''
        ODE for the liquid volume. 
        '''
        # The liquid volume is the dependent variable
        V_L = y #y[0]

        # Updates liquid filling
        self.LF = V_L / self.V

        # Computes total heat ingress to the liquid
        Q_L_tot = self.Q_L_in + self.Q_b + self.Q_VL(self.cryogen.T_V)

        # Calculates latent heat of vaporisation
        dH_LV = self.cryogen.h_V - self.cryogen.h_L

        # Returns RHS of ODE that governs the liquid volume
        return -1 / self.cryogen.rho_L * (Q_L_tot/dH_LV)
    
    def sys_temperature(self, t, y):
        '''
        Method of lines for the vapour temperature
        '''
        # Define vapour temperature as de dependent variable
        T = y

        # Print time if self.print_time is True
        # print("t = %.3f s" % t)

        # Grid and properties initialization
        L_dry = self.l*(1-self.LF) # Dry height of the tank
        
        # Update average vapour temperature using Simpson's rule
        self.cryogen.Tv_avg = simps(T, self.z_grid)

        # Update vapour temperature
        self.cryogen.T_V = T

        # Update average vapour density, thermal conductivity
        # and specific heat capacity using the Simpson's rule

        self.cryogen.update_rho_V(self.z_grid, T)
        self.cryogen.update_k_V(self.z_grid, T)
        self.cryogen.update_cp_V(self.z_grid, T)

        # Advective velocity
        v_z = self.v_z

        # Interface velocity
        v_int = v_z * (self.cryogen.rho_V_avg/self.cryogen.rho_L)

        # Vapour thermal diffusivity
        alpha = self.cryogen.k_V_avg/(self.cryogen.rho_V_avg*self.cryogen.cp_V_avg) 

        # Uniform spacing
        dz = (self.z_grid[1] - self.z_grid[0])*L_dry

        # Number of grid points
        n = len(self.z_grid) 

        # Initialise temperature change vector
        dT = np.zeros(n) 

        # Compute the differences
        dT_dz = (T[1:-1] - T[:-2]) / dz

        # Compute the second derivatives
        d2T_dz2 = (T[:-2] - 2*T[1:-1] + T[2:]) / dz**2

        # Compute the wall heating
        S_wall = (4*self.U_V*self.d_o/self.d_i**2) * (self.T_air - T[1:-1])

        # Update dT
        dT[1:-1] = alpha*d2T_dz2 - (v_z-v_int) * dT_dz + (alpha/self.cryogen.k_V_avg) * S_wall

        # DIFFERENTIAL BOUNDARY CONDITIONS
        # In the vapour-liquid interface the
        # temperature is constant for isobaric evaporation
        dT[0] = 0

        # 2nd order extrapolation
        if self.roof_BC == "Robin":
            dT[-1] = (4*dT[-2] - dT[-3])/(3 + 2*self.U_roof * dz)
        else:
            # Neumann boundary condition
            dT[-1] = (4*dT[-2] - dT[-3])/3
            #dT[-1] = dT[-2]
        
        return dT
    
    def sys_isobaric(self, t, y):
        '''
        Constructs liquid volume + vapour temperature subsystem
        '''
        # Liquid volume derivative
        # dV = self.sys_liq_volume(self, t, y[0])
        dV = self.sys_liq_volume(t, y[0])
        # ODE system with nodal vapour temperature derivatives
        #dT_V =  self.sys_temperature(self, t, y[1:])
        dT_V =  self.sys_temperature(t, y[1:])

        # Return right hand side of the ODE system
        return np.append(dV, dT_V)

    def evap_rate(self):
        '''
        Calculates evaporation rate
        '''
        if self.sol is None:
            print("No solution available, perform tank.evaporate() to construct solution")
            return
        else:
            # Calculates latent heat of vaporisation
            dH_LV = self.cryogen.h_V - self.cryogen.h_L

            # Extracts liquid volume
            V_L = self.sol.y[0]

            Q_L_in = 4 * V_L * self.d_o/(self.d_i**2) * self.U_L * (self.T_air - self.cryogen.T_sat)
            return 1 /dH_LV * (self.Q_b + Q_L_in) 
    
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
    
    def Q_VL(self, T_V):
        '''
        Calculate vapour to liquid heat transfer rate
        using the Fourier's law
        '''

        # Temperature gradient at the interface
        dz = (self.z_grid[1] - self.z_grid[0])*self.l_V
        dTdz_i = (-3 * T_V[0] + 4 * T_V[1] - T_V[2])/(2*dz)
        
        return self.cryogen.k_V_avg * self.A_T * dTdz_i


    @property
    def v_z(self):
        """Update advective velocity with respect to tank liquid filling"""
        # Initial evaporation rate kg/s
        BL_0 = (self.Q_L_in + self.Q_b) / ((self.cryogen.h_V - self.cryogen.h_L))
        v_z = 4 * BL_0 / (self.cryogen.rho_V_avg * np.pi * self.d_i ** 2)
        return v_z

    @property
    def b_l_dot(self):
        """Returns evaporation rate in kg/s
        """
        return self.v_z * self.A_T * self.cryogen.rho_V_avg

    @property
    def Q_b(self):
        if self.Q_b_fixed is None:
            "If Q_b_fixed is not set, calculate"
            return self.U_L * self.A_T * (self.T_air - self.cryogen.T_sat)
        else:
            return self.Q_b_fixed
