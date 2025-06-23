import numpy as np
# from ..cryogens.cryogen import Cryogen
from ..cryogens.cryogen import Cryogen

# Open source Coolprop module
import CoolProp.CoolProp as CP

# Numerical integration of ODE systems
from scipy.integrate import solve_ivp

# Simpson's rule for integration with 2nd order accuracy
from scipy.integrate import simpson

# Linear interpolant to reconstruct solutions
from scipy.interpolate import interp1d

# Plotting routines
from . import plots

# Datetime module for time management
import datetime

class Tank:
    """ Class to be used as a container for the
    evaporation of pure cryogens"""

    def __init__(self, d_i, d_o, V, LF=0.97):
        """
        Constructor for the Tank class.
        This initializes a tank object with its geometric properties, 
        cryogen state, simulation parameters, and data storage for 
        simulation results.
        
        Inputs:
        -------
            d_i : Internal diameter of the tank [m].
            d_o : External diameter of the tank [m].
            V   : Volume of the tank [m^3].
            LF  : Initial liquid filling fraction (default is 0.97).}

        --------
        """

        # Compulsory parameters
        self.d_i     = d_i                   # [m] Tank internal diameter
        self.d_o     = d_o                   # [m] Tank external diameter
        self.V       = V                     # [m^3] Tank volume
        self.A_T     = np.pi * d_i ** 2 / 4  # [m^2] cross section area
        self.l       = V / self.A_T          # [m] Tank height
        self.LF      = LF                    # Initial liquid filling
        self.cryogen = Cryogen()             # Empty Cryogen, see Cryogen class

        # Simulation control

        # Initialise dimensionless grid with 100 nodes as default
        self.z_grid = np.linspace(0, 1, 100)
        self.r_grid = np.linspace(0, 1, 3)

        # Solution object
        self.sol = None

        # Time interval in seconds to plot vapour temperature profile
        self.time_interval = 3600

        # Store integrated quantities as dictionaries
        self.data = {'Time' : [], 'Tv_avg' : [], 'rho_V_avg': [], 'Q_VL'      : [],
                     'Q_L'  : [], 'Q_V'    : [], 'Q_Vw'     : [], 'Q_env_w'   : [],
                     'Q_w_L': [], 'Q_tot'  : [], 'V_L'      : [], 'B_L'       : [],
                     'BOG'  : [], 'dV_L'   : [], 'Tw_avg'   : [], 'drho_V_avg': [],
                     'T_w_raw':[], 'T_V_raw': [],}
        pass

    def set_EnvironmentalProps(self, T_avg_day = None, T_range_day = None, h_env = 15, p_anual = None, start_date = None):
        """
        Set environmental temperature and properties.

        You can set a daily temperature cycle (T_avg_day, T_range_day), an annual temperature profile (p_anual), or both.
        - For a daily cycle, provide T_avg_day (K) and T_range_day (K).
        - For an annual profile, provide p_anual as a list of polynomial coefficients [a0, a1, a2] for temperature as a function of day of year.
        - To combine both, provide both daily and annual parameters.
        - Optionally, set h_env (W/m^2K) and start_date (YYYY-MM-DD) for annual phase alignment.

        Inputs:
        -------
            T_avg_day   : Average daily environmental temperature [K]
            T_range_day : Daily temperature range (T_max - T_min) [K]
            h_env       : Convective heat transfer coefficient [W/m^2K]
            p_anual     : Polynomial coefficients [a0, a1, a2] for annual temperature profile. 
                          The first coefficient a0 is the quadratic term, a1 is the linear term, and a2 is the constant term (array).
            start_date  : Start date for simulation (YYYY-MM-DD) for annual phase alignment (string)

        Returns:
        --------
            None

        """
        # Set environmental temperature properties for the day
        self.T_env_avg_day = T_avg_day          # [K]
        self.range_env_day = T_range_day        # [K]
        self.freq_env_day  = 2*np.pi/(24*3600)  # [Hz]
        self.poly_annual   = p_anual            # [a0, a1, a2] coefficients for the annual temperature polynomial
        self.h_env         = h_env              # [W/m^2K]
        self.start_date    = start_date         # [datetime] Start date of the simulation

        # Case 1: Just day temperature is set
        if p_anual is None:
            self.T_env  = lambda t: self.T_env_avg_day + 0.5*self.range_env_day*np.sin(self.freq_env_day*t)

        # Case 2: Just annual temperature is set
        elif T_range_day is None or T_avg_day is None:
            
            if self.start_date is None:
                self.T_env  = lambda t: 273.15 + (self.poly_annual[0]*((t/86400) % 365)**2 + self.poly_annual[1]*((t/86400) % 365) + self.poly_annual[2])

            else:
                reference_dt   = datetime.datetime(2024, 6, 20)
                start_date_dt  = datetime.datetime.strptime(self.start_date, '%Y-%m-%d')
                days_offset    = (start_date_dt - reference_dt).days
                self.T_env     = lambda t: 273.15 + (self.poly_annual[0]*(((t/86400) - days_offset) % 365)**2 + self.poly_annual[1]*(((t/86400) - days_offset) % 365) + self.poly_annual[2])

        # Case 3: Both day and annual temperatures are set
        else:
            if self.start_date is None:
                self.T_env  = lambda t: 273.15 +(self.poly_annual[0]*((t/86400) % 365)**2 + self.poly_annual[1]*((t/86400) % 365) + self.poly_annual[2]) + 0.5*self.range_env_day*np.sin(self.freq_env_day*t)  
            else:
                reference_dt   = datetime.datetime(2024, 6, 20)
                start_date_dt  = datetime.datetime.strptime(self.start_date, '%Y-%m-%d')
                days_offset    = (start_date_dt - reference_dt).days
                self.T_env     = lambda t: 273.15 + (self.poly_annual[0]*(((t/86400) - days_offset) % 365)**2 + self.poly_annual[1]*(((t/86400) - days_offset) % 365) + self.poly_annual[2]) + 0.5*self.range_env_day*np.sin(self.freq_env_day*t)
        
        # Fix the initial wall temperature to the environmental temperature
        if self.T_init:
            self.T_w = np.ones(len(self.r_grid)) * self.T_env(0)

        pass



    # @property
    # def Ra(self):
    #     # Descripcion
    #     return (g * self.cryogen.beta_L * delta_T * self.l**3) / (self.cryogen.mu/self.cryogen.rho_L * self.alpha )
    

    # def h_i_base(self, k_w):
    #     g = self.cryogen.g
    #     self.k_w = k_w
    #     self.alpha = self.k_w/(self.cryogen.rho_l*self.cryogen.cp_L)
    #     self.Pr = self.cryogen.Pr       
    #     f1 = (1 + (0.492*self.Pr)**(9/16))**(-16/9)
    #     # print(f'Raf1: {Ra*f1:.6e}')
    #     # Heat emission at lower surface 
    #     Nu_b = 0.6 * (self.Ra * f1)**(1/5)
    #     # print(f'Nusselt number: {Nu_b:.3e}')
    #     h_b = Nu_b * k_w / self.d_i 
    #     print(f'Internal heat transfer coefficient at the tank base, h_b: {h_b:.3e}', 'Wm^-2K^-1')
    #     return h_b
        

    
    def set_HeatTransProps(self, U_L, U_V, T_air, Q_b_fixed=None, Q_roof=0, eta_w = 0, k_w = 0.1, rho_w = 50, cp_w = 1000, h_L = 0, T_init = True):  
        """
        Set separately tank heat transfer properties
        
        Inputs:
        -------
            U_L       : liquid phase overall heat transfer coefficient / W m^-2 K ^-1
            U_V       : vapour phase overall heat transfer coefficient / W m^-2 K ^-1
            T_air     : Temperature of the surroundings / K
            Q_b_fixed : Fixed bottom heat ingress if specified 
            Q_roof    : Heat ingress through the roof / W
            eta_w     : Wall heat partitioning fraction
            k_w       : Wall thermal conductivity / W m^-1 K^-1
            rho_w     : Wall density / kg m^-3
            cp_w      : Wall specific heat capacity / J kg^-1 K^-1
            h_L       : Liquid phase heat transfer coefficient / W m^-2 K^-1
            T_init    : Set initial wall temperature as environmental temperature / True or False
        
        Returns:
        --------
            None

        """
        # Tank parameters
        self.U_L     = U_L  
        self.U_V     = U_V 
        self.Q_roof  = Q_roof  
        self.T_air   = T_air  
        self.H_L     = h_L

        # Wall properties
        self.k_w     = k_w
        self.rho_w   = rho_w
        self.cp_w    = cp_w
        self.alpha_w = k_w/(rho_w*cp_w)
        self.T_init  = T_init

        # Environmental property
        self.h_env = 0

        # The walls and roof materials are the same, hence, it is assumed
        # that their heat transfer coefficients are the same
        self.U_roof = U_V         

        # By default, the roof is thermally insulated
        self.roof_BC = "Neumann"  

        # In large scale applications, the tank bottom is heated by en
        # electric element at a constant rate to prevent ground freezing. 
        self.Q_b_fixed = Q_b_fixed 

        # Wall heat partitioning fraction
        self.eta_w = eta_w   

        # Set initial wall temperature
        if self.T_init:
            self.T_w = np.ones(len(self.r_grid)) * T_air
        else:
            self.T_w = np.ones(len(self.r_grid)) * self.cryogen.T_sat

        # Define T_env as constant by default
        self.T_env = lambda t: self.T_air

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
        Tv_0[-1] = ((2 * self.U_roof * (1-self.eta_w) * dz * self.T_env(0)/self.cryogen.k_V_avg +
                    4 * Tv_0[-2] - Tv_0[-3])/(3 + 2 * self.U_roof * (1-self.eta_w) * dz * self.cryogen.k_V_avg))

        # Initial wall temperature
        Tw_0 = np.ones(len(self.r_grid)) * self.T_w[0]

        dr = (self.r_grid[1] - self.r_grid[0]) * (self.d_o - self.d_i) * 0.5

        Tw_0[0]  = (self.H_L * Tv_0[0] + Tw_0[1] * (4 * self.k_w / (2 * dr)) - Tw_0[2] * (self.k_w / (2 * dr))) / (3 * self.k_w / (2 * dr) + self.H_L)
        Tw_0[-1] = (self.h_env * self.T_env(0) + Tw_0[-2] * (4 * self.k_w / (2 * dr)) - Tw_0[-3] * (self.k_w / (2 * dr))) / (3 * self.k_w / (2 * dr) + self.h_env)
        
        # Concatenate initial conditions in a single vector
        IC = np.concatenate([[VL_0], Tv_0, Tw_0])

        # Integrate
        sol = solve_ivp(self.sys_isobaric, (0, t_f), IC, t_eval = t_eval, method='RK45', atol=1e-6, rtol=1e-6)        

        # Set tank solution object with the volume and vapour temperature profiles
        # as a function of time
        self.sol = sol

        # Reconstruct integrated quantities
        self._reconstruct()
    
    def sys_liq_volume(self, t, y):
        '''
        ODE for the liquid volume. 
        '''
        # The liquid volume is the dependent variable
        V_L = y #y[0]

        # Updates liquid filling
        self.LF = V_L / self.V

        # Computes total heat ingress to the liquid
        Q_L_tot = self.Q_L_in(self.T_w) + self.Q_b(t) + self.Q_VL(self.cryogen.T_V) + self.Q_w_i(t)

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
        self.cryogen.Tv_avg = simpson(T, x = self.z_grid)

        # Update vapour temperature
        self.cryogen.T_V = T

        # Update average vapour density, thermal conductivity
        # and specific heat capacity using the Simpson's rule
        self.cryogen.update_rho_V(self.z_grid, T)
        self.cryogen.update_k_V(self.z_grid, T)
        self.cryogen.update_cp_V(self.z_grid, T)

        # Advective velocity
        v_z = self.v_z(t)

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

        # Compute the wall heating considering the wall heat partitioning.
        # (1-eta_w) is the fraction of the external vapour heat ingress
        # that is transferred in the vapour 
        S_wall = (4*self.U_V*self.d_o/self.d_i**2) * (self.T_env(t) - T[1:-1]) * (1-self.eta_w)

        # Update dT
        dT[1:-1] = alpha*d2T_dz2 - (v_z-v_int) * dT_dz + (alpha/self.cryogen.k_V_avg) * S_wall

        # DIFFERENTIAL BOUNDARY CONDITIONS
        # In the vapour-liquid interface the
        # temperature is constant for isobaric evaporation
        dT[0] = 0

        # 2nd order extrapolation
        # Assumes that the wall heat flow partitioning also applies at the tank roof
        if self.roof_BC == "Robin":
            dT[-1] = (4*dT[-2] - dT[-3])/(3 + 2 * self.U_roof * (1-self.eta_w) * dz)
        else:
            # Neumann boundary condition
            dT[-1] = (4*dT[-2] - dT[-3])/3
        
        return dT

    def sys_wall(self, t, y):
        '''
        Method of lines for the wall temperature
        '''
        # Liquid temperature
        T_L = y[0]

        # Wall temperature
        T = y[len(self.z_grid):]

        # Update wall temperature
        self.T_w = T

        # Dimension grid 
        r  = self.r_grid * (self.d_o - self.d_i) * 0.5 + self.d_i*0.5

        # Uniform spacing
        dr = (self.r_grid[1] - self.r_grid[0]) * (self.d_o - self.d_i) * 0.5
        
        # Number of grid points
        n = len(self.r_grid) 

        # Boundary corrections
        T[0]  = (self.H_L * T_L + T[1] * (4 * self.k_w / (2 * dr)) - T[2] * (self.k_w / (2 * dr))) / (3 * self.k_w / (2 * dr) + self.H_L)

        T[-1] = (self.h_env * self.T_env(t) + T[-2] * (4 * self.k_w / (2 * dr)) - T[-3] * (self.k_w / (2 * dr))) / (3 * self.k_w / (2 * dr) + self.h_env)

        # Initialise temperature change vector
        dT = np.zeros(n) 

        # Compute the differences
        dT_dr = (T[2:] - T[:-2]) / (2 * dr)

        # Compute the second derivatives
        d2T_dr2 = (T[2:] - 2*T[1:-1] + T[:-2]) / (dr**2)

        # Update dT
        dT[1:-1] = (self.alpha_w /r[1:-1]) * (dT_dr + r[1:-1] * d2T_dr2)

        # Boundary conditions
        dT[0]  = self.alpha_w * ( (self.H_L / (self.k_w * r[0])) * (T[0] - T_L) + (2*T[0] - 5*T[1] + 4*T[2] - T[3]) / dr**2 )
        dT[-1] = self.alpha_w * ( (self.h_env / (self.k_w * r[-1])) * (self.T_env(t) - T[-1]) +  (2*T[-1] - 5*T[-2] + 4*T[-3] - T[-4]) / dr**2)

        # # # Boundary conditions (extrapolation)
        # dT[-1] = 2*dT[-2] - dT[-3]
        # dT[0]  = 2*dT[1]  - dT[2]

        # Boundary conditions
        # dT_V_wall = 0  # Assuming no change in liquid temperature at the wall because is isobaric
        # dT[0]     = (1/(self.H_L + 3*self.k_w/(2*dr))) * (4*self.k_w/(2*dr) * dT[1] - self.k_w/(2*dr) * dT[2] + self.H_L * dT_V_wall)
        # dT[-1]    = (1/(self.h_env + 3*self.k_w/(2*dr))) * (4*self.k_w/(2*dr) * dT[-2] - self.k_w/(2*dr) * dT[-3] + self.h_env * self.dT_env(t))


        return dT

    def sys_isobaric(self, t, y):
        '''
        Constructs liquid volume + vapour temperature + wall temperature subsystem
        '''
        # Liquid volume derivative
        dV = self.sys_liq_volume(t, y[0])

        # ODE system with nodal vapour temperature derivatives
        dT_V =  self.sys_temperature(t, y[1:len(self.z_grid)+1])

        # Wall temperature derivative
        dT_w = self.sys_wall(t, y[1:])

        # Return right hand side of the ODE system
        return np.concatenate([[dV], dT_V, dT_w])

    def evap_rate(self, t):
        '''
        Calculates evaporation rate after completing an evaporation
        '''
        if self.sol is None:
            print("No solution available, perform tank.evaporate() to construct solution")
            return
        else:
            # Calculates latent heat of vaporisation
            dH_LV = self.cryogen.h_V - self.cryogen.h_L

            return 1 /dH_LV * (self.Q_b(t) + self.data['Q_L'] + self.data['Q_VL'] + self.data['Q_Vw']) 
        
    def Q_VL(self, T_V):
        '''
        Calculate vapour to liquid heat transfer rate
        using the Fourier's law
        '''
        # Temperature gradient at the interface
        dz = (self.z_grid[1] - self.z_grid[0])*self.l_V
        dTdz_i = (-3 * T_V[0] + 4 * T_V[1] - T_V[2])/(2*dz)
        
        return self.cryogen.k_V_avg * self.A_T * dTdz_i
    
    def Q_L_in(self, T_w):
        """ Liquid heat ingress through the walls in W """

        # Temperature gradient at the internal wall
        dr = (self.r_grid[1] - self.r_grid[0])*(self.d_o - self.d_i)*0.5
        dTdr_i = (3 * T_w[0] - 4 * T_w[1] + T_w[2])/(2*dr)

        return - self.k_w * self.A_L * dTdr_i

    def Q_w_i(self, t):
        """ Heat transferred directly to the vapour-liquid interface
        through the tank wall in contact to the vapour / W """
        return self.U_V * self.A_V * self.eta_w * (self.T_env(t) - self.cryogen.Tv_avg)
    
    def Q_env_w(self, T_w):
        """ Heat transferred from the environment to the external wall of the tank / W """
        # Temperature gradient at the internal wall
        dr     = (self.r_grid[1] - self.r_grid[0])*(self.d_o - self.d_i)*0.5
        dTdr_o = ( 3 * T_w[-1] - 4 * T_w[-2] + T_w[-3]) / (2*dr) 
        return - dTdr_o * self.k_w * self.A_L
    
    def Q_b(self, t):
        """Calculate bottom heat ingress"""
        if self.Q_b_fixed is None:
            "If Q_b_fixed is not set, calculate"
            return self.U_L * self.A_T * (self.T_env(t) - self.cryogen.T_sat)
        return self.Q_b_fixed
    
    def v_z(self, t):
        """Calculate advective velocity with respect to tank liquid filling"""
        # Initial evaporation rate kg/s
        BL_0 = (self.Q_L_in(self.T_w) + self.Q_b(t) + self.Q_w_i(t)) / ((self.cryogen.h_V - self.cryogen.h_L))
        v_z = 4 * BL_0 / (self.cryogen.rho_V_sat * np.pi * self.d_i ** 2)
        return v_z   
    
    # Plotting routines
    def plot_tv(self, t_unit='s', hour = 0):
        '''
        Plots vapour temperature profile after running a solution
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_tv() again')  

        # Produce vapour temperature plot
        plots.plot_tv(self, t_unit, hour)
        return

    def plot_tw(self, t_unit='s', hour = 0):
        '''
        Plots wall temperature profile after running a solution
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_tv() again')  

        # Produce vapour temperature plot
        plots.plot_tw(self, t_unit, hour)
        return
    
    def plot_V_L(self, unit='m3', t_unit='s'):
        '''
        Plots liquid volume as a function of time

        Inputs:
            unit: 'm3', 'L', 'mL'
            
        Returns:
            None
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_V_L() again')  

        # Produce liquid volume plot
        plots.plot_V_L(self, unit, t_unit)
    
    def plot_BOG(self, unit='kg/h', t_unit='s'):
        '''
        Plots boil off gas rate as a function of time

        Inputs:
            unit: 'kg/h', 'kg/s', 'g/s'
        Returns:
            None
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_BOG() again')  

        # Produce liquid volume plot
        plots.plot_BOG(self, unit, t_unit)
    
    def plot_Q(self, unit='kW', t_unit='s'):
        '''
        Plots vapour to liquid heat transfer rate

        Inputs:
            tank: Tank object with a sol object 
            unit: [Q_V, Q_L, Q_VL], units. kW or W
        
        Returns:
            None:
        '''
        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_QVL() again') 
        
        # Produce liquid volume plot
        plots.plot_Q(self, unit, t_unit)
    
    def plot_tv_BOG(self, t_unit='s'):
        '''
        Plots average vapour temperature and boil-off gas temperature
        as a function of time

        Inputs:
            t_unit(string): Time units among 's', 'min', 'h', 'd', 'w'.
        Returns:
            None
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_BOG() again')  

        # Produce liquid volume plot
        plots.plot_tv_BOG(self, t_unit)

    def plot_Q_w(self, unit='kW', t_unit = 's'):
        '''
        Plots wall heat ingress as a function of time

        Inputs:
            unit: 'kW', 'W'
        Returns:
            None
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_Q_w() again')  

        # Produce liquid volume plot
        plots.plot_Q_w(self, unit, t_unit)

    def plot_T_w_avg(self, t_unit = 's'):
        '''
        Plots average wall temperature as a function of time

        Inputs:
            t_unit: 's', 'h', 'd', 'w'
        Returns:
            None
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_Q_w() again')  

        # Produce liquid volume plot
        plots.plot_T_w_avg(self, t_unit)

    # Reconstruction of results
    def _reconstruct(self):
        '''
        Reconstructs integrated quantities such as the vapour
        to liquid heat transfer rate (Q_VL), the liquid heat ingress
        (Q_Lin), bottom heat ingress and the vapour heat ingress (Q_V).

        It also reconstructs local variables of interest, such as T_BOG.
        '''  
        Tv_avg    = []
        Tw_avg    = []
        rho_V_avg = []
        T_BOG     = []
        Q_VL      = []
        Q_L_in    = []
        Q_env_w   = []

        # Extract time-steps in seconds
        self.data['Time'] = self.sol.t

        # Reconstruct liquid length for heat transfer calculations
        l_L = self.sol.y[0] / self.A_T

        for i in range(len(self.sol.t)):
            # Get the temperature at this time step
            T_v = self.sol.y[1:len(self.z_grid) + 1, i]
            T_w = self.sol.y[len(self.z_grid) + 1:, i]
            # Calculate and append Q_VL

            # Update vapour thermal conductivity
            self.cryogen.update_k_V(self.z_grid, T_v)
            
            # Calculate vapour temperature gradient from the vapour length
            # at the desired timestep
            dz     = (self.z_grid[1] - self.z_grid[0])* ((self.l - l_L[i]))
            dTdz_i = (3 * T_v[0] - 4 * T_v[1] + T_v[2])/(2*dz)    
            
            # Append Q_VL calculated using the Fourier's law
            Q_VL.append(-self.cryogen.k_V_avg * self.A_T * dTdz_i)

            # Average vapour temperature
            Tv_avg.append(simpson(T_v, x = self.z_grid))

            # Average vapour density
            self.cryogen.update_rho_V(self.z_grid, T_v)
            rho_V_avg.append(self.cryogen.rho_V_avg)

            # BOG temperature
            T_BOG.append(T_v[-1])

            # Average wall temperature
            Tw_avg.append(simpson(T_w, x = self.r_grid))

            # Calculate wall temperature gradient 
            dr     = (self.r_grid[1] - self.r_grid[0])*(self.d_o - self.d_i)*0.5
            dTdr_i = ( 3 * T_w[0] - 4 * T_w[1] + T_w[2])/(2*dr)
            dTdr_o = ( 3 * T_w[-1] - 4 * T_w[-2] + T_w[-3])/(2*dr) 

            # Calculate wall areas
            A_int = np.pi * self.d_i * l_L[i]
            A_ext = np.pi * self.d_o * l_L[i]

            # Append Q_LW and Q_Wenv calculated using the Fourier's law
            Q_L_in.append(-self.k_w * dTdr_i * A_int) # REVISAR Heat transferred to the liquid from the internal wall
            Q_env_w.append(-self.k_w * dTdr_o * A_ext) # REVISAR Heat transferred to the wall from the environment

        
        # Extrapolate average vapour density for t = 0
        rho_V_avg[0] = self.interpolate(self.sol.t, rho_V_avg)

        # Vectorise
        self.data['V_L']       = self.sol.y[0]
        self.data['Tv_avg']    = np.array(Tv_avg)
        self.data['rho_V_avg'] = np.array(rho_V_avg)
        self.data['T_BOG']     = np.array(T_BOG)
        self.data['Tw_avg']    = np.array(Tw_avg)

        # Reconstruct liquid and vapour heat ingresses.
        # Note that A_L, A_V are not used from the tank
        # but reconstructed from the liquid volume

        # The driving force of Q_V is the average temperature
        Q_V = self.U_V * (np.pi * self.d_o * (self.l - l_L)) *( self.T_env(self.sol.t) - self.data['Tv_avg']) 
        
        # Store reconstructed heat ingresses in the tank object
        self.data['Q_VL']     = np.array(Q_VL)
        self.data['Q_L']      = np.array(Q_L_in)
        self.data['Q_V']      = np.array(Q_V)
        self.data['Q_Vw']     = np.array(Q_V) * self.eta_w
        self.data['Q_env_w']  = np.array(Q_env_w)
        self.data['Q_w_L']    = np.array(Q_L_in) * -1
        self.data['Q_tot']    = self.data['Q_L'] + self.data['Q_Vw'] + self.data['Q_VL'] +  self.Q_b(self.sol.t)
        
        # Evaporation rate in kg/s
        self.data['B_L'] = np.array(self.evap_rate(self.sol.t))
        
        # Average vapour density time derivative
        self.data['drho_V_avg'] = self.dydt(self.sol.t, np.array(rho_V_avg))

        # Liquid volume time derivative
        self.data['dV_L'] = self.dydt(self.sol.t, self.sol.y[0])

        # BOG rate calculation
        self.data['BOG'] = (self.data['B_L']
            + self.data['rho_V_avg'] * self.data['dV_L']
            - (self.V - self.data['V_L']) * self.data['drho_V_avg'])

        # Wall temperature raw
        self.data['T_w_raw'] = self.sol.y[len(self.z_grid) + 1:, :]
        self.data['t_raw']   = self.sol.t
        self.data['T_V_raw'] = self.sol.y[1:len(self.z_grid) + 1, :]

        return
    
    def dydt(self, t, y):
        '''
        Calculate state variable derivative finite
        difference approximation dydt. from
        the numerical solution y evaluated in the
        discretized time domain t 
        
        Inputs:
            y: Dependent variable 
            t: Discretized time domain
        
        Returns:
            dy: State variable derivative
        '''
        dy = np.zeros(len(t))
        # dydt|t=0 = 0
        # Assume equally spaced dt
        dt = t[1] - t[0]
        # Second order central differences
        dy[1:-1] = (y[2:] - y[0:-2])/(2*dt)
        # Second order backward differences
        dy[-1] = (3*y[-1] - 4*y[-2] + y[-3])/(2*dt)
        return dy

    def interpolate(self, t, y, t_int = 0):
        '''
        Interpolate/extrapolate dependent variable 
        in the index 0 for purposes of
        reconstructing solutions

        Inputs:
            t_int: interpolation time / s
        '''
        # Create a linear interpolant
        linear_interp = interp1d(t[1:], y[1:],
                                  kind='linear', fill_value='extrapolate')
        return linear_interp(t_int)


    # Properties
    @property
    def l_V(self):
        """Update liquid filling and vapour length"""
        return self.l * (1 - self.LF)  # [m] sets vapour length

    @property
    def A_L(self):
        """Tank wall area in contact with the liquid"""
        return np.pi * self.d_o * self.l * self.LF

    @property
    def A_V(self):
        """Tank wall area in contact with the vapour"""
        return np.pi * self.d_o * self.l * (1-self.LF)

    @property
    def Q_L_in_(self):
        """ Liquid heat ingress through the walls
        in W, considering the environmental temperature as constant """
        return self.U_L * self.A_L * (self.T_air - self.cryogen.T_sat)
    
    @property
    def Q_wi_(self):
        """ Heat transferred directly to the vapour-liquid interface
        through the tank wall in contact to the vapour / W, considering the environmental temperature as constant """
        return self.U_V * self.A_V * self.eta_w * (self.T_air - self.cryogen.Tv_avg)
    
    @property
    def v_z_(self):
        """Update advective velocity with respect to tank liquid filling,
        considering the environmental temperature as constant"""
        # Initial evaporation rate kg/s
        BL_0 = (self.Q_L_in_ + self.Q_b_ + self.Q_wi_) / ((self.cryogen.h_V - self.cryogen.h_L))
        v_z = 4 * BL_0 / (self.cryogen.rho_V_sat * np.pi * self.d_i ** 2)
        return v_z

    @property
    def b_l_dot(self):
        """Returns evaporation rate in kg/s,
        considering the environmental temperature as constant"""
        return self.v_z_ * self.A_T * self.cryogen.rho_V_avg

    @property
    def Q_b_(self):
        """Calculate bottom heat ingress, considering the environmental temperature as constant"""
        if self.Q_b_fixed is None:
            "If Q_b_fixed is not set, calculate"
            return self.U_L * self.A_T * (self.T_air - self.cryogen.T_sat)
        else:
            return self.Q_b_fixed
    
    @property
    def tau(self):
        """Provides a conservative estimate of the 
        duration of the transient period
        of rapid vapour heating, considering the environmental temperature as constant"""
        return 2 * self.l_V/self.v_z_
