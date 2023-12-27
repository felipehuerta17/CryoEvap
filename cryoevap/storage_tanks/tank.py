import numpy as np
# from ..cryogens.cryogen import Cryogen
from ..cryogens.cryogen import Cryogen

# Open source Coolprop module
import CoolProp.CoolProp as CP

# Numerical integration of ODE systems
from scipy.integrate import solve_ivp

# Simpson's rule for integration with 2nd order accuracy
from scipy.integrate import simps

# Linear interpolant to reconstruct solutions
from scipy.interpolate import interp1d

# Plotting routines
from . import plots

class Tank:
    """ Class to be used as a container for the
    evaporation of pure cryogens"""

    def __init__(self, d_i, d_o, V, LF=0.97):
        """ Class constructor """
        # Compulsory parameters
        self.d_i = d_i  # [m] Tank internal diameter
        self.d_o = d_o  # [m] Tank external diameter
        self.V = V  # [m^3] Tank volume
        self.A_T = np.pi * d_i ** 2 / 4  # [m^2] cross section area
        self.l = V / self.A_T  # [m] Tank height
        self.LF = LF # Initial liquid filling
        self.cryogen = Cryogen()  # Empty Cryogen, see Cryogen class

        # Simulation control

        # Initialise dimensionless grid with 100 nodes as default
        self.z_grid = np.linspace(0, 1, 100)

        # Solution object
        self.sol = None

        # Time interval in seconds to plot vapour temperature profile
        self.time_interval = 3600

        # Store integrated quantities as dictionaries
        self.data = {'Time':[], 'Tv_avg': [], 'rho_V_avg': [],
                    'Q_VL': [],'Q_L': [], 'Q_V': [],
                    'V_L': [], 'B_L': [], 'BOG': [],
                    'drho_V_avg': [], 'dV_L': []}
        pass

    def set_HeatTransProps(self, U_L, U_V, T_air, Q_b_fixed=None, Q_roof=0, eta_w = 0.97):
        """Set separately tank heat transfer properties
        
        Inputs:
            U_L: liquid phase overall heat transfer coefficient / W m^-2 K ^-1
            U_V: vapour phase overall heat transfer coefficient / W m^-2 K ^-1
            T_air: Temperature of the surroundings / K
            Q_b_fixed: Fixed bottom heat ingress if specified 
        
        Returns:
            None
        """
        # Tank parameters
        self.U_L = U_L  
        self.U_V = U_V 
        self.Q_roof = Q_roof  
        self.T_air = T_air 

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

        # Integrate
        sol = solve_ivp(self.sys_isobaric, (0, t_f), IC, t_eval = t_eval, method='Radau', atol=1e-9, rtol=1e-6)

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
        Q_L_tot = self.Q_L_in + self.Q_b + self.Q_VL(self.cryogen.T_V) + self.Q_wi

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

        # Compute the wall heating considering the wall heat partitioning.
        # (1-eta_w) is the fraction of the external vapour heat ingress
        # that is transferred in the vapour 
        S_wall = (4*self.U_V*self.d_o/self.d_i**2) * (self.T_air - T[1:-1]) * (1-self.eta_w)

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
        Calculates evaporation rate after completing an evaporation
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
            return 1 /dH_LV * (self.Q_b + Q_L_in + self.data['Q_VL']) 
        
    def Q_VL(self, T_V):
        '''
        Calculate vapour to liquid heat transfer rate
        using the Fourier's law
        '''

        # Temperature gradient at the interface
        dz = (self.z_grid[1] - self.z_grid[0])*self.l_V
        dTdz_i = (-3 * T_V[0] + 4 * T_V[1] - T_V[2])/(2*dz)
        
        return self.cryogen.k_V_avg * self.A_T * dTdz_i
    
    # Plotting routines
    def plot_tv(self):
        '''
        Plots vapour temperature profile after running a solution
        '''

        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_tv() again')  

        # Produce vapour temperature plot
        plots.plot_tv(self)
        return
    
    def plot_V_L(self, unit='m3'):
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
        plots.plot_V_L(self, unit)
    
    def plot_BOG(self, unit='kg/h'):
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
        plots.plot_BOG(self, unit)
    
    def plot_Q(self, unit='kW'):
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
        plots.plot_Q(self, unit)

    
    def _reconstruct(self):
        '''
        Reconstructs integrated quantities such as the vapour
        to liquid heat transfer rate (Q_VL), the liquid heat ingress
        (Q_Lin) and the vapour heat ingress (Q_V)
        '''  
        Q_VL = []
        Tv_avg = []
        rho_V_avg = []

        # Extract time-steps in seconds
        self.data['Time'] = self.sol.t

        for i in range(0, len(self.sol.t)):
            # Get the temperature at this time step
            T_v = self.sol.y[1:, i]

            # Calculate and append Q_VL
            Q_VL.append(self.Q_VL(T_v))

            # Average vapour temperature
            Tv_avg.append(simps(T_v, self.z_grid))

            # Average vapour density
            self.cryogen.update_rho_V(self.z_grid, T_v)
            rho_V_avg.append(self.cryogen.rho_V_avg)
        
        # Extrapolate average vapour density for t = 0
        rho_V_avg[0] = self.interpolate(self.sol.t, rho_V_avg)
        rho_V_avg = np.array(rho_V_avg)

        # Vectorise
        self.data['V_L'] = self.sol.y[0]
        self.data['Tv_avg'] = np.array(Tv_avg)
        self.data['rho_V_avg'] = rho_V_avg
        self.data['Q_VL'] = np.array(Q_VL)

        # Reconstruct liquid and vapour heat ingresses
        l_L = self.sol.y[0] / self.A_T

        # Reconstruct: note that A_L, A_V are not used from the tank
        # but reconstructed from the liquid volume
        Q_L = self.U_L * (np.pi * self.d_o * l_L) * (self.T_air - self.cryogen.T_sat)

        # The driving force of Q_V is the average temperature
        Q_V = self.U_V * (np.pi * self.d_o * (self.l - l_L)) *( self.T_air - self.data['Tv_avg'])
        
        # Store reconstructed heat ingresses in the tank object
        self.data['Q_L'] = np.array(Q_L)
        self.data['Q_V'] = np.array(Q_V)
        self.data['Q_Vw'] = np.array(Q_V) * self.eta_w

        # Evaporation rate in kg/s
        self.data['B_L'] = self.evap_rate()
        
        # Average vapour density time derivative
        self.data['drho_V_avg'] = self.dydt(self.sol.t, rho_V_avg)

        # Liquid volume time derivative
        self.data['dV_L'] = self.dydt(self.sol.t, self.sol.y[0])

        # BOG rate calculation
        self.data['BOG'] = (self.data['B_L']
            + self.data['rho_V_avg'] * self.data['dV_L']
            - (self.V - self.data['V_L']) * self.data['drho_V_avg'])
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
    def Q_L_in(self):
        """ Liquid heat ingress through the walls
        in W """
        return self.U_L * self.A_L * (self.T_air - self.cryogen.T_sat)
    
    @property
    def Q_wi(self):
        """ Heat transferred directly to the vapour-liquid interface
        through the tank wall in contact to the vapour / W """
        return self.U_V * self.A_V * self.eta_w * (self.T_air - self.cryogen.Tv_avg)
    
    @property
    def v_z(self):
        """Update advective velocity with respect to tank liquid filling"""
        # Initial evaporation rate kg/s
        BL_0 = (self.Q_L_in + self.Q_b) / ((self.cryogen.h_V - self.cryogen.h_L))
        v_z = 4 * BL_0 / (self.cryogen.rho_V_sat * np.pi * self.d_i ** 2)
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
    
    @property
    def tau(self):
        '''Provides a conservative estimate of the 
        duration of the transient period
        of rapid vapour heating'''
        return self.l_V/self.v_z
