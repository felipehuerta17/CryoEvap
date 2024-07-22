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

from scipy.optimize import fsolve

# Plotting routines
from . import plots

class Tank:
    """ Class to be used as a container for the
    evaporation of pure cryogens"""

    def __init__(self, d_i, d_o, V, vapour_geometry,liquid_geometry,LF=0.97):
        """ Class constructor """
        # Compulsory parameters
        self.d_i = d_i  # [m] Tank internal diameter
        self.d_o = d_o  # [m] Tank external diameter
        self.V = V  # [m^3] Tank volume
        self.LF = LF # Initial liquid filling
        self.Geo_v = vapour_geometry
        self.Geo_l = liquid_geometry

        if self.Geo_l=="cylindrical":
            self.A_T = np.pi * d_i ** 2 / 4  # [m^2] cross section area
            self.l = V / self.A_T  # [m] Tank height
        elif self.Geo_l=="spherical":
            self.l = d_i
            self.z = np.roots([-np.pi/3,np.pi*self.l/2,0,-self.V*self.LF])[1]
            assert(self.z<self.l and self.z>=0 and isinstance(self.z,float))
            self.A_T = abs(np.pi*(2*self.d_i/2 * self.z - self.z**2))
        elif self.Geo_l=="horizontal cylinder":
            self.l = V/(np.pi*(d_i/2)**2)
            self.z = fsolve(self.h_cylinder_vol,[self.LF*self.d_i])[0]
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

    def set_HeatTransProps(self, U_L, U_V, T_air, Q_b_fixed=None, Q_roof=0, eta_w = 0):
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
        Tv_0[-1] = ((2 * self.U_roof * (1-self.eta_w) * dz * self.T_air/self.cryogen.k_V_avg +
                    4 * Tv_0[-2] - Tv_0[-3])/(3 + 2 * self.U_roof * (1-self.eta_w)
                                              * dz * self.cryogen.k_V_avg))

        # Concatenate initial conditions in a single vector
        IC = np.append(VL_0, Tv_0)

        self.vz0 = self.v_z

        # Integrate
        sol = solve_ivp(self.sys_isobaric, (0, t_f), IC, t_eval = t_eval, method='RK45', atol=1e-6, rtol=1e-6) #change to RK45

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

        if self.Geo_l == "spherical":
            self.z = np.roots([-np.pi/3,np.pi*self.l/2,0,-self.V*self.LF])[1]
            assert(self.z<self.l and self.z>=0 and isinstance(self.z,float))
            self.A_T = abs(np.pi*(2*self.l/2 * self.z - self.z**2))
        elif self.Geo_l == "horizontal cylinder":
            self.z = fsolve(self.h_cylinder_vol,[self.z])[0]
            self.A_T = self.l*np.sqrt(abs(2*self.z*self.d_i/2 - self.z**2))

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
        if self.Geo_l == "cylindrical":
            L_dry = self.l*(1-self.LF) # Dry height of the tank
        elif self.Geo_l == "spherical" or self.Geo_l == "horizontal cylinder":
            L_dry = self.l - self.z

        # Update average vapour temperature using Simpson's rule
        if self.Geo_v == "spherical" or self.Geo_v == "horizontal cylinder":
            h_grid = self.z_grid*(L_dry)+self.z - (self.z_grid[1]-self.z_grid[0])*L_dry/2
            #h_grid = self.z_grid*(L_dry)+self.z
            radius = np.sqrt(abs(2*h_grid*self.l/2 - h_grid**2))
            self.cryogen.Tv_avg = simps(T*radius, self.z_grid)/simps(radius,self.z_grid)
        elif self.Geo_v == "cylindrical":
            radius = None
            self.cryogen.Tv_avg = simps(T,self.z_grid)

        # Update vapour temperature
        self.cryogen.T_V = T
        #print(self.cryogen.rho_V_avg)
        # Update average vapour density, thermal conductivity
        # and specific heat capacity using the Simpson's rule
        self.cryogen.update_rho_V(self.z_grid, T, radius)
        self.cryogen.update_k_V(self.z_grid, T, radius)
        self.cryogen.update_cp_V(self.z_grid, T, radius)
        #print(self.cryogen.rho_V_avg)
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

        if self.Geo_v == "cylindrical":
            #S_wall for a cylinder
            S_wall = (4*self.U_V*self.d_o/self.d_i**2) * (self.T_air - T[1:-1]) * (1-self.eta_w)

            # Update dT
            dT[1:-1] = alpha*d2T_dz2 - (v_z-v_int) * dT_dz + (alpha/self.cryogen.k_V_avg) * S_wall


        elif self.Geo_v == "spherical":
            z = self.z_grid[1:-1]*L_dry + self.z
            #S_wall for a sphere
            S_wall = 2*self.U_V*(self.T_air - T[1:-1]) * (1-self.eta_w) / (self.cryogen.rho_V_avg * self.cryogen.cp_V_avg *
                                                                           np.sqrt(abs(2*(self.l/2)*(z) - (z)**2)))
            #shape factor
            #shape = 2*(abs(self.l/2 - self.z_grid[1:-1]*L_dry)) * (-self.cryogen.rho_V_avg*self.cryogen.cp_V_avg*v_z*T[1:-1] 
            #            + self.cryogen.k_V_avg*dT_dz)/(self.cryogen.rho_V_avg*self.cryogen.cp_V_avg*abs((self.z_grid[1:-1]*L_dry)**2 - 2*(self.z_grid[1:-1]*L_dry)*self.l/2))
            
            # shape factor conduction only
            # shape = 2*(abs(self.l/2 - self.z_grid[1:-1]*L_dry)) * (
            #             + self.cryogen.k_V_avg*dT_dz)/(self.cryogen.rho_V_avg*self.cryogen.cp_V_avg*abs((self.z_grid[1:-1]*L_dry)**2 - 2*(self.z_grid[1:-1]*L_dry)*self.l/2))

            v_z = (self.v_z/self.z)*(z)*np.exp(-(z)*self.d_i/2)
            
            r = np.sqrt(abs((z)*2*self.l/2 - (z)**2))
            dr_dz = (self.d_i/2 - z)/np.sqrt(abs(2*(z)*self.d_i/2 - (z)**2))
            
            #Update dT
            #dT[1:-1] = alpha*d2T_dz2 - (v_z-v_int) * dT_dz + S_wall + shape
            dT[1:-1] = alpha*d2T_dz2 - v_z*dT_dz/(r**2) + 2*dr_dz*alpha*dT_dz/r + S_wall
            #if any(dT<0):
            #    print(dT," ",alpha*d2T_dz2," ",(v_z-v_int)*dT_dz," ",S_wall," ",shape," ",t)


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

            return 1 /dH_LV * (self.Q_b + self.data['Q_L'] + self.data['Q_VL'] + self.data['Q_Vw'])
        
    def h_cylinder_vol(self,z):
            return [self.l*(((z-self.d_i/2)*np.sqrt(-z*(z-self.d_i)))/2 + np.arcsin((z-self.d_i/2)/(self.d_i/2))*(self.d_i/2)**2 +
                       ((self.d_i/2)**2)*np.pi/2) - self.LF*self.V]
        
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

    def plot_l_L(self):
        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_lL() again')
        plots.plot_l_L(self)

    def plot_A_T(self):
        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_AT() again')
        plots.plot_A_T(self)
    
    def plot_LF(self):
        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_LF() again')
        plots.plot_LF(self)
    
    def plot_rho_V_avg(self):
        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_rhoVavg() again')
        plots.plot_rho_V_avg(self)

    def plot_vz(self):
        if self.sol is None:
            raise TypeError('The solution object tank.sol does not exist.\n'
                            'Run tank.evaporate(t) to generate a solution\n'
                            'and thereafter run tank.plot_rhoVavg() again')
        plots.plot_vz(self)

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

        # Reconstruct liquid length for heat transfer calculations
        if self.Geo_l=="cylindrical":
            l_L = self.sol.y[0] / self.A_T
        elif self.Geo_l =="spherical":
            l_L = np.ones_like(self.sol.y[0])
            for i in range(len(self.sol.y[0])):
                l_L[i] = np.roots([-np.pi/3,np.pi*self.l/2,0,-self.sol.y[0][i]])[1]
                assert(l_L[i]<self.l and l_L[i]>=0)
        elif self.Geo_l == "horizontal cylinder":
            l_L = np.ones_like(self.sol.y[0])
            for i in range(len(self.sol.y[0])):
                self.LF = self.sol.y[0][i]/self.V
                l_L[i] = fsolve(self.h_cylinder_vol,[self.LF*self.d_i])[0]

        vz_avg = []

        for i in range(0, len(self.sol.t)):
            # Get the temperature at this time step
            T_v = self.sol.y[1:, i]

            # Calculate and append Q_VL
            
            # Update vapour thermal conductivity
            if self.Geo_v == "spherical" or self.Geo_v == "horizontal cylinder":
                #h_grid = self.z_grid*(self.l-l_L[i])+l_L[i]
                h_grid = self.z_grid*(self.l-l_L[i])+l_L[i] - (self.z_grid[1]-self.z_grid[0])*(self.l-l_L[i])/2
                radius = np.sqrt(abs(2*h_grid*self.l/2 - h_grid**2))
            else:
                radius = None
            self.cryogen.update_k_V(self.z_grid, T_v, radius)


            # Calculate vapour temperature gradient from the vapour length
            # at the desired tiemstep
            dz = (self.z_grid[1] - self.z_grid[0])* ((self.l - l_L[i]))
            dTdz_i = (-3 * T_v[0] + 4 * T_v[1] - T_v[2])/(2*dz)

            # Append Q_VL calculated using the Fourier's law
            Q_VL.append(self.cryogen.k_V_avg * self.A_T * dTdz_i)

            # Average vapour temperature
            if radius is None:
                Tv_avg.append(simps(T_v, self.z_grid))
            else:
                Tv_avg.append(simps(T_v*radius, self.z_grid)/simps(radius,self.z_grid))
                zed = self.z_grid*(self.d_i-l_L[i]) + l_L[i]
                vz = (self.vz0/l_L[i])*zed*np.exp(-(zed-l_L[i])*self.d_i/2)
                vz_avg.append(simps(vz*radius, self.z_grid)/simps(radius,self.z_grid))
            

            self.cryogen.update_rho_V(self.z_grid, T_v, radius)

            # Average vapour density
            rho_V_avg.append(self.cryogen.rho_V_avg)

        # Extrapolate average vapour density for t = 0
        rho_V_avg[0] = self.interpolate(self.sol.t, rho_V_avg)
        rho_V_avg = np.array(rho_V_avg)

        # Vectorise
        self.data["z"] = l_L
        if self.Geo_l=='spherical':
            self.data["A_T"] = np.pi*abs(2*l_L*self.l/2 - l_L**2) #sphere only
            self.data['vz_avg'] = np.array(vz_avg)
        self.data['V_L'] = self.sol.y[0]
        self.data["LF"] = self.sol.y[0]/self.V
        self.data['Tv_avg'] = np.array(Tv_avg)
        self.data['dTV_avg'] = self.dydt(self.sol.t,np.array(Tv_avg))
        self.data['rho_V_avg'] = rho_V_avg
        self.data['Q_VL'] = np.array(Q_VL)

        # Reconstruct liquid and vapour heat ingresses.
        # Note that A_L, A_V are not used from the tank
        # but reconstructed from the liquid volume

        # The driving force of Q_V is the average temperature
        if self.Geo_v == "cylindrical" and self.Geo_l == "cylindrical":
            Q_L = self.U_L * (np.pi * self.d_o * l_L) * (self.T_air - self.cryogen.T_sat)
            Q_V = self.U_V * (np.pi * self.d_o * (self.l - l_L)) *( self.T_air - self.data['Tv_avg'])
        elif self.Geo_v == "spherical" and self.Geo_l == "spherical":
            Q_L = self.U_L * (np.pi * self.d_o * (l_L + (self.d_o-self.d_i)/2)) * (self.T_air - self.cryogen.T_sat)
            Q_V = self.U_V * (np.pi * self.d_o * (self.l - l_L + (self.d_o - self.d_i)/2)) *( self.T_air - self.data['Tv_avg'])
        elif self.Geo_v == "spherical" and self.Geo_l == "cylindrical":
            Q_L = self.U_L * (np.pi * self.d_o * l_L) * (self.T_air - self.cryogen.T_sat)
            Q_V = self.U_V*((np.pi*self.d_o**2)/2 + ((self.d_i/2 - l_L)*np.pi*self.d_o))*(self.T_air - self.data['Tv_avg'])
        elif self.Geo_l == "horizontal cylinder" and self.Geo_v == "horizontal cylinder":
            SA_L = ((((self.d_o/2)**2)*((2*np.arccos(1-self.z/(self.d_o/2))))-np.sin(2*np.arccos(1-self.z/(self.d_o/2))))+
                              (self.d_o/2)*np.arccos(1-self.z/(self.d_o/2))*self.l)
            Q_L = self.U_L * SA_L* (self.T_air - self.cryogen.T_sat)
            Q_V = self.U_V*(2*np.pi*(self.d_o/2)**2 + 2*np.pi*(self.d_o/2)*self.l - SA_L)*(self.T_air-self.data["Tv_avg"])

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
        if self.Geo_l=="cylindrical":
            return self.l * (1 - self.LF)  # [m] sets vapour length
        elif self.Geo_l =="spherical" or self.Geo_l == "horizontal cylinder":
            return self.l - self.z

    @property
    def A_L(self):
        """Tank wall area in contact with the liquid"""
        if self.Geo_l == "cylindrical":
            return np.pi * self.d_o * self.l * self.LF
        elif self.Geo_l == "spherical":
            return np.pi*self.d_o*(self.z + (self.d_o - self.d_i)/2)
        elif self.Geo_l == "horizontal cylinder":
            return ((((self.d_o/2)**2)*((2*np.arccos(1-self.z/(self.d_o/2))))-np.sin(2*np.arccos(1-self.z/(self.d_o/2))))+
                              (self.d_o/2)*np.arccos(1-self.z/(self.d_o/2))*self.l)

    @property
    def A_V(self):
        """Tank wall area in contact with the vapour"""
        if self.Geo_l == "cylindrical":
            return np.pi * self.d_o * self.l * (1-self.LF)
        elif self.Geo_l == "spherical":
            return np.pi*self.d_o * (self.d_o-(self.z + (self.d_o - self.d_i)/2))
        elif self.Geo_l == "horizontal cylinder":
            return np.pi*(self.d_o/2)**2 + 2*np.pi*(self.d_o/2)*self.l - ((((self.d_o/2)**2)*((2*np.arccos(1-self.z/(self.d_o/2))))-
                    np.sin(2*np.arccos(1-self.z/(self.d_o/2))))+(self.d_o/2)*np.arccos(1-self.z/(self.d_o/2))*self.l)

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
        BL_0 = (self.Q_L_in + self.Q_b + self.Q_wi) / ((self.cryogen.h_V - self.cryogen.h_L))
        if self.Geo_v == "cylindrical":
            v_z = 4 * BL_0 / (self.cryogen.rho_V_sat * np.pi * self.d_i ** 2)
        elif self.Geo_v == "spherical":
            v_z = BL_0/(self.cryogen.rho_V_sat*np.pi*abs(2*self.z*self.d_i/2 - self.z**2))
        elif self.Geo_v == "horizontal cylinder":
            v_z = BL_0/(self.cryogen.rho_V_sat*self.l*np.sqrt(2*self.z*self.d_i/2 - self.z**2))
        return v_z

    @property
    def b_l_dot(self):
        """Returns evaporation rate in kg/s
        """
        return self.v_z * self.A_T * self.cryogen.rho_V_avg


    @property
    def Q_b(self): #needs to be changed for spherical
        if self.Geo_l =="cylindrical":
            if self.Q_b_fixed is None:
                "If Q_b_fixed is not set, calculate"
                return self.U_L * self.A_T * (self.T_air - self.cryogen.T_sat)
            else:
                return self.Q_b_fixed
        elif self.Geo_l == "spherical" or self.Geo_l == "horizontal cylinder":
            return self.Q_b_fixed

    @property
    def tau(self):
        '''Provides a conservative estimate of the 
        duration of the transient period
        of rapid vapour heating'''
        return self.l_V/self.v_z
