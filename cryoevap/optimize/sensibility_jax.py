import jax
import jax.numpy as jnp
from jax.example_libraries import optimizers
from diffrax import diffeqsolve, ODETerm, Tsit5, SaveAt, PIDController, BacksolveAdjoint, DirectAdjoint
import pandas as pd
import os
import subprocess
import matplotlib.pyplot as plt
jax.config.update("jax_enable_x64", True)

# Import polyfit function from CoolProp
folder   = '../cryoevap/cryogens/Coeffs/'

try:
    cp_V_df  = pd.read_csv(folder + 'coeffs_cpV.csv')
    k_V_df   = pd.read_csv(folder + 'coeffs_kV.csv')
    rho_V_df = pd.read_csv(folder + 'coeffs_rhoV.csv')
except FileNotFoundError:
    fitting_script = os.path.join(folder, 'Coolprop_fitting.py')
    subprocess.run(['python', fitting_script], check=True)
    cp_V_df  = pd.read_csv(folder + 'coeffs_cpV.csv')
    k_V_df   = pd.read_csv(folder + 'coeffs_kV.csv')
    rho_V_df = pd.read_csv(folder + 'coeffs_rhoV.csv')

class Sensibility_jax:
    """
    Class Sensibility_jax
    ------------------
    The Sensibility_jax class is designed to perform JAX-accelerated optimization on a cryogenic storage tank system.
    It uses automatic differentiation for gradient-based optimization of the tank's aspect ratio to minimize
    boil-off rate (BOR). The class leverages JAX for both performance and differentiability of the entire
    simulation pipeline.

    Parameters
    ----------
    tank_obj : Tank
        An instance of the Tank class representing the cryogenic tank to be optimized.
        Must have all required heat transfer properties and a defined cryogen.

    Attributes
    ----------
    tank : Tank
        The provided Tank object to be optimized.
    params : dict
        Dictionary containing all the parameters needed for the simulation, extracted from the tank object.
    time : float
        Simulation time in seconds, set during optimization.
    opt_state : OptState
        JAX optimizer state after optimization.
    optimal_aspect_ratio : float
        The computed optimal aspect ratio after optimization.
    train_loss : list
        History of the objective function values (BOR) during optimization.

    Methods
    -------
    make_params(self, tank)
        Static method to extract and format all necessary parameters from the tank object.
    cp_V_fun(T, p), k_V_fun(T, p), rho_V_fun(T, p)
        JIT-compiled functions to calculate vapor specific heat capacity, thermal conductivity,
        and density using polynomial coefficients.
    sys_isobaric_jax(t, y, args)
        JIT-compiled ODE system representing the isobaric evaporation process.
    BOR(V_L, t)
        JIT-compiled function to calculate the boil-off rate from liquid volume and time.
    evaporate(aspect_ratio)
        Simulates the evaporation process for a given aspect ratio using diffrax solvers.
    objective_function(aspect_ratio)
        Calculates the boil-off rate for a given aspect ratio, used as the optimization objective.
    optimize(verbose=True, t_final=3600*24, max_iter=100, lr=1e-1, x0=1.0)
        Performs gradient-based optimization to find the aspect ratio that minimizes boil-off rate.
        Returns the optimal aspect ratio found.
    plot_loss_history()
        Plots the convergence history of the optimization process.
    plot_surface_respone(a_array, t_final)
        Plots the response surface of the boil-off rate as a function of the aspect ratio.
    """

    def __init__(self, tank_obj):
        self.tank   = tank_obj
        self.params = self.make_params(self, self.tank)

    @staticmethod
    def make_params(self, tank):
        cryo = tank.cryogen
        return {
            "V":      jnp.array(tank.V),
            "d_i":    jnp.array(tank.d_i),
            "d_o":    jnp.array(tank.d_o),
            "LF":     jnp.array(tank.LF),
            "eta_w":  jnp.array(tank.eta_w),
            "T_air":  jnp.array(tank.T_air),
            "U_L":    jnp.array(tank.U_L),
            "U_V":    jnp.array(tank.U_V),
            "dz":     jnp.array(tank.z_grid[1] - tank.z_grid[0]),
            "z_grid": jnp.array(tank.z_grid),
            "T_sat":  jnp.array(cryo.T_sat),
            "h_L":    jnp.array(cryo.h_L),
            "h_V":    jnp.array(cryo.h_V),
            "rho_L":  jnp.array(cryo.rho_L),
            'time_interval': jnp.array(tank.time_interval),
            "cp_V_poly":  jnp.array(cp_V_df[cryo.name].values, dtype=jnp.float64),
            "k_V_poly":   jnp.array(k_V_df[cryo.name].values, dtype=jnp.float64),
            "rho_V_poly": jnp.array(rho_V_df[cryo.name].values, dtype=jnp.float64),
            "q_b_fixed": jnp.array(tank.q_b_fixed)
        }

    @staticmethod
    @jax.jit
    def cp_V_fun(T, p):
        """
        Calculate the specific heat capacity of the vapor at temperature T using polynomial coefficients p.
        
        Parameters
        ----------
        T : jnp.ndarray
            Temperature array.
        p : jnp.ndarray
            Polynomial coefficients for specific heat capacity.
        
        Returns
        -------
        jnp.ndarray
            Specific heat capacity evaluated at temperature T.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def k_V_fun(T, p):
        """
        Calculate the thermal conductivity of the vapor at temperature T using polynomial coefficients p.
        
        Parameters
        ----------
        T : jnp.ndarray
            Temperature array.
        p : jnp.ndarray
            Polynomial coefficients for thermal conductivity.
        
        Returns
        -------
        jnp.ndarray
            Thermal conductivity evaluated at temperature T.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def rho_V_fun(T, p):
        """
        Calculate the density of the vapor at temperature T using polynomial coefficients p.
        
        Parameters
        ----------
        T : jnp.ndarray
            Temperature array.
        p : jnp.ndarray
            Polynomial coefficients for density.
        
        Returns
        -------
        jnp.ndarray
            Density evaluated at temperature T.
        """
        return jnp.polyval(p, T)
   
    @staticmethod
    def q_b_fun(p):
        if p["q_b_fixed"] is None:
            "If q_b_fixed is not set, calculate"
            return p["U_L"] * (p["T_air"] - p["T_sat"])
        else:
            return p["q_b_fixed"]
    
    #####################################
    # Geometric Aspect Ratio Optimization
    #####################################

    @staticmethod
    @jax.jit
    def sys_isobaric_jax(t, y, args):
        """
        ODE system representing the isobaric evaporation process.
        
        Parameters
        ----------
        t : float
            Current time (not used in this system).
        y : jnp.ndarray
            State vector containing liquid volume and temperatures.
        args : tuple
            Tuple containing aspect ratio and parameters dictionary.
        
        Returns
        -------
        jnp.ndarray
            Derivatives of the state vector.
        """

        aspect_ratio, p = args
        d_i = ((4 * p["V"])/(jnp.pi * aspect_ratio))**(1/3)
        d_o = d_i + 0.02 

        V_L = y[0]
        T_V = y[1:]

        A_T   = jnp.pi * d_i**2 / 4
        l     = p["V"] / A_T
        l_dry = l * (1 - p['LF'])
        dz    = p["dz"] * l_dry

        k_V   = jnp.mean(Sensibility_jax.k_V_fun(T_V, p["k_V_poly"]))
        cp_V  = jnp.mean(Sensibility_jax.cp_V_fun(T_V, p["cp_V_poly"]))
        rho_V = jnp.mean(Sensibility_jax.rho_V_fun(T_V, p["rho_V_poly"]))

        A_V = jnp.pi * d_o * l * (1 - p['LF'])
        A_L = jnp.pi * d_o * l * p['LF']

        Q_Lin = p["U_L"] * A_L * (p["T_air"] - p["T_sat"])
        Q_VL  = k_V * A_T * (-3 * T_V[0] + 4 * T_V[1] - T_V[2]) / (2 * dz)
        Q_b   = Sensibility_jax.q_b_fun(p) * A_T
        Q_wi  = p["U_V"] * A_V * p["eta_w"] * (p["T_air"] - jnp.mean(T_V))

        BL_0  = (Q_Lin + Q_b + Q_wi) / (p["h_V"] - p["h_L"])
        v_z   = 4 * BL_0 / (rho_V * jnp.pi * d_i ** 2)
        v_int = v_z * (rho_V / p["rho_L"])

        alpha = k_V / (rho_V * cp_V)

        dT      = jnp.zeros_like(T_V)
        dT_dz   = (T_V[1:-1] - T_V[:-2]) / dz
        d2T_dz2 = (T_V[:-2] - 2 * T_V[1:-1] + T_V[2:]) / (dz ** 2)

        S_wall = (4 * p["U_V"] * d_o / (d_i ** 2)) * (p["T_air"] - T_V[1:-1]) * (1 - p["eta_w"])

        dT = dT.at[1:-1].set(alpha * d2T_dz2 - (v_z - v_int) * dT_dz + (alpha / k_V) * S_wall)

        dT = dT.at[0].set(0.0)
        dT = dT.at[-1].set((4 * dT[-2] - dT[-3]) / 3)

        dV = (-1 / p["rho_L"]) * (Q_Lin + Q_b + Q_wi + Q_VL) / (p["h_V"] - p["h_L"])

        return jnp.concatenate([jnp.array([dV], dtype=jnp.float64), dT])

    @staticmethod
    @jax.jit
    def BOR(V_L, t):
        """
        Calculate the Boil-Off Rate (BOR) per day from liquid volume and time.

        Parameters
        ----------
        V_L : jnp.ndarray
            Array of liquid volumes at different time steps.
        t : jnp.ndarray
            Array of time steps corresponding to the liquid volumes.
        Returns
        -------
        float
            Boil-Off Rate (BOR) per day, calculated as the fraction of initial liquid volume lost per day.
        """

        return (1.0 - V_L[-1] / V_L[0]) * (86400.0 / t[-1])

    def evaporate(self, aspect_ratio, params):
        """
        Simulates the evaporation process for a given aspect ratio using diffrax solvers.
        
        Parameters
        ----------
        aspect_ratio : float
            Aspect ratio to be used in the simulation, defined as the ratio of height to diameter.
        
        Returns
        -------
        sol : diffrax.Solution
            Solution object containing the results of the ODE integration.
        """

        VL_0 = jnp.array(params['V'] * params['LF'], dtype=jnp.float64)
        Tv_0 = jnp.ones(len(params['z_grid']), dtype=jnp.float64) * params['T_sat']
        IC   = jnp.concatenate([jnp.array([VL_0], dtype=jnp.float64), Tv_0])

        term = ODETerm(self.sys_isobaric_jax)

        sol = diffeqsolve(
            term,
            solver=Tsit5(),
            t0=0,
            t1=self.time,
            dt0=0.01,
            y0=IC,
            args=(aspect_ratio, params),
            saveat=SaveAt(ts=jnp.arange(0, self.time + 1, params['time_interval'])),
            max_steps=1000000,
            stepsize_controller=PIDController(rtol=1e-8, atol=1e-8),
            adjoint=DirectAdjoint()
        )
        return sol

    def objective_function(self, aspect_ratio, params):
        """
        Objective function to minimize the Boil-Off Rate (BOR) with respect to the aspect ratio.
        This function simulates the evaporation process for a given aspect ratio and calculates the BOR.
        
        Parameters
        ----------
        aspect_ratio : float
            Aspect ratio to be used in the simulation, defined as the ratio of height to diameter.
        
        Returns
        -------
        float
            The Boil-Off Rate (BOR) calculated from the simulation results.
        """

        sol   = self.evaporate(aspect_ratio, params)
        V_L   = sol.ys[:, 0]
        t     = sol.ts
        return self.BOR(V_L, t)

    def calculate_sensibility_param(self, param_name, aspect_ratio):
        """
        Calculate the gradient of the objective function with respect to a specific parameter.
        
        Parameters
        ----------
        param_name : str
            Name of the parameter to calculate the sensitivity for.
            Available parameters:
            - "V": Tank volume
            - "d_i": Inner diameter
            - "d_o": Outer diameter
            - "LF": Liquid fraction
            - "eta_w": Wall efficiency
            - "T_air": Air temperature
            - "U_L": Liquid heat transfer coefficient
            - "U_V": Vapor heat transfer coefficient
            - "dz": Grid spacing
            - "z_grid": Height grid array
            - "T_sat": Saturation temperature
            - "h_L": Liquid enthalpy
            - "h_V": Vapor enthalpy
            - "rho_L": Liquid density
            - "time_interval": Time interval for saving results
            - "cp_V_poly": Polynomial coefficients for vapor specific heat
            - "k_V_poly": Polynomial coefficients for vapor thermal conductivity
            - "rho_V_poly": Polynomial coefficients for vapor density
            - "q_b_fixed": Fixed bottom heat flux
        aspect_ratio : float
            The aspect ratio at which to evaluate the gradient.
        
        Returns
        -------
        float or jnp.ndarray
            The gradient of the objective function with respect to the specified parameter.
        """
        
        param_names = list(self.params.keys())
        
        if param_name not in param_names:
            raise ValueError(f"Parameter '{param_name}' not found. Available parameters: {param_names}")
        
        grads_params = jax.grad(self.objective_function, argnums=1, allow_int=True)(aspect_ratio, self.params)
        
        return grads_params[param_name]
