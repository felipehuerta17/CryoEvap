import jax
import jax.numpy as jnp
from diffrax import diffeqsolve, ODETerm, Tsit5, SaveAt, PIDController, DirectAdjoint
import pandas as pd
import os
import functools

def load_coolprop_coeffs(folder='../cryoevap/cryogens/Coeffs/'):
    """
    Loads CoolProp polynomial coefficients for density, heat capacity, and thermal conductivity.
    
    Parameters
    ----------
    folder : str
        The path to the folder containing the CSV coefficient files.
        
    Returns
    -------
    dict
        A dictionary containing the DataFrames with the polynomial coefficients.
    """
    cp_V_df  = pd.read_csv(os.path.join(folder, 'coeffs_cpV.csv'))
    k_V_df   = pd.read_csv(os.path.join(folder, 'coeffs_kV.csv'))
    rho_V_df = pd.read_csv(os.path.join(folder, 'coeffs_rhoV.csv'))
    return {'cp_V': cp_V_df, 'k_V': k_V_df, 'rho_V': rho_V_df}

class TankOptimizerJAX:
    """
    Class TankOptimizerJAX
    ----------------------
    The TankOptimizerJAX class is designed to perform JAX-accelerated optimization and sensitivity
    analysis on a cryogenic storage tank system. It uses automatic differentiation for gradient-based 
    sensitivities and provides rapid grid-search minimization of boil-off rate (BOR). The class leverages 
    JAX for both performance and differentiability of the entire simulation pipeline.

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
        Dictionary containing all the parameters needed for the simulation, formatted as JAX arrays.

    Methods
    -------
    _build_params(tank, coeffs)
        Extracts and formats all necessary parameters from the tank object.
    cp_V_fun(T, p), k_V_fun(T, p), rho_V_fun(T, p)
        JIT-compiled functions to calculate vapor properties using polynomial coefficients.
    sys_isobaric_jax(t, y, args)
        JIT-compiled system of ordinary differential equations (ODEs) describing the thermodynamic state.
    BOR(V_L, t)
        Calculates the Boil-Off Rate.
    thermal_aspect_ratio(aspect_ratio, p, T_V)
        Calculates the thermal aspect ratio of the tank.
    simulate(aspect_ratio, params, t_final)
        Runs the Diffrax ODE solver for the specified aspect ratio.
    optimize(t_final, ...)
        Searches for the optimal aspect ratio to minimize BOR using grid refinement.
    generate_surface_response_data(a_array, lf_array, t_final)
        Generates a DataFrame of BOR and thermal aspect ratio across a grid of geometries and fillings.
    calculate_sensibility_param(param_name, aspect_ratio, evap_time)
        Uses jax.grad to compute the exact derivative of BOR with respect to a target parameter.
    """
    def __init__(self, tank_obj):
        """
        Initializes the optimizer by loading the necessary polynomial properties and building the JAX parameters.
        """
        self.tank = tank_obj
        coeffs = load_coolprop_coeffs()
        self.params = self._build_params(self.tank, coeffs)

    def _build_params(self, tank, coeffs):
        """
        Extracts all scalar properties and configurations from the Tank object and places them into 
        a JAX-compatible dictionary.
        
        Parameters
        ----------
        tank : Tank
            The physical tank object.
        coeffs : dict
            The thermodynamic polynomial coefficients.
            
        Returns
        -------
        dict
            JAX dictionary of simulation parameters.
        """
        cryo = tank.cryogen
        
        # Check if tank has U_b
        if hasattr(tank, 'U_b') and tank.U_b is not None:
            U_b_val = tank.U_b
        elif hasattr(tank, 'q_b_fixed') and tank.q_b_fixed is not None:
            U_b_val = tank.q_b_fixed / (tank.T_air - cryo.T_sat)
        else:
            U_b_val = tank.U_L

        return {
            "V":      jnp.array(tank.V, dtype=jnp.float64),
            "d_i":    jnp.array(tank.d_i, dtype=jnp.float64),
            "d_o":    jnp.array(tank.d_o, dtype=jnp.float64),
            "LF":     jnp.array(tank.LF, dtype=jnp.float64),
            "eta_w":  jnp.array(tank.eta_w, dtype=jnp.float64),
            "T_air":  jnp.array(tank.T_air, dtype=jnp.float64),
            "U_L":    jnp.array(tank.U_L, dtype=jnp.float64),
            "U_V":    jnp.array(tank.U_V, dtype=jnp.float64),
            "U_b":    jnp.array(U_b_val, dtype=jnp.float64),
            "dz":     jnp.array(tank.z_grid[1] - tank.z_grid[0], dtype=jnp.float64),
            "z_grid": jnp.array(tank.z_grid, dtype=jnp.float64),
            "T_sat":  jnp.array(cryo.T_sat, dtype=jnp.float64),
            "h_L":    jnp.array(cryo.h_L, dtype=jnp.float64),
            "h_V":    jnp.array(cryo.h_V, dtype=jnp.float64),
            "rho_L":  jnp.array(cryo.rho_L, dtype=jnp.float64),
            'time_interval': float(tank.time_interval),
            "cp_V_poly":  jnp.array(coeffs['cp_V'][cryo.name].values, dtype=jnp.float64),
            "k_V_poly":   jnp.array(coeffs['k_V'][cryo.name].values, dtype=jnp.float64),
            "rho_V_poly": jnp.array(coeffs['rho_V'][cryo.name].values, dtype=jnp.float64),
        }

    @staticmethod
    @jax.jit
    def cp_V_fun(T, p):
        """
        Calculates vapor specific heat capacity using polynomial coefficients.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def k_V_fun(T, p):
        """
        Calculates vapor thermal conductivity using polynomial coefficients.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def rho_V_fun(T, p):
        """
        Calculates vapor density using polynomial coefficients.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def sys_isobaric_jax(t, y, args):
        """
        JIT-compiled system of ordinary differential equations (ODEs) describing the tank's isobaric state.
        Calculates derivatives of liquid volume and vapor temperature distribution.
        
        Parameters
        ----------
        t : float
            Current time.
        y : jax.numpy.ndarray
            State vector containing liquid volume at index 0, followed by vapor node temperatures.
        args : tuple
            Tuple containing the aspect_ratio and JAX parameter dictionary.
            
        Returns
        -------
        jax.numpy.ndarray
            Array of state derivatives (dV/dt, dT_V/dt).
        """
        aspect_ratio, p = args
        d_i = ((4 * p["V"]) / (jnp.pi * aspect_ratio)) ** (1/3)
        d_o = d_i + 0.02 

        V_L = y[0]
        T_V = y[1:]

        A_T   = jnp.pi * d_i**2 / 4
        l     = p["V"] / A_T
        l_dry = l * (1 - p['LF'])
        dz    = p["dz"] * l_dry

        k_V   = jnp.mean(TankOptimizerJAX.k_V_fun(T_V, p["k_V_poly"]))
        cp_V  = jnp.mean(TankOptimizerJAX.cp_V_fun(T_V, p["cp_V_poly"]))
        rho_V = jnp.mean(TankOptimizerJAX.rho_V_fun(T_V, p["rho_V_poly"]))

        A_V = jnp.pi * d_o * l * (1 - p['LF'])
        A_L = jnp.pi * d_o * l * p['LF']

        Q_Lin = p["U_L"] * A_L * (p["T_air"] - p["T_sat"])
        Q_VL  = k_V * A_T * (-3 * T_V[0] + 4 * T_V[1] - T_V[2]) / (2 * dz)
        Q_b   = p["U_b"] * A_T * (p["T_air"] - p["T_sat"])
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
        Calculates the Boil-Off Rate (BOR) based on initial and final liquid volumes over time.
        
        Parameters
        ----------
        V_L : jax.numpy.ndarray
            Array of liquid volumes over the evaluated time steps.
        t : jax.numpy.ndarray
            Array of time step values.
            
        Returns
        -------
        float
            Daily boil-off rate fraction.
        """
        return (1.0 - V_L[-1] / V_L[0]) * (86400.0 / t[-1])

    @staticmethod
    @jax.jit
    def thermal_aspect_ratio(aspect_ratio, p, T_V):
        """
        Calculates the thermal aspect ratio (fraction of side heat leak vs total heat leak).
        
        Parameters
        ----------
        aspect_ratio : float
            Current geometric aspect ratio being tested.
        p : dict
            JAX-compatible parameter dictionary.
        T_V : jax.numpy.ndarray
            Vapor temperature profile.
            
        Returns
        -------
        float
            The calculated thermal aspect ratio.
        """
        d_i = ((4 * p["V"]) / (jnp.pi * aspect_ratio)) ** (1/3)
        d_o = d_i + 0.02 

        A_T   = jnp.pi * d_i**2 / 4
        l     = p["V"] / A_T

        A_V = jnp.pi * d_o * l * (1 - p['LF'])
        A_L = jnp.pi * d_o * l * p['LF']

        Q_Lin = p["U_L"] * A_L * (p["T_air"] - p["T_sat"])
        Q_wi  = p["U_V"] * A_V * p["eta_w"] * (p["T_air"] - jnp.mean(T_V))
        Q_side = Q_Lin + Q_wi
        Q_b   = p["U_b"] * A_T * (p["T_air"] - p["T_sat"])

        Q_total = Q_side + Q_b
        return Q_side / Q_total

    def simulate(self, aspect_ratio, params, t_final):
        """
        Runs the Diffrax ODE solver (Tsit5) to simulate tank behavior over time.
        
        Parameters
        ----------
        aspect_ratio : float
            The aspect ratio to simulate.
        params : dict
            JAX parameter dictionary for the tank.
        t_final : float
            Simulation duration in seconds.
            
        Returns
        -------
        diffrax.Solution
            The solution object containing simulated time steps and state traces.
        """
        VL_0 = jnp.array(params['V'] * params['LF'], dtype=jnp.float64)
        Tv_0 = jnp.ones(len(params['z_grid']), dtype=jnp.float64) * params['T_sat']
        IC   = jnp.concatenate([jnp.array([VL_0], dtype=jnp.float64), Tv_0])

        term = ODETerm(TankOptimizerJAX.sys_isobaric_jax)

        sol = diffeqsolve(
            term,
            solver=Tsit5(),
            t0=0,
            t1=t_final,
            dt0=0.01,
            y0=IC,
            args=(aspect_ratio, params),
            saveat=SaveAt(ts=jnp.arange(0, t_final + 1, self.tank.time_interval)),
            max_steps=1000000,
            stepsize_controller=PIDController(rtol=1e-8, atol=1e-8),
            adjoint=DirectAdjoint()
        )
        return sol

    def _objective_bor(self, aspect_ratio, params, t_final):
        """
        Objective function calculating BOR directly from simulation for differentiation purposes.
        """
        sol = self.simulate(aspect_ratio, params, t_final)
        V_L = sol.ys[:, 0]
        t   = sol.ts
        return TankOptimizerJAX.BOR(V_L, t)

    def _objective_bor_and_tar(self, aspect_ratio, params, t_final):
        """
        Objective function wrapper that returns both BOR and Thermal Aspect Ratio.
        """
        sol = self.simulate(aspect_ratio, params, t_final)
        V_L = sol.ys[:, 0]
        T_V = sol.ys[:, 1:]
        t   = sol.ts
        bor = TankOptimizerJAX.BOR(V_L, t)
        tar = TankOptimizerJAX.thermal_aspect_ratio(aspect_ratio, params, T_V[-1])
        return bor, tar

    def optimize(self, t_final, coarse_samples=100, fine_samples=500, ar_min=0.05, ar_max=1.0, refinement_window=0.2):
        """
        Searches for the geometric aspect ratio that minimizes the Boil-Off Rate.
        Uses a two-step approach: coarse grid search followed by fine resolution refinement.
        
        Parameters
        ----------
        t_final : float
            Simulation time in seconds.
        coarse_samples : int, optional
            Number of points for the initial coarse search. Defaults to 100.
        fine_samples : int, optional
            Number of points for the refined high-resolution search. Defaults to 500.
        ar_min : float, optional
            Lower bound of geometric aspect ratio. Defaults to 0.05.
        ar_max : float, optional
            Upper bound of geometric aspect ratio. Defaults to 1.0.
        refinement_window : float, optional
            Width of the refined interval centered around the coarse optimal point. Defaults to 0.2.
            
        Returns
        -------
        tuple
            (Optimal Geometric Aspect Ratio, Optimal Thermal Aspect Ratio, Minimum BOR)
        """
        aspect_ratios = jnp.linspace(ar_min, ar_max, coarse_samples)
        
        # We define a helper, but since we don't want to recompile, we could rely on a static method or just jit over self.
        bor_values, tar_values = self._evaluate_all(aspect_ratios, self.params, t_final)
        min_idx = jnp.argmin(bor_values)
        coarse_optimal = aspect_ratios[min_idx]
        refined_min = jnp.maximum(ar_min, coarse_optimal - refinement_window/2)
        refined_max = jnp.minimum(ar_max, coarse_optimal + refinement_window/2)
        
        refined_aspect_ratios = jnp.linspace(refined_min, refined_max, fine_samples)
        refined_bor_values, refined_tar_values = self._evaluate_all(refined_aspect_ratios, self.params, t_final)
        
        refined_min_idx = jnp.argmin(refined_bor_values)
        optimal_ar = refined_aspect_ratios[refined_min_idx]
        optimal_tar = refined_tar_values[refined_min_idx]
        min_bor = refined_bor_values[refined_min_idx]

        return optimal_ar, optimal_tar, min_bor

    @functools.partial(jax.jit, static_argnums=(0,3))
    def _evaluate_all(self, a_array, params, t_final):
        """
        JIT-compiled vector mapping function that evaluates _objective_bor_and_tar over an array of aspect ratios.
        """
        return jax.vmap(lambda a: self._objective_bor_and_tar(a, params, t_final))(a_array)

    def generate_surface_response_data(self, a_array, lf_array, t_final):
        """
        Iterates over a range of Liquid Filling (LF) levels to generate a response surface.
        For each LF, it evaluates BOR and Thermal AR across the provided aspect ratios.
        
        Parameters
        ----------
        a_array : jax.numpy.ndarray
            Array of aspect ratios to evaluate.
        lf_array : jax.numpy.ndarray
            Array of liquid filling fractions to iterate across.
        t_final : float
            Simulation time in seconds.
            
        Returns
        -------
        pd.DataFrame
            DataFrame containing evaluated points with columns: LF, Geometric_AR, Thermal_AR, BOR.
        """
        results = []
        for lf in lf_array:
            lf_float = float(lf)
            local_params = self.params.copy()
            local_params['LF'] = jnp.array(lf_float, dtype=jnp.float64)
            
            bor_values, tar_values = self._evaluate_all(a_array, local_params, t_final)
            
            for a, bor, tar in zip(a_array, bor_values, tar_values):
                results.append({
                    'LF': lf_float,
                    'Geometric_AR': float(a),
                    'Thermal_AR': float(tar),
                    'BOR': float(bor)
                })
        return pd.DataFrame(results)

    def calculate_sensibility_param(self, param_name, aspect_ratio, evap_time):
        """
        Calculates the local sensitivity (gradient) of the BOR objective function with respect 
        to a specified thermodynamic/tank parameter using JAX's automatic differentiation.
        
        Parameters
        ----------
        param_name : str
            The name of the key in self.params to compute the gradient for (e.g. "U_b").
        aspect_ratio : float
            The geometric aspect ratio at which to evaluate the gradient.
        evap_time : float
            Simulation time in seconds.
            
        Returns
        -------
        float
            The sensitivity gradient (d(BOR) / d(param)).
        """
        grad_fun = jax.grad(self._objective_bor, argnums=1, allow_int=True)
        grads_params = grad_fun(aspect_ratio, self.params, evap_time)
        return grads_params[param_name]
