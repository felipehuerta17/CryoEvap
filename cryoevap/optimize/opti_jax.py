"""
Class Opti_jax: JAX-accelerated optimization of cryogenic storage tanks.
Optimizes aspect ratio and parameters to minimize Boil-Off Rate (BOR).
"""

import os
import subprocess
import matplotlib.pyplot as plt
import jax
import jax.numpy as jnp
from jax.example_libraries import optimizers
from diffrax import diffeqsolve, ODETerm, Tsit5, SaveAt, PIDController, BacksolveAdjoint
import pandas as pd

jax.config.update("jax_enable_x64", True)

# Load polynomial coefficients for cryogenic fluid properties
current_dir = os.path.dirname(os.path.abspath(__file__))
coeffs_dir = os.path.join(current_dir, '..', 'cryogens', 'Coeffs')

try:
    cp_V_df  = pd.read_csv(os.path.join(coeffs_dir, 'coeffs_cpV.csv'), index_col=0)
    k_V_df   = pd.read_csv(os.path.join(coeffs_dir, 'coeffs_kV.csv'), index_col=0)
    rho_V_df = pd.read_csv(os.path.join(coeffs_dir, 'coeffs_rhoV.csv'), index_col=0)
except Exception:
    fitting_script = os.path.join(coeffs_dir, 'Coolprop_fitting.py')
    if os.path.exists(fitting_script):
        subprocess.run(['python', fitting_script], check=True)
        cp_V_df  = pd.read_csv(os.path.join(coeffs_dir, 'coeffs_cpV.csv'), index_col=0)
        k_V_df   = pd.read_csv(os.path.join(coeffs_dir, 'coeffs_kV.csv'), index_col=0)
        rho_V_df = pd.read_csv(os.path.join(coeffs_dir, 'coeffs_rhoV.csv'), index_col=0)
    else:
        cp_V_df, k_V_df, rho_V_df = None, None, None


class Opti_jax:
    """
    The Opti_jax class performs JAX-accelerated aspect ratio and thermal aspect ratio
    optimization on a cryogenic storage tank system to minimize Boil-Off Rate (BOR).
    """

    def __init__(self, tank_obj):
        self.tank = tank_obj
        self.params = self.make_params(self.tank)
        self.time = 3600.0 * 24.0
        self.opt_state = None
        self.optimal_aspect_ratio = None
        self.optimal_thermal_aspect_ratio = None
        self.train_loss = []

    @staticmethod
    def make_params(tank):
        cryo = tank.cryogen
        name = cryo.name.lower()
        if cp_V_df is not None and name not in cp_V_df.columns:
            name = 'hydrogen'

        cp_poly = jnp.array(cp_V_df[name].values, dtype=jnp.float64) if cp_V_df is not None else jnp.zeros(7)
        k_poly = jnp.array(k_V_df[name].values, dtype=jnp.float64) if k_V_df is not None else jnp.zeros(5)
        rho_poly = jnp.array(rho_V_df[name].values, dtype=jnp.float64) if rho_V_df is not None else jnp.zeros(7)

        return {
            "V": tank.V,
            "d_i": tank.d_i,
            "d_o": tank.d_o,
            "LF": tank.LF,
            "eta_w": getattr(tank, 'eta_w', 0.5),
            "T_air": getattr(tank, 'T_air', 288.15),
            "U_L": getattr(tank, 'U_L', 0.002),
            "U_V": getattr(tank, 'U_V', 0.002),
            "dz": tank.z_grid[1] - tank.z_grid[0],
            "z_grid": tank.z_grid,
            "T_sat": cryo.T_sat,
            "h_L": cryo.h_L,
            "h_V": cryo.h_V,
            "rho_L": cryo.rho_L,
            'time_interval': getattr(tank, 'time_interval', 3600.0),
            "cp_V_poly": cp_poly,
            "k_V_poly": k_poly,
            "rho_V_poly": rho_poly,
            "q_b_fixed": getattr(tank, 'q_b_fixed', None)
        }

    @staticmethod
    @jax.jit
    def cp_V_fun(T, p):
        return jnp.polyval(p["cp_V_poly"], T)

    @staticmethod
    @jax.jit
    def k_V_fun(T, p):
        return jnp.polyval(p["k_V_poly"], T)

    @staticmethod
    @jax.jit
    def rho_V_fun(T, p):
        return jnp.polyval(p["rho_V_poly"], T)

    @staticmethod
    def q_b_fun(p):
        if p["q_b_fixed"] is None:
            return p["U_L"] * (p["T_air"] - p["T_sat"])
        else:
            return p["q_b_fixed"]

    @staticmethod
    @jax.jit
    def sys_isobaric_jax(t, y, args):
        aspect_ratio, p = args
        d_i = ((4 * p["V"]) / (jnp.pi * aspect_ratio))**(1/3)
        d_o = d_i + 0.02

        V_L = y[0]
        T_V = y[1:]

        A_T = jnp.pi * d_i**2 / 4.0
        l = p["V"] / A_T
        l_dry = l * (1.0 - p['LF'])
        dz = p["dz"] * l_dry

        k_V = jnp.mean(Opti_jax.k_V_fun(T_V, p))
        cp_V = jnp.mean(Opti_jax.cp_V_fun(T_V, p))
        rho_V = jnp.mean(Opti_jax.rho_V_fun(T_V, p))

        A_V = jnp.pi * d_o * l * (1.0 - p['LF'])
        A_L = jnp.pi * d_o * l * p['LF']

        Q_Lin = p["U_L"] * A_L * (p["T_air"] - p["T_sat"])
        Q_VL = k_V * A_T * (-3.0 * T_V[0] + 4.0 * T_V[1] - T_V[2]) / (2.0 * dz)
        Q_b = Opti_jax.q_b_fun(p) * A_T
        Q_wi = p["U_V"] * A_V * p["eta_w"] * (p["T_air"] - jnp.mean(T_V))

        BL_0 = (Q_Lin + Q_b + Q_wi) / (p["h_V"] - p["h_L"])
        v_z = 4.0 * BL_0 / (rho_V * jnp.pi * d_i**2)
        v_int = v_z * (rho_V / p["rho_L"])

        alpha = k_V / (rho_V * cp_V)

        dT_dz = (T_V[1:-1] - T_V[:-2]) / dz
        d2T_dz2 = (T_V[:-2] - 2.0 * T_V[1:-1] + T_V[2:]) / (dz**2)

        S_wall = (4.0 * p["U_V"] * d_o / (d_i**2)) * (p["T_air"] - T_V[1:-1]) * (1.0 - p["eta_w"])

        dT = jnp.zeros_like(T_V)
        dT = dT.at[1:-1].set(alpha * d2T_dz2 - (v_z - v_int) * dT_dz + (alpha / k_V) * S_wall)
        dT = dT.at[0].set(0.0)
        dT = dT.at[-1].set((4.0 * dT[-2] - dT[-3]) / 3.0)

        dV = (-1.0 / p["rho_L"]) * (Q_Lin + Q_b + Q_wi + Q_VL) / (p["h_V"] - p["h_L"])

        return jnp.concatenate([jnp.array([dV], dtype=jnp.float64), dT])

    @staticmethod
    @jax.jit
    def BOR(V_L, t):
        return (1.0 - V_L[-1] / V_L[0]) * (86400.0 / t[-1])

    def evaporate(self, aspect_ratio):
        VL_0 = jnp.array(self.params['V'] * self.params['LF'], dtype=jnp.float64)
        Tv_0 = jnp.ones(len(self.params['z_grid']), dtype=jnp.float64) * self.params['T_sat']
        IC = jnp.concatenate([jnp.array([VL_0], dtype=jnp.float64), Tv_0])

        term = ODETerm(self.sys_isobaric_jax)

        sol = diffeqsolve(
            term,
            solver=Tsit5(),
            t0=0.0,
            t1=self.time,
            dt0=0.01,
            y0=IC,
            args=(aspect_ratio, self.params),
            saveat=SaveAt(ts=jnp.arange(0.0, self.time + 1.0, self.params['time_interval'])),
            max_steps=1000000,
            stepsize_controller=PIDController(rtol=1e-8, atol=1e-8),
            adjoint=BacksolveAdjoint()
        )
        return sol

    def objective_function(self, aspect_ratio):
        a_eff = jnp.exp(aspect_ratio)
        sol = self.evaporate(a_eff)
        V_L = sol.ys[:, 0]
        t = sol.ts
        return self.BOR(V_L, t)

    def optimize(self, verbose=True, t_final=3600*24, max_iter=100, lr=1e-1, x0=1.0, lr_2=1e-3, iter_decay_lr=100):
        train_loss = []
        self.time = t_final

        opt_init, opt_update, get_params = optimizers.adam(lr)
        opt_state = opt_init(jnp.array(x0, dtype=jnp.float64))

        for i in range(max_iter):
            params = get_params(opt_state)
            loss_val, grads = jax.value_and_grad(self.objective_function)(params)
            opt_state = opt_update(i, grads, opt_state)
            train_loss.append(loss_val)

            if verbose:
                print(f"Iter {i}, Loss {loss_val:.6e}, Param {float(jnp.exp(params)):.4g}")

            if i == iter_decay_lr:
                opt_init, opt_update, get_params = optimizers.adam(lr_2)
                opt_state = opt_init(get_params(opt_state))

        self.opt_state = opt_state
        self.optimal_aspect_ratio = float(jnp.exp(get_params(opt_state)))
        self.train_loss = train_loss
        return self.optimal_aspect_ratio

    def optimize_grid_with_refinement(self, verbose=True, t_final=3600*24,
                                      coarse_samples=100, fine_samples=100,
                                      aspect_ratio_min=0.2, aspect_ratio_max=3.0,
                                      refinement_window=0.2):
        self.time = t_final
        if verbose:
            print(f"Phase 1: Coarse grid search with {coarse_samples} samples...")

        aspect_ratios = jnp.linspace(aspect_ratio_min, aspect_ratio_max, coarse_samples)
        bor_values = jax.vmap(lambda a: self.objective_function(jnp.log(a)))(aspect_ratios)

        min_idx = jnp.argmin(bor_values)
        coarse_optimal = float(aspect_ratios[min_idx])

        if verbose:
            print(f"Coarse search optimal: {coarse_optimal:.6f}, BOR: {bor_values[min_idx]:.6e}")

        refined_min = max(aspect_ratio_min, coarse_optimal - refinement_window/2.0)
        refined_max = min(aspect_ratio_max, coarse_optimal + refinement_window/2.0)

        refined_aspect_ratios = jnp.linspace(refined_min, refined_max, fine_samples)
        refined_bor_values = jax.vmap(lambda a: self.objective_function(jnp.log(a)))(refined_aspect_ratios)

        refined_min_idx = jnp.argmin(refined_bor_values)
        optimal_aspect_ratio = float(refined_aspect_ratios[refined_min_idx])
        min_bor = float(refined_bor_values[refined_min_idx])

        self.optimal_aspect_ratio = optimal_aspect_ratio
        self.train_loss = bor_values.tolist() + refined_bor_values.tolist()

        if verbose:
            print(f"Refined optimal aspect ratio: {optimal_aspect_ratio:.6f}, Minimum BOR: {min_bor:.6e}")

        return optimal_aspect_ratio
