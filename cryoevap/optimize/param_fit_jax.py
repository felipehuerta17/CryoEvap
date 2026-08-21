"""
Parameter fitting module accelerated with JAX and Diffrax for cryogenic storage tanks.
Supports cylindrical (vertical) and horizontal tanks (flat and hemispherical heads).
"""

import os
import subprocess
import time
from copy import deepcopy

import jax
import jax.numpy as jnp
from diffrax import diffeqsolve, ODETerm, Tsit5, SaveAt, PIDController
from equinox import EquinoxRuntimeError
import numpy as np
import pandas as pd
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

# Enable 64-bit precision for high accuracy ODE integration
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


class ParameterFittingJAX:
    """
    High-performance parameter fitting and simulation class using JAX and Diffrax.
    
    Fits heat transfer parameters (overall heat transfer coefficient U and wall heat partitioning eta_w)
    against experimental or benchmark time series data (average vapor temperature and BOG rate).
    """

    def __init__(self, tank_obj=None, **kwargs):
        if tank_obj is not None:
            self.tank = tank_obj
            self.params = self.make_params(self.tank)
        else:
            self.tank = None
            self.params = kwargs

        # Experimental datasets
        self.Tv_exp = None
        self.BOG_exp = None
        self.scales = {'Tv_avg': 1.0, 'BOG': 1.0}
        self.weights = {'Tv_avg': 1.0, 'BOG': 1.0}
        self.optimization_history = []
        self.last_results = None

    @staticmethod
    def _solve_h_horiz_np(V_L, d_i, L, head_type, V_tank):
        R_i = d_i / 2.0
        is_hemi = 1.0 if head_type == 1 else 0.0
        LF = np.clip(V_L / V_tank, 1e-5, 1.0 - 1e-5)
        h = d_i * LF
        for _ in range(8):
            h_c = np.clip(h, 1e-6, d_i - 1e-6)
            arg = np.clip((R_i - h_c) / R_i, -1.0, 1.0)
            sqrt_t = np.sqrt(np.maximum(2.0 * R_i * h_c - h_c**2, 1e-12))
            V_c = L * (R_i**2 * np.arccos(arg) - (R_i - h_c) * sqrt_t)
            V_s = (1.0 / 3.0) * np.pi * h_c**2 * (3.0 * R_i - h_c) * is_hemi
            V_calc = V_c + V_s
            A_c = 2.0 * L * sqrt_t
            A_s = np.pi * (2.0 * R_i * h_c - h_c**2) * is_hemi
            A_T = np.maximum(A_c + A_s, 1e-6)
            h = np.clip(h_c - (V_calc - V_L) / A_T, 1e-6, d_i - 1e-6)
        return float(h), float(A_T)

    @staticmethod
    def make_params(tank):
        """Extract simulation parameters from a Tank object into a JAX-ready dictionary."""
        cryo = tank.cryogen
        name = cryo.name.lower()
        if cp_V_df is not None and name not in cp_V_df.columns:
            name = 'hydrogen'

        d_i = float(tank.d_i)
        d_o = float(tank.d_o)
        V = float(tank.V)
        LF = float(tank.LF)
        L = getattr(tank, 'L', None)
        if L is None:
            L = V / (np.pi * d_i**2 / 4.0)
        else:
            L = float(L)

        head_type_str = getattr(tank, 'head_type', 'flat')
        head_type = 1 if head_type_str == 'hemispherical' else 0

        geom_v = getattr(tank, 'Geo_v', 'cylindrical')
        geom_l = getattr(tank, 'Geo_l', 'cylindrical')

        dz_grid = float(tank.z_grid[1] - tank.z_grid[0])
        n_nodes = len(tank.z_grid)
        z_grid = jnp.linspace(0.0, 1.0, n_nodes, dtype=jnp.float64)

        cp_poly = jnp.array(cp_V_df[name].values, dtype=jnp.float64) if cp_V_df is not None else jnp.zeros(7)
        k_poly = jnp.array(k_V_df[name].values, dtype=jnp.float64) if k_V_df is not None else jnp.zeros(5)
        rho_poly = jnp.array(rho_V_df[name].values, dtype=jnp.float64) if rho_V_df is not None else jnp.zeros(5)

        # Precompute initial liquid height and area for horizontal tank
        VL_0 = V * LF
        z_L_0, A_T_0 = ParameterFittingJAX._solve_h_horiz_np(VL_0, d_i, L, head_type, V)

        return {
            'V': V,
            'd_i': d_i,
            'd_o': d_o,
            'L': L,
            'LF': LF,
            'V_L_0': float(VL_0),
            'z_L_0': float(z_L_0),
            'A_T_0': float(A_T_0),
            'eta_w': float(getattr(tank, 'eta_w', 0.5)),
            'T_air': float(getattr(tank, 'T_air', 288.15)),
            'U_L': float(getattr(tank, 'U_L', 0.002)),
            'U_V': float(getattr(tank, 'U_V', 0.002)),
            'Q_b_fixed': float(tank.Q_b_fixed) if getattr(tank, 'Q_b_fixed', None) is not None else -1.0,
            'T_sat': float(cryo.T_sat),
            'h_L': float(cryo.h_L),
            'h_V': float(cryo.h_V),
            'rho_L': float(cryo.rho_L),
            'rho_V_sat': float(cryo.rho_V_sat) if hasattr(cryo, 'rho_V_sat') else float(cryo.rho_V),
            'dz_grid': dz_grid,
            'n_nodes': n_nodes,
            'z_grid': z_grid,
            'geom_v': 1 if geom_v == 'horizontal' else 0,
            'geom_l': 1 if geom_l == 'horizontal' else 0,
            'head_type': head_type,
            'time_interval': float(getattr(tank, 'time_interval', 60.0)),
            'cp_V_poly': cp_poly,
            'k_V_poly': k_poly,
            'rho_V_poly': rho_poly
        }

    # Property evaluation functions
    @staticmethod
    def cp_V_fun(T, p):
        return jnp.polyval(p['cp_V_poly'], T)

    @staticmethod
    def k_V_fun(T, p):
        return jnp.polyval(p['k_V_poly'], T)

    @staticmethod
    def rho_V_fun(T, p):
        return jnp.polyval(p['rho_V_poly'], 1.0 / T)

    # Fast 2-step unrolled Newton solver for liquid level
    @staticmethod
    def solve_h_horizontal(V_L, p):
        d_i = p['d_i']
        L = p['L']
        R_i = d_i / 2.0
        is_hemi = p['head_type']
        
        # Initial guess from precomputed z_L_0
        h_0 = jnp.clip(p['z_L_0'] + (V_L - p['V_L_0']) / p['A_T_0'], 1e-6, d_i - 1e-6)

        def step(h_c):
            arg = jnp.clip((R_i - h_c) / R_i, -1.0, 1.0)
            sqrt_t = jnp.sqrt(jnp.maximum(2.0 * R_i * h_c - h_c**2, 1e-12))
            V_c = L * (R_i**2 * jnp.arccos(arg) - (R_i - h_c) * sqrt_t)
            V_s = (1.0 / 3.0) * jnp.pi * h_c**2 * (3.0 * R_i - h_c) * is_hemi
            V_calc = V_c + V_s
            A_c = 2.0 * L * sqrt_t
            A_s = jnp.pi * (2.0 * R_i * h_c - h_c**2) * is_hemi
            A_T = jnp.maximum(A_c + A_s, 1e-6)
            return jnp.clip(h_c - (V_calc - V_L) / A_T, 1e-6, d_i - 1e-6)

        h_1 = step(h_0)
        h_2 = step(h_1)
        return h_2

    # Isobaric ODE right-hand side in JAX
    @staticmethod
    def sys_isobaric_jax(t, y, args):
        p, U, eta_w = args
        V_L = y[0]
        T_V = y[1:]
        z_grid = p['z_grid']
        
        is_horizontal = (p['geom_l'] == 1)

        # ----------------------------------------------------
        # 1. Vertical Cylindrical Tank
        # ----------------------------------------------------
        def vertical_branch():
            A_T = jnp.pi * p['d_i']**2 / 4.0
            l_tank = p['V'] / A_T
            LF = jnp.clip(V_L / p['V'], 1e-6, 1.0)
            L_dry = jnp.maximum(l_tank * (1.0 - LF), 1e-4)
            dz = p['dz_grid'] * L_dry

            k_V = jnp.mean(ParameterFittingJAX.k_V_fun(T_V, p))
            cp_V = jnp.mean(ParameterFittingJAX.cp_V_fun(T_V, p))
            rho_V = jnp.mean(ParameterFittingJAX.rho_V_fun(T_V, p))
            alpha = k_V / (rho_V * cp_V)

            Tv_avg = jnp.mean(T_V)

            A_L = jnp.pi * p['d_o'] * l_tank * LF
            A_V = jnp.pi * p['d_o'] * l_tank * (1.0 - LF) + jnp.pi * p['d_o']**2 / 4.0

            Q_Lin = U * A_L * (p['T_air'] - p['T_sat'])
            Q_b = jnp.where(p['Q_b_fixed'] >= 0.0, p['Q_b_fixed'], U * A_T * (p['T_air'] - p['T_sat']))
            Q_wi = U * A_V * eta_w * (p['T_air'] - Tv_avg)
            Q_VL = k_V * A_T * (-3.0 * T_V[0] + 4.0 * T_V[1] - T_V[2]) / (2.0 * dz)

            BL_0 = (Q_Lin + Q_b + Q_wi) / (p['h_V'] - p['h_L'])
            v_z = 4.0 * BL_0 / (rho_V * jnp.pi * p['d_i']**2)
            v_int = v_z * (rho_V / p['rho_L'])

            dT_dz = (T_V[1:-1] - T_V[:-2]) / dz
            d2T_dz2 = (T_V[:-2] - 2.0 * T_V[1:-1] + T_V[2:]) / (dz**2)
            S_wall = (4.0 * U * p['d_o'] / (p['d_i']**2)) * (p['T_air'] - T_V[1:-1]) * (1.0 - eta_w)

            dT = jnp.zeros_like(T_V)
            dT_inner = alpha * d2T_dz2 - (v_z - v_int) * dT_dz + (alpha / k_V) * S_wall
            dT = dT.at[1:-1].set(dT_inner)
            dT = dT.at[0].set(0.0)
            dT = dT.at[-1].set((4.0 * dT[-2] - dT[-3]) / 3.0)

            dV = (-1.0 / p['rho_L']) * (Q_Lin + Q_b + Q_wi + Q_VL) / (p['h_V'] - p['h_L'])
            return jnp.concatenate([jnp.array([dV], dtype=jnp.float64), dT])

        # ----------------------------------------------------
        # 2. Horizontal Tank
        # ----------------------------------------------------
        def horizontal_branch():
            z_L = ParameterFittingJAX.solve_h_horizontal(V_L, p)
            L_dry = jnp.maximum(p['d_i'] - z_L, 1e-4)
            dz = p['dz_grid'] * L_dry

            # Physical height of each node in vapor space
            h_grid = z_grid * L_dry + z_L
            half_chord_grid = jnp.sqrt(jnp.maximum(p['d_i'] * h_grid - h_grid**2, 0.0))
            weight_cylinder = 2.0 * p['L'] * half_chord_grid
            weight_heads = jnp.pi * half_chord_grid**2 * p['head_type']
            weights = jnp.maximum(weight_cylinder + weight_heads, 1e-6)

            # Average vapor temperature with weights
            Tv_avg = jnp.sum(T_V * weights) / jnp.sum(weights)

            k_V = jnp.mean(ParameterFittingJAX.k_V_fun(T_V, p))
            cp_V = jnp.mean(ParameterFittingJAX.cp_V_fun(T_V, p))
            rho_V = jnp.mean(ParameterFittingJAX.rho_V_fun(T_V, p))
            alpha = k_V / (rho_V * cp_V)

            # Interfacial area at liquid surface
            half_chord_int = jnp.sqrt(jnp.maximum(p['d_i'] * z_L - z_L**2, 0.0))
            A_T = 2.0 * p['L'] * half_chord_int + jnp.pi * half_chord_int**2 * p['head_type']
            A_T = jnp.maximum(A_T, 1e-6)

            # Wetted and vapor wall areas
            R_o = p['d_o'] / 2.0
            wall_thick = (p['d_o'] - p['d_i']) / 2.0
            h_o = jnp.clip(z_L + wall_thick, 0.0, 2.0 * R_o)
            theta_o = jnp.arccos(jnp.clip((R_o - h_o) / R_o, -1.0, 1.0))
            A_lateral_wet = p['d_o'] * p['L'] * theta_o

            sqrt_term_o = jnp.sqrt(jnp.maximum(2.0 * R_o * h_o - h_o**2, 0.0))
            A_heads_flat = 2.0 * (R_o**2 * theta_o - (R_o - h_o) * sqrt_term_o)
            A_heads_hemi = 2.0 * jnp.pi * R_o * h_o
            A_heads_wet = jnp.where(p['head_type'] == 1, A_heads_hemi, A_heads_flat)

            A_total_flat = jnp.pi * p['d_o'] * p['L'] + 2.0 * jnp.pi * R_o**2
            A_total_hemi = jnp.pi * p['d_o'] * p['L'] + 4.0 * jnp.pi * R_o**2
            A_total = jnp.where(p['head_type'] == 1, A_total_hemi, A_total_flat)

            A_L = A_lateral_wet + A_heads_wet
            A_V = jnp.maximum(A_total - A_L, 1e-6)

            Q_Lin = U * A_L * (p['T_air'] - p['T_sat'])
            Q_b = jnp.where(p['Q_b_fixed'] >= 0.0, p['Q_b_fixed'], 0.0)
            Q_wi = U * A_V * eta_w * (p['T_air'] - Tv_avg)
            Q_VL = k_V * A_T * (-3.0 * T_V[0] + 4.0 * T_V[1] - T_V[2]) / (2.0 * dz)

            BL_0 = (Q_Lin + Q_b + Q_wi) / (p['h_V'] - p['h_L'])
            v_z0 = BL_0 / (p['rho_V_sat'] * A_T)
            v_int = v_z0 * (rho_V / p['rho_L'])

            # Velocity profile at intermediate nodes
            z_mid = h_grid[1:-1]
            sqrt_int = jnp.sqrt(jnp.maximum(p['d_i'] * z_L - z_L**2, 0.0))
            sqrt_mid = jnp.sqrt(jnp.maximum(p['d_i'] * z_mid - z_mid**2, 1e-7))
            v_z = v_z0 * (sqrt_int / sqrt_mid)

            # Area variation and geometric convection
            A_cs = 2.0 * p['L'] * sqrt_mid
            dAcs_dz = p['L'] * (p['d_i'] - 2.0 * z_mid) / sqrt_mid
            geom_term = alpha * (dAcs_dz / jnp.maximum(A_cs, 1e-12)) * (T_V[1:-1] - T_V[:-2]) / dz

            dA_dz = -(4.0 * sqrt_mid + p['d_i'] * p['L'] / sqrt_mid)
            dV_dz = -2.0 * p['L'] * sqrt_mid
            S_wall = U * (p['T_air'] - T_V[1:-1]) * (1.0 - eta_w) * (dA_dz / dV_dz)

            dT_dz = (T_V[1:-1] - T_V[:-2]) / dz
            d2T_dz2 = (T_V[:-2] - 2.0 * T_V[1:-1] + T_V[2:]) / (dz**2)

            dT = jnp.zeros_like(T_V)
            dT_inner = alpha * d2T_dz2 + geom_term - (v_z - v_int) * dT_dz + (alpha / k_V) * S_wall
            dT = dT.at[1:-1].set(dT_inner)
            dT = dT.at[0].set(0.0)
            dT = dT.at[-1].set((4.0 * dT[-2] - dT[-3]) / 3.0)

            dV = (-1.0 / p['rho_L']) * (Q_Lin + Q_b + Q_wi + Q_VL) / (p['h_V'] - p['h_L'])
            return jnp.concatenate([jnp.array([dV], dtype=jnp.float64), dT])

        return jax.lax.cond(is_horizontal, horizontal_branch, vertical_branch)

    # Fast simulation driver with Diffrax
    def simulate(self, theta, t_final=5.0*3600.0, time_interval=60.0):
        """
        Simulate tank evaporation for parameter vector theta = [U, eta_w].
        
        Returns a DataFrame with columns: ['time_h', 'Tv_avg_K', 'BOG_kg_h', 'V_L', 'LF']
        """
        # Keep U and eta_w as JAX arrays (not Python floats) so Equinox's
        # filter_jit treats them as dynamic (traced) arguments. Otherwise they
        # are treated as static, forcing a full XLA recompilation on every
        # single call with a different theta (extremely slow during fitting).
        U = jnp.asarray(theta[0], dtype=jnp.float64)
        eta_w = jnp.asarray(theta[1], dtype=jnp.float64)

        p = self.params
        VL_0 = jnp.array(p['V'] * p['LF'], dtype=jnp.float64)
        Tv_0 = jnp.ones(p['n_nodes'], dtype=jnp.float64) * p['T_sat']
        IC = jnp.concatenate([jnp.array([VL_0]), Tv_0])

        t_eval = jnp.arange(0.0, t_final + 1.0, time_interval)
        term = ODETerm(self.sys_isobaric_jax)

        sol = diffeqsolve(
            term,
            solver=Tsit5(),
            t0=0.0,
            t1=t_final,
            dt0=0.1,
            y0=IC,
            args=(p, U, eta_w),
            saveat=SaveAt(ts=t_eval),
            max_steps=2000000,
            stepsize_controller=PIDController(rtol=1e-5, atol=1e-5)
        )

        t_arr = np.array(sol.ts)
        y_arr = np.array(sol.ys)
        VL_arr = y_arr[:, 0]
        Tv_grid = y_arr[:, 1:]

        # Post-process Tv_avg, rho_V_avg, heat transfers, and BOG matching Tank._reconstruct
        n_nodes = p['n_nodes']
        z_grid = np.linspace(0.0, 1.0, n_nodes)
        is_horiz = (p['geom_l'] == 1)

        if is_horiz:
            # Vectorized liquid height
            z_L_arr = p['z_L_0'] + (VL_arr - p['V_L_0']) / p['A_T_0']
            L_dry_arr = np.maximum(p['d_i'] - z_L_arr, 1e-4) # shape (N_t,)
            
            # h_nodes shape: (N_t, N_z)
            h_nodes = z_L_arr[:, None] + z_grid[None, :] * L_dry_arr[:, None]
            half_c = np.sqrt(np.maximum(p['d_i'] * h_nodes - h_nodes**2, 0.0))
            w = half_c * 2.0 * p['L'] + np.pi * half_c**2 * p['head_type']
            w = np.maximum(w, 1e-6)

            # Tv_avg and rho_V_avg vectorized over time
            Tv_avg_arr = np.sum(Tv_grid * w, axis=1) / np.sum(w, axis=1)
            rho_nodes = np.polyval(p['rho_V_poly'], 1.0 / np.maximum(Tv_grid, 1.0))
            rho_V_avg_arr = np.sum(rho_nodes * w, axis=1) / np.sum(w, axis=1)

            # Areas
            R_o = p['d_o'] / 2.0
            wall_t = (p['d_o'] - p['d_i']) / 2.0
            h_o = np.clip(z_L_arr + wall_t, 0.0, 2.0 * R_o)
            th_o = np.arccos(np.clip((R_o - h_o) / R_o, -1.0, 1.0))
            A_lat_wet = p['d_o'] * p['L'] * th_o
            if p['head_type'] == 1:
                A_h_wet = 2.0 * np.pi * R_o * h_o
                A_tot = np.pi * p['d_o'] * p['L'] + 4.0 * np.pi * R_o**2
            else:
                sqrt_t = np.sqrt(np.maximum(2.0 * R_o * h_o - h_o**2, 0.0))
                A_h_wet = 2.0 * (R_o**2 * th_o - (R_o - h_o) * sqrt_t)
                A_tot = np.pi * p['d_o'] * p['L'] + 2.0 * np.pi * R_o**2
            A_L_arr = A_lat_wet + A_h_wet
            A_V_arr = np.maximum(A_tot - A_L_arr, 1e-6)

            half_c_int = np.sqrt(np.maximum(p['d_i'] * z_L_arr - z_L_arr**2, 0.0))
            A_T_arr = 2.0 * p['L'] * half_c_int + np.pi * half_c_int**2 * p['head_type']

            dz_arr = (z_grid[1] - z_grid[0]) * L_dry_arr
            k_v_nodes = np.polyval(p['k_V_poly'], Tv_grid)
            k_v_avg = np.mean(k_v_nodes, axis=1)
            dTdz_arr = (-3.0 * Tv_grid[:, 0] + 4.0 * Tv_grid[:, 1] - Tv_grid[:, 2]) / (2.0 * dz_arr)
            Q_VL_arr = k_v_avg * A_T_arr * dTdz_arr

            Q_L_arr = U * A_L_arr * (p['T_air'] - p['T_sat'])
            Q_V_arr = U * A_V_arr * (p['T_air'] - Tv_avg_arr)
            Q_b_val = 0.0 if p['Q_b_fixed'] < 0.0 else p['Q_b_fixed']
        else:
            Tv_avg_arr = np.mean(Tv_grid, axis=1)
            rho_nodes = np.polyval(p['rho_V_poly'], 1.0 / np.maximum(Tv_grid, 1.0))
            rho_V_avg_arr = np.mean(rho_nodes, axis=1)

            A_T = np.pi * p['d_i']**2 / 4.0
            l_tank = p['V'] / A_T
            l_L_arr = VL_arr / A_T
            l_V_arr = np.maximum(l_tank - l_L_arr, 1e-4)

            A_L_arr = np.pi * p['d_o'] * l_L_arr
            A_V_arr = np.pi * p['d_o'] * l_V_arr + np.pi * p['d_o']**2 / 4.0

            dz_arr = (z_grid[1] - z_grid[0]) * l_V_arr
            k_v_nodes = np.polyval(p['k_V_poly'], Tv_grid)
            k_v_avg = np.mean(k_v_nodes, axis=1)
            dTdz_arr = (-3.0 * Tv_grid[:, 0] + 4.0 * Tv_grid[:, 1] - Tv_grid[:, 2]) / (2.0 * dz_arr)
            Q_VL_arr = k_v_avg * A_T * dTdz_arr

            Q_L_arr = U * A_L_arr * (p['T_air'] - p['T_sat'])
            Q_V_arr = U * A_V_arr * (p['T_air'] - Tv_avg_arr)
            Q_b_val = p['Q_b_fixed'] if p['Q_b_fixed'] >= 0.0 else (U * A_T * (p['T_air'] - p['T_sat']))

        Q_Vw_arr = Q_V_arr * eta_w

        # Derivatives with 2nd order finite differences matching Tank.dydt
        def dydt(t, y):
            dy = np.zeros(len(t))
            dt = t[1] - t[0]
            dy[0] = (-3.0 * y[0] + 4.0 * y[1] - y[2]) / (2.0 * dt)
            dy[1:-1] = (y[2:] - y[:-2]) / (2.0 * dt)
            dy[-1] = (3.0 * y[-1] - 4.0 * y[-2] + y[-3]) / (2.0 * dt)
            return dy

        dVL_dt = dydt(t_arr, VL_arr)
        drhoV_dt = dydt(t_arr, rho_V_avg_arr)

        dH_LV = p['h_V'] - p['h_L']
        B_L = (Q_L_arr + Q_b_val + Q_VL_arr + Q_Vw_arr) / dH_LV

        BOG_kg_s = B_L + rho_V_avg_arr * dVL_dt - (p['V'] - VL_arr) * drhoV_dt
        BOG_kg_h = np.maximum(BOG_kg_s * 3600.0, 0.0)

        results = pd.DataFrame({
            'time_h': t_arr / 3600.0,
            'Tv_avg_K': Tv_avg_arr,
            'BOG_kg_h': BOG_kg_h,
            'V_L': VL_arr,
            'LF': VL_arr / p['V']
        })
        self.last_results = results
        return results

    # Data loading and scaling methods
    def set_experimental_data(self, Tv_exp=None, BOG_exp=None, scale_method='std'):
        """Set experimental dataframes for Tv_avg and BOG."""
        if Tv_exp is not None:
            self.Tv_exp = Tv_exp.copy()
            self.scales['Tv_avg'] = self.calculate_scale(self.Tv_exp['value'], method=scale_method)

        if BOG_exp is not None:
            self.BOG_exp = BOG_exp.copy()
            self.scales['BOG'] = self.calculate_scale(self.BOG_exp['value'], method=scale_method)

    def from_wang_merida(self, excel_path, case='vertical', t_max_h=5.0, scale_method='std'):
        """
        Load Wang and Merida experimental data from Excel sheet (Figure 12 LH2 benchmark).
        
        case : 'vertical' or 'horizontal'
        """
        sheets = pd.read_excel(excel_path, sheet_name=None)

        def prep_curve(df, t_col, val_col):
            c = df[[t_col, val_col]].copy().dropna()
            c.columns = ['time_h', 'value']
            c['time_h'] = pd.to_numeric(c['time_h'], errors='coerce')
            c['value'] = pd.to_numeric(c['value'], errors='coerce')
            c = c.dropna().sort_values('time_h')
            if t_max_h is not None:
                c = c[c['time_h'] <= t_max_h]
            return c.groupby('time_h', as_index=False).mean().reset_index(drop=True)

        if case == 'vertical':
            Tv = prep_curve(sheets['12a_T'], 't_Tv_vert_h', 'Tv_vert_K')
            BOG = prep_curve(sheets['12e_BOG'], 't_vert_h', 'BOG_vert_kg_h')
        else:
            Tv = prep_curve(sheets['12a_T'], 't_Tv_horiz_h', 'Tv_horiz_K')
            BOG = prep_curve(sheets['12e_BOG'], 't_horiz_h', 'BOG_horiz_kg_h')

        self.set_experimental_data(Tv_exp=Tv, BOG_exp=BOG, scale_method=scale_method)
        print(f"Loaded Wang & Merida ({case}): {len(Tv)} Tv points, {len(BOG)} BOG points.")
        return Tv, BOG

    @staticmethod
    def calculate_scale(values, method='std'):
        values = np.asarray(values, dtype=float)
        if method == 'std':
            scale = float(np.std(values))
        elif method == 'max':
            scale = float(np.max(np.abs(values)))
        elif method == 'mean':
            scale = float(np.mean(np.abs(values)))
        else:
            scale = 1.0
        return scale if np.isfinite(scale) and scale > 1e-12 else 1.0

    # Residuals & Objective computation
    def calculate_residuals(self, theta, t_final=None, time_interval=60.0,
                             steady_state=False, steady_state_frac=0.2):
        """
        Calculate weighted normalized residual vector for theta = [U, eta_w].

        Parameters
        ----------
        steady_state : bool
            If True, only the last `steady_state_frac` fraction of the
            experimental time window (by time) is used to build the
            residuals. Useful when the reduced-order model's transient
            response does not match the reference case (e.g. Wang & Merida)
            but the asymptotic steady-state Tv_avg / BOG values are the
            quantities of interest.
        steady_state_frac : float
            Fraction (0-1) of the experimental time span, measured from the
            end, considered "steady-state". Only used if steady_state=True.
        """
        if t_final is None:
            max_t = 0.0
            if self.Tv_exp is not None:
                max_t = max(max_t, self.Tv_exp['time_h'].max() * 3600.0)
            if self.BOG_exp is not None:
                max_t = max(max_t, self.BOG_exp['time_h'].max() * 3600.0)
            t_final = max(max_t, 3600.0)

        try:
            results = self.simulate(theta, t_final=t_final, time_interval=time_interval)
        except EquinoxRuntimeError:
            # The ODE solver diverged / exceeded max_steps for this particular
            # (U, eta_w) combination. This can happen transiently while the
            # optimizer probes points near the parameter bounds. Instead of
            # crashing the whole fit, return a heavily penalized residual so
            # least_squares treats this as a very poor fit and steers away.
            n_res = 0
            if self.Tv_exp is not None:
                n_res += len(self.Tv_exp)
            if self.BOG_exp is not None:
                n_res += len(self.BOG_exp)
            return np.full(n_res, 1e6), {'diverged': True}

        t_model = results['time_h'].values

        residuals_list = []
        details = {}

        def steady_window(df):
            if not steady_state:
                return df
            t_cut = (1.0 - steady_state_frac) * df['time_h'].max()
            window = df[df['time_h'] >= t_cut]
            # Guard against an empty/singleton window (e.g. very sparse data)
            return window if len(window) >= 1 else df.iloc[[-1]]

        if self.Tv_exp is not None:
            Tv_fit = steady_window(self.Tv_exp)
            t_exp_tv = Tv_fit['time_h'].values
            val_exp_tv = Tv_fit['value'].values
            val_model_tv = np.interp(t_exp_tv, t_model, results['Tv_avg_K'].values)
            res_tv = (val_model_tv - val_exp_tv) / self.scales['Tv_avg']
            w_res_tv = res_tv * np.sqrt(self.weights['Tv_avg'] / len(res_tv))
            residuals_list.append(w_res_tv)
            details['Tv_model'] = val_model_tv
            details['residuals_Tv'] = res_tv
            details['t_exp_Tv_fit'] = t_exp_tv

        if self.BOG_exp is not None:
            BOG_fit = steady_window(self.BOG_exp)
            t_exp_bog = BOG_fit['time_h'].values
            val_exp_bog = BOG_fit['value'].values
            val_model_bog = np.interp(t_exp_bog, t_model, results['BOG_kg_h'].values)
            res_bog = (val_model_bog - val_exp_bog) / self.scales['BOG']
            w_res_bog = res_bog * np.sqrt(self.weights['BOG'] / len(res_bog))
            residuals_list.append(w_res_bog)
            details['BOG_model'] = val_model_bog
            details['residuals_BOG'] = res_bog
            details['t_exp_BOG_fit'] = t_exp_bog

        res_vector = np.concatenate(residuals_list)
        return res_vector, details

    def objective(self, theta, t_final=None, time_interval=60.0,
                  steady_state=False, steady_state_frac=0.2):
        """Evaluate total least squares objective J(theta)."""
        res_vector, _ = self.calculate_residuals(
            theta, t_final=t_final, time_interval=time_interval,
            steady_state=steady_state, steady_state_frac=steady_state_frac)
        return float(np.sum(res_vector**2))

    # Fast optimization routines
    def fit_least_squares(self, theta0, bounds=([1e-4, 0.0], [0.05, 1.0]),
                          t_final=None, time_interval=60.0, verbose=1,
                          ftol=1e-4, xtol=1e-4, gtol=1e-4, max_nfev=30,
                          steady_state=False, steady_state_frac=0.2):
        """
        Fit parameters theta = [U, eta_w] using Scipy least_squares with JAX Diffrax engine.

        If steady_state=True, only the last `steady_state_frac` fraction of
        the experimental time window is used to build the residuals, so the
        fit targets the asymptotic steady-state Tv_avg / BOG values instead
        of the full (possibly model-mismatched) transient.
        """
        self.optimization_history = []
        start_time = time.time()

        def residual_fn(th):
            res_vec, _ = self.calculate_residuals(
                th, t_final=t_final, time_interval=time_interval,
                steady_state=steady_state, steady_state_frac=steady_state_frac)
            cost = float(np.sum(res_vec**2))
            self.optimization_history.append({'U': float(th[0]), 'eta_w': float(th[1]), 'cost': cost})
            if verbose >= 2:
                print(f"Eval {len(self.optimization_history):02d}: U={th[0]:.8f}, eta_w={th[1]:.5f}, J={cost:.6f}")
            return res_vec

        res = least_squares(
            residual_fn,
            theta0,
            bounds=bounds,
            method='trf',
            ftol=ftol,
            xtol=xtol,
            gtol=gtol,
            max_nfev=max_nfev,
            verbose=verbose
        )

        elapsed = time.time() - start_time
        th_opt = res.x
        cost_opt = 2.0 * res.cost

        if verbose >= 1:
            print(f"Optimization finished in {elapsed:.2f} s ({res.nfev} evaluations).")
            print(f"Optimal U     = {th_opt[0]:.10f} W/(m^2 K)")
            print(f"Optimal eta_w = {th_opt[1]:.6f}")
            print(f"Optimal J     = {cost_opt:.8f}")

        return {
            'x': th_opt,
            'cost': cost_opt,
            'success': res.success,
            'message': res.message,
            'nfev': res.nfev,
            'elapsed_time': elapsed,
            'history': self.optimization_history
        }

    # Plotting routines
    def plot_fit(self, theta, title='Parameter Fit Comparison', save_path=None,
                 steady_state=False, steady_state_frac=0.2):
        """Plot model results against experimental data for Tv_avg and BOG.

        If steady_state=True, the region of the experimental time window that
        was actually used to build the fit residuals (the trailing
        `steady_state_frac` of the time span) is shaded, and points outside
        that window are shown in a lighter, non-fitted style.
        """
        res_vector, details = self.calculate_residuals(
            theta, steady_state=steady_state, steady_state_frac=steady_state_frac)
        results = self.last_results

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5), dpi=150)

        def shade_steady_window(ax, df):
            if not steady_state or df is None:
                return
            t_cut = (1.0 - steady_state_frac) * df['time_h'].max()
            ax.axvspan(t_cut, df['time_h'].max(), color='grey', alpha=0.15,
                       label='Steady-state fit window')

        # 1. Temperature
        if self.Tv_exp is not None:
            shade_steady_window(ax1, self.Tv_exp)
            ax1.plot(self.Tv_exp['time_h'], self.Tv_exp['value'], 'o', color='lightcoral',
                     label='Experimental (excluded)', markersize=6)
            t_fit = details.get('t_exp_Tv_fit')
            if t_fit is not None:
                mask = self.Tv_exp['time_h'].isin(t_fit)
                ax1.plot(self.Tv_exp['time_h'][mask], self.Tv_exp['value'][mask], 'ro',
                         label='Experimental (fitted)', markersize=6)
        ax1.plot(results['time_h'], results['Tv_avg_K'], 'b-', label='JAX Model', linewidth=2)
        ax1.set_xlabel('Time (h)')
        ax1.set_ylabel('Average Vapor Temp (K)')
        ax1.set_title('Vapor Temperature')
        ax1.grid(True, linestyle='--', alpha=0.6)
        ax1.legend()

        # 2. BOG
        if self.BOG_exp is not None:
            shade_steady_window(ax2, self.BOG_exp)
            ax2.plot(self.BOG_exp['time_h'], self.BOG_exp['value'], 'o', color='lightgreen',
                     label='Experimental (excluded)', markersize=6)
            t_fit = details.get('t_exp_BOG_fit')
            if t_fit is not None:
                mask = self.BOG_exp['time_h'].isin(t_fit)
                ax2.plot(self.BOG_exp['time_h'][mask], self.BOG_exp['value'][mask], 'go',
                         label='Experimental (fitted)', markersize=6)
        ax2.plot(results['time_h'], results['BOG_kg_h'], 'g-', label='JAX Model', linewidth=2)
        ax2.set_xlabel('Time (h)')
        ax2.set_ylabel('BOG Rate (kg/h)')
        ax2.set_title('Boil-Off Gas (BOG)')
        ax2.grid(True, linestyle='--', alpha=0.6)
        ax2.legend()

        fig.suptitle(f"{title} (U={theta[0]:.6f} W/m^2K, eta_w={theta[1]:.4f})")
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
        return fig, (ax1, ax2)
