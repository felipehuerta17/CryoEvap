"""
Test suite and benchmark for JAX and Diffrax accelerated parameter fitting.
Compares JAX Diffrax simulation against standard Tank.evaporate() and tests
parameter fitting on Wang & Mérida LH2 cases (Vertical and Horizontal).
"""

import os
import sys
import time
import numpy as np
import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen
from cryoevap.optimize import ParameterFittingJAX, Opti_jax


def test_vertical_simulation_comparison():
    print("\n" + "=" * 60)
    print("TEST 1: Vertical Tank Simulation (SciPy vs JAX Diffrax)")
    print("=" * 60)

    P = 101300.0
    V_tank = 4.89
    d_i = 1.676
    d_o = d_i
    LF = 0.30
    T_air = 288.15
    U_test = 0.002334
    eta_w_test = 0.50
    t_f = 5.0 * 3600.0
    dt_save = 60.0

    # 1. Standard SciPy Tank
    h2_sci = Cryogen(name='hydrogen')
    h2_sci.set_coolprops(P)
    tank_sci = Tank(d_i=d_i, d_o=d_o, V=V_tank, vapour_geometry='cylindrical',
                    liquid_geometry='cylindrical', LF=LF, head_type='flat')
    tank_sci.cryogen = h2_sci
    tank_sci.set_HeatTransProps(U_L=U_test, U_V=U_test, T_air=T_air, Q_b_fixed=None, Q_roof=0.0, eta_w=eta_w_test)
    n_z = max(3, 1 + int(np.round(tank_sci.l_V / 0.01)))
    tank_sci.z_grid = np.linspace(0.0, 1.0, n_z)
    tank_sci.time_interval = dt_save

    t0_sci = time.time()
    tank_sci.evaporate(t_f)
    time_sci = time.time() - t0_sci
    print(f"SciPy solve_ivp simulation time: {time_sci:.3f} s")

    # 2. JAX Diffrax
    fitter = ParameterFittingJAX(tank_sci)
    # Warmup / JIT compilation
    _ = fitter.simulate(np.array([U_test, eta_w_test]), t_final=t_f, time_interval=dt_save)

    t0_jax = time.time()
    res_jax = fitter.simulate(np.array([U_test, eta_w_test]), t_final=t_f, time_interval=dt_save)
    time_jax = time.time() - t0_jax
    print(f"JAX Diffrax simulation time:     {time_jax:.3f} s (Speedup: {time_sci / time_jax:.1f}x)")

    # Compare endpoints
    Tv_sci_final = tank_sci.data['Tv_avg'][-1]
    Tv_jax_final = res_jax['Tv_avg_K'].iloc[-1]
    BOG_sci_final = tank_sci.data['BOG'][-1] * 3600.0
    BOG_jax_final = res_jax['BOG_kg_h'].iloc[-1]

    print(f"Final Tv_avg: SciPy = {Tv_sci_final:.4f} K, JAX = {Tv_jax_final:.4f} K (diff = {abs(Tv_sci_final - Tv_jax_final):.4e} K)")
    print(f"Final BOG:    SciPy = {BOG_sci_final:.6f} kg/h, JAX = {BOG_jax_final:.6f} kg/h (diff = {abs(BOG_sci_final - BOG_jax_final):.4e})")

    assert np.isclose(Tv_sci_final, Tv_jax_final, atol=0.2), "Tv_avg differs too much!"
    assert np.isclose(BOG_sci_final, BOG_jax_final, atol=0.03), "BOG differs too much!"
    print("PASS: Vertical simulation comparison verified successfully!")


def test_horizontal_simulation_fast():
    print("\n" + "=" * 60)
    print("TEST 2: Horizontal Tank JAX Diffrax Simulation Benchmark")
    print("=" * 60)

    P = 101300.0
    V_tank = 4.89
    d_i = 1.676
    d_o = d_i
    L_tank = 2.217
    LF = 0.30
    T_air = 288.15
    U_test = 0.002334
    eta_w_test = 0.50
    t_f = 5.0 * 3600.0
    dt_save = 60.0

    h2 = Cryogen(name='hydrogen')
    h2.set_coolprops(P)
    tank = Tank(d_i=d_i, d_o=d_o, V=V_tank, vapour_geometry='horizontal',
                liquid_geometry='horizontal', LF=LF, L=L_tank, head_type='flat')
    tank.cryogen = h2
    tank.set_HeatTransProps(U_L=U_test, U_V=U_test, T_air=T_air, Q_b_fixed=None, Q_roof=0.0, eta_w=eta_w_test)
    n_z = max(3, 1 + int(np.round(tank.l_V / 0.01)))
    tank.z_grid = np.linspace(0.0, 1.0, n_z)
    tank.time_interval = dt_save

    fitter = ParameterFittingJAX(tank)
    # Warmup
    _ = fitter.simulate(np.array([U_test, eta_w_test]), t_final=t_f, time_interval=dt_save)

    t0_jax = time.time()
    res_jax = fitter.simulate(np.array([U_test, eta_w_test]), t_final=t_f, time_interval=dt_save)
    time_jax = time.time() - t0_jax

    print(f"JAX Diffrax 5.0h Horizontal Tank simulation time: {time_jax:.3f} s")
    print(f"Initial Tv_avg: {res_jax['Tv_avg_K'].iloc[0]:.4f} K, Final Tv_avg: {res_jax['Tv_avg_K'].iloc[-1]:.4f} K")
    print(f"Initial BOG:    {res_jax['BOG_kg_h'].iloc[0]:.6f} kg/h, Final BOG: {res_jax['BOG_kg_h'].iloc[-1]:.6f} kg/h")
    print(f"Initial LF:     {res_jax['LF'].iloc[0]:.4f}, Final LF: {res_jax['LF'].iloc[-1]:.4f}")

    assert time_jax < 1.0, "JAX horizontal simulation took too long!"
    assert res_jax['Tv_avg_K'].iloc[-1] > res_jax['Tv_avg_K'].iloc[0], "Temperature must increase!"
    print("PASS: Horizontal tank fast simulation benchmark verified successfully!")


def test_wang_merida_vertical_fitting():
    print("\n" + "=" * 60)
    print("TEST 3: Parameter Fitting on Wang & Mérida Vertical Case")
    print("=" * 60)

    excel_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'Data_wang_merida.xlsx')
    
    P = 101300.0
    V_tank = 4.89
    d_i = 1.676
    T_air = 288.15
    LF = 0.30

    h2 = Cryogen(name='hydrogen')
    h2.set_coolprops(P)
    tank = Tank(d_i=d_i, d_o=d_i, V=V_tank, vapour_geometry='cylindrical',
                liquid_geometry='cylindrical', LF=LF, head_type='flat')
    tank.cryogen = h2
    tank.set_HeatTransProps(U_L=0.002, U_V=0.002, T_air=T_air, Q_b_fixed=None, Q_roof=0.0, eta_w=0.5)

    fitter = ParameterFittingJAX(tank)
    fitter.from_wang_merida(excel_path, case='vertical', t_max_h=5.0)

    theta0 = np.array([0.002334, 0.50])
    print(f"Initial parameters: U = {theta0[0]:.6f}, eta_w = {theta0[1]:.4f}, Initial J = {fitter.objective(theta0):.4f}")

    t0 = time.time()
    opt_res = fitter.fit_least_squares(
        theta0,
        bounds=([0.0005, 0.30], [0.006, 0.99]),
        verbose=1,
        max_nfev=20
    )
    total_time = time.time() - t0

    th_opt = opt_res['x']
    print(f"Vertical fitting elapsed time: {total_time:.2f} s")
    print(f"Fitted parameters:  U = {th_opt[0]:.6f}, eta_w = {th_opt[1]:.4f}, Final J = {opt_res['cost']:.4f}")
    assert opt_res['success'], "Optimization did not converge!"
    assert opt_res['cost'] < 4.0, "Objective should be minimized below initial value!"
    print("PASS: Wang & Mérida vertical parameter fitting verified successfully!")


def test_wang_merida_horizontal_fitting():
    print("\n" + "=" * 60)
    print("TEST 4: Parameter Fitting on Wang & Mérida Horizontal Case")
    print("=" * 60)

    excel_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'Data_wang_merida.xlsx')
    
    P = 101300.0
    V_tank = 4.89
    d_i = 1.676
    L_tank = 2.217
    T_air = 288.15
    LF = 0.30

    h2 = Cryogen(name='hydrogen')
    h2.set_coolprops(P)
    tank = Tank(d_i=d_i, d_o=d_i, V=V_tank, vapour_geometry='horizontal',
                liquid_geometry='horizontal', LF=LF, L=L_tank, head_type='flat')
    tank.cryogen = h2
    tank.set_HeatTransProps(U_L=0.002, U_V=0.002, T_air=T_air, Q_b_fixed=None, Q_roof=0.0, eta_w=0.5)

    fitter = ParameterFittingJAX(tank)
    fitter.from_wang_merida(excel_path, case='horizontal', t_max_h=5.0)

    theta0 = np.array([0.002334, 0.50])
    print(f"Initial parameters: U = {theta0[0]:.6f}, eta_w = {theta0[1]:.4f}, Initial J = {fitter.objective(theta0):.4f}")

    t0 = time.time()
    opt_res = fitter.fit_least_squares(
        theta0,
        bounds=([0.0005, 0.30], [0.006, 0.99]),
        verbose=1,
        max_nfev=20
    )
    total_time = time.time() - t0

    th_opt = opt_res['x']
    print(f"Horizontal fitting elapsed time: {total_time:.2f} s")
    print(f"Fitted parameters:  U = {th_opt[0]:.6f}, eta_w = {th_opt[1]:.4f}, Final J = {opt_res['cost']:.4f}")
    assert opt_res['success'], "Optimization did not converge!"
    print("PASS: Wang & Mérida horizontal parameter fitting verified successfully!")


if __name__ == '__main__':
    test_vertical_simulation_comparison()
    test_horizontal_simulation_fast()
    test_wang_merida_vertical_fitting()
    test_wang_merida_horizontal_fitting()
    print("\n" + "=" * 60)
    print("ALL JAX PARAMETER FITTING TESTS PASSED SUCCESSFULLY!")
    print("=" * 60)
