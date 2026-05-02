"""
Sensitivity analysis: external convective heat transfer coefficient h_env.

Runs the LF=0.95 ammonia tank simulation with the wall model at three values
of h_env (baseline, -50 %, +100 %) and saves transient and 3-week pickle
datasets so Figs. 2 and 3 can be re-plotted with the sensitivity bands.

Only the wall Robin BC is altered: U_L, U_V and U_b are kept at their
nominal values because the perlite resistance dominates and they shift by
< 0.5 % across the swept range.
"""

import os
import pickle
import sys

import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from cryoevap.cryogens import Cryogen
from cryoevap.storage_tanks import Tank


def build_tank(h_env_value):
    """Construct the LF=0.95 large ammonia tank for a given h_env."""
    e = 10.6972 * 2 - 10.24 * 2          # wall thickness / m
    V_tank = 88023.6952                  # tank volume / m^3
    a = 0.5                              # aspect ratio
    d_i = ((4 * V_tank) / (np.pi * a)) ** (1 / 3)
    d_o = d_i + e
    LF = 0.95

    tank = Tank(d_i, d_o, V_tank, LF)

    # Heat transfer / wall properties (kept at the paper baseline).
    U_L = 8.86344e-02
    U_V = 8.86344e-02
    U_b = 8.80129e-02
    h_L = 135.08
    k_w, cp_w, rho_w = 0.0411, 900.0, 60.0
    T_air = 5.3286 + 273.15
    P = 101325

    tank.set_HeatTransProps(
        U_L, U_V, T_air,
        Q_b_fixed=None, Q_roof=0.0, eta_w=0.70,
        k_w=k_w, rho_w=rho_w, cp_w=cp_w, h_L=h_L, T_init=True,
    )
    tank.U_b = U_b

    cryogen = Cryogen(name="ammonia")
    cryogen.set_coolprops(P)
    tank.cryogen = cryogen

    # Spatial grid (matches the published runs).
    dz, dr = 0.05, 0.01
    n_z = 1 + int(np.round(tank.l_V / dz, 0))
    n_r = 1 + int(np.round((tank.d_o - tank.d_i) / (2 * dr), 0))
    tank.z_grid = np.linspace(0, 1, n_z)
    tank.r_grid = np.linspace(0, 1, n_r)
    tank.U_roof = 0.0

    # Environmental BC: only h_env is varied.
    T_env_avg = 5.3286 + 273.15
    T_range = 15.0
    tank.set_EnvironmentalProps(
        T_avg_day=T_env_avg, T_range_day=T_range, h_env=h_env_value, p_anual=None,
    )

    return tank


def run_case(label, h_env_value, data_dir):
    """Run the transient (24 h) and 3-week simulations for one h_env value."""
    print(f"\n=== Case '{label}' : h_env = {h_env_value:.4f} W m^-2 K^-1 ===")

    # ---- Transient ( t = 5*tau = 24 h, dt = 1728 s ) ----
    tank = build_tank(h_env_value)
    tank.time_interval = 1728
    transient_time = 86400
    print(f"  - transient: integrating {transient_time/3600:.1f} h ...")
    tank.evaporate(transient_time)
    out = os.path.join(data_dir, f"large_tank_data_LF95_transient_{label}.pkl")
    with open(out, "wb") as fh:
        pickle.dump(tank.data, fh)
    print(f"    saved {out}")

    # ---- 3 weeks ( dt = 3600 s ) ----
    tank = build_tank(h_env_value)
    tank.time_interval = 3600
    long_time = 3600 * 24 * 7 * 3
    print(f"  - long term: integrating {long_time/3600/24:.1f} d ...")
    tank.evaporate(long_time)
    out = os.path.join(data_dir, f"large_tank_data_LF95_3weeks_{label}.pkl")
    with open(out, "wb") as fh:
        pickle.dump(tank.data, fh)
    print(f"    saved {out}")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(here, "Data")
    os.makedirs(data_dir, exist_ok=True)

    h_base = 14.839  # baseline used in the paper
    cases = {
        "h_low":  0.5 * h_base,
        "h_base": 1.0 * h_base,
        "h_high": 2.0 * h_base,
    }
    for label, h_value in cases.items():
        run_case(label, h_value, data_dir)


if __name__ == "__main__":
    main()
