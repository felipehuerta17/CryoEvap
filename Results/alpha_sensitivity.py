"""
Sensitivity analysis: wall thermal diffusivity alpha_w = k_w / (rho_w * cp_w).

Runs the LF=0.95 ammonia tank simulation with the wall model at three values
of k_w (baseline, -8.8 %, +8.8 %), yielding alpha_w variations of +/-8.8 %.
Both the wall PDE and the Robin BC use k_w, so varying k_w captures the full
physical effect, including the change in overall heat-transfer coefficient
(perlite resistance dominates; U scales approximately as k_w).

k values correspond to the range of perlite thermal conductivity across the
operating temperature window (240-300 K) reported by the Perlite Institute.

U_L, U_V are rescaled proportionally to k_w (perlite resistance ~ 99 % of
the total), while h_env and all geometry are kept at the paper baseline.
"""

import os
import pickle
import sys

import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from cryoevap.cryogens import Cryogen
from cryoevap.storage_tanks import Tank

# Baseline paper parameters
K_W_BASE  = 0.0411          # W m^-1 K^-1 (baseline thermal conductivity)
CP_W      = 900.0            # J kg^-1 K^-1 (kept constant)
RHO_W     = 60.0             # kg m^-3 (kept constant)
U_L_BASE  = 8.86344e-02      # W m^-2 K^-1
U_V_BASE  = 8.86344e-02      # W m^-2 K^-1
U_B_BASE  = 8.80129e-02      # W m^-2 K^-1
H_L       = 135.08           # W m^-2 K^-1 (liquid-wall film coefficient)
H_ENV     = 14.839           # W m^-2 K^-1 (environmental convection, fixed)
T_AIR     = 5.3286 + 273.15  # K
T_RANGE   = 15.0             # K (diurnal temperature swing)
P         = 101325           # Pa

# Sensitivity: +-8.8 % variation in k_w (and hence alpha_w and U)
DELTA = 0.088


def build_tank(k_w_value):
    """Construct the LF=0.95 large ammonia tank for a given k_w."""
    e      = 10.6972 * 2 - 10.24 * 2
    V_tank = 88023.6952
    a      = 0.5
    d_i    = ((4 * V_tank) / (np.pi * a)) ** (1 / 3)
    d_o    = d_i + e
    LF     = 0.95

    tank = Tank(d_i, d_o, V_tank, LF)

    # Scale U proportionally to k_w (perlite resistance dominates)
    scale = k_w_value / K_W_BASE
    U_L   = U_L_BASE * scale
    U_V   = U_V_BASE * scale
    U_b   = U_B_BASE * scale

    tank.set_HeatTransProps(
        U_L, U_V, T_AIR,
        Q_b_fixed=None, Q_roof=0.0, eta_w=0.70,
        k_w=k_w_value, rho_w=RHO_W, cp_w=CP_W, h_L=H_L, T_init=True,
    )
    tank.U_b = U_b

    cryogen = Cryogen(name="ammonia")
    cryogen.set_coolprops(P)
    tank.cryogen = cryogen

    dz, dr = 0.05, 0.01
    n_z = 1 + int(np.round(tank.l_V / dz, 0))
    n_r = 1 + int(np.round((tank.d_o - tank.d_i) / (2 * dr), 0))
    tank.z_grid = np.linspace(0, 1, n_z)
    tank.r_grid = np.linspace(0, 1, n_r)
    tank.U_roof = 0.0

    tank.set_EnvironmentalProps(
        T_avg_day=T_AIR, T_range_day=T_RANGE, h_env=H_ENV, p_anual=None,
    )

    return tank


def run_case(label, k_w_value, data_dir):
    """Run transient (24 h) and 3-week simulations for one k_w value."""
    alpha = k_w_value / (RHO_W * CP_W)
    print(f"\n=== Case '{label}' : k_w = {k_w_value:.5f} W m^-1 K^-1 "
          f"(alpha = {alpha*1e6:.4f} x10^-6 m^2/s) ===")

    # ---- Transient ( t = 5*tau ≈ 24 h, dt = 1728 s ) ----
    tank = build_tank(k_w_value)
    tank.time_interval = 1728
    transient_time = 86400
    print(f"  transient: integrating {transient_time/3600:.1f} h ...")
    tank.evaporate(transient_time)
    out = os.path.join(data_dir, f"large_tank_data_LF95_transient_alpha_{label}.pkl")
    with open(out, "wb") as fh:
        pickle.dump(tank.data, fh)
    print(f"  saved {out}")

    # ---- 3 weeks ( dt = 3600 s ) ----
    tank = build_tank(k_w_value)
    tank.time_interval = 3600
    long_time = 3600 * 24 * 7 * 3
    print(f"  long-term: integrating {long_time/3600/24:.1f} d ...")
    tank.evaporate(long_time)
    out = os.path.join(data_dir, f"large_tank_data_LF95_3weeks_alpha_{label}.pkl")
    with open(out, "wb") as fh:
        pickle.dump(tank.data, fh)
    print(f"  saved {out}")


def main():
    here     = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(here, "Data")
    os.makedirs(data_dir, exist_ok=True)

    cases = {
        "alpha_low":  K_W_BASE * (1 - DELTA),
        "alpha_base": K_W_BASE,
        "alpha_high": K_W_BASE * (1 + DELTA),
    }
    for label, k_w in cases.items():
        run_case(label, k_w, data_dir)


if __name__ == "__main__":
    main()
