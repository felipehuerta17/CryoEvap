#!/usr/bin/env python3
"""
fig3a_Qavg_sweep.py
-------------------
Parametric sweep over aspect ratio (H/D) for isobaric evaporation of LN2
in a medium-scale tank (800 m³, 720 h, LF = 0.95).

For each aspect ratio the simulation is run and the time-averaged
contributions of the four heat ingress channels are extracted:
    Q_L   – wall → liquid  (lateral wall area in contact with liquid)
    Q_b   – bottom         (fixed bottom heat flux × cross-section area)
    Q_VL  – vapour → interface (Fourier conduction across vapour to interface)
    Q_Vw  – wall → interface   (η_w fraction of vapour-side wall heat)

The figure produced mirrors fig3a.png (LH2 non-isobaric reference):
  • stacked-area (left axis): fractional contributions, normalised to Q_tot
  • line (right axis, grey): absolute Q_tot in W

Outputs
-------
  Results/Data/fig3a_Qavg_nitrogen_medium.csv  – raw data per AR step
  Results/Figures/fig3a_Qavg_nitrogen_medium.png
  Results/Figures/fig3a_Qavg_nitrogen_medium.svg
"""

import os
import sys

# Always use the source-tree cryoevap (parent of Results/) instead of any
# potentially stale installed version in the venv.
_project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, _project_root)

import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # non-interactive backend for headless / WSL runs
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

# ── Physical and tank parameters (nitrogen medium tank) ───────────────────────
T_air  = 298.15   # K   – ambient temperature
P      = 100000   # Pa  – tank operating pressure (≈ 1 atm)
V_tank = 800.0    # m³  – tank volume
LF     = 0.95     # –   – initial liquid filling fraction
U_L    = 0.026    # W m⁻² K⁻¹ – overall HTC, liquid-wetted wall
U_V    = 0.026    # W m⁻² K⁻¹ – overall HTC, vapour-wetted wall
eta_w  = 0.90     # –   – wall heat partitioning fraction (→ interface)
thick  = 0.02     # –   – d_o = d_i × (1 + thick)
dz     = 0.1      # m   – nominal vertical grid spacing

evap_time     = 720 * 3600  # s  (720 h)
time_interval = 3600        # s  (record every 1 h)

# q_b_fixed: constant bottom heat flux (W m⁻²).
# Initialised per AR loop since it depends on T_sat which is fixed by P.
# Value = 1.2 × U_L × (T_air – T_sat), matching the nitrogen notebook.

# ── Aspect ratio sweep ────────────────────────────────────────────────────────
# Use a denser grid near low AR where Q_tot changes rapidly.
ar_lo  = np.linspace(0.05, 0.50, 20)   # dense near the sharp fall-off
ar_hi  = np.linspace(0.55, 2.50, 30)   # sparser toward the optimum and beyond
ar_array = np.unique(np.concatenate([ar_lo, ar_hi]))

# ── Output paths ──────────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
out_csv    = os.path.join(script_dir, 'Data',    'fig3a_Qavg_nitrogen_medium.csv')
out_png    = os.path.join(script_dir, 'Figures', 'fig3a_Qavg_nitrogen_medium.png')
out_svg    = os.path.join(script_dir, 'Figures', 'fig3a_Qavg_nitrogen_medium.svg')
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
os.makedirs(os.path.dirname(out_png), exist_ok=True)

# ── Main sweep ────────────────────────────────────────────────────────────────
records = []
t_start = time.time()
n_ar    = len(ar_array)

print(f"LN2 isobaric AR sweep: {n_ar} points, {evap_time/3600:.0f} h per run")
print(f"{'#':>4}  {'AR':>6}  {'Q_L':>10}  {'Q_b':>10}  {'Q_VL':>10}  "
      f"{'Q_Vw':>10}  {'Q_tot':>10}  {'BOR%/d':>8}  {'dt(s)':>6}")

for i, a in enumerate(ar_array):
    t0 = time.time()

    # ── Geometry ──────────────────────────────────────────────────────────────
    d_i = ((4.0 * V_tank) / (np.pi * a)) ** (1.0 / 3.0)
    d_o = d_i * (1.0 + thick)

    # ── Fresh tank ────────────────────────────────────────────────────────────
    tank = Tank(d_i, d_o, V_tank, LF)

    # ── Cryogen ───────────────────────────────────────────────────────────────
    cryogen = Cryogen(name='nitrogen')
    cryogen.set_coolprops(P)
    tank.cryogen = cryogen

    # ── Heat transfer ─────────────────────────────────────────────────────────
    q_b_fixed = 1.2 * U_L * (T_air - cryogen.T_sat)
    tank.set_HeatTransProps(U_L, U_V, T_air,
                            q_b_fixed=q_b_fixed, Q_roof=0, eta_w=eta_w)
    tank.U_roof = 0   # Neumann (insulated roof) – must be set after set_HeatTransProps

    # ── Vapour grid ───────────────────────────────────────────────────────────
    n_z = max(4, 1 + int(np.round(tank.l_V / dz, 0)))
    tank.z_grid = np.linspace(0, 1, n_z)
    tank.time_interval = time_interval

    # ── Simulate ──────────────────────────────────────────────────────────────
    tank.evaporate(evap_time)

    # ── Extract time-averaged heat ingresses ──────────────────────────────────
    Q_L_mean  = float(np.mean(tank.data['Q_L']))
    Q_b_val   = float(tank.Q_b)                  # geometry-fixed scalar
    Q_VL_mean = float(np.mean(tank.data['Q_VL']))
    Q_Vw_mean = float(np.mean(tank.data['Q_Vw']))
    Q_tot     = Q_L_mean + Q_b_val + Q_VL_mean + Q_Vw_mean
    bor       = tank.BOR() * 100.0               # convert to %/day

    records.append({
        'AR':    a,
        'Q_L':   Q_L_mean,
        'Q_b':   Q_b_val,
        'Q_VL':  Q_VL_mean,
        'Q_Vw':  Q_Vw_mean,
        'Q_tot': Q_tot,
        'BOR':   bor,
    })

    dt = time.time() - t0
    elapsed = time.time() - t_start
    print(f"{i+1:4d}  {a:6.3f}  {Q_L_mean:10.1f}  {Q_b_val:10.1f}  "
          f"{Q_VL_mean:10.1f}  {Q_Vw_mean:10.1f}  {Q_tot:10.1f}  "
          f"{bor:8.4f}  {dt:6.1f}   [total {elapsed:.0f}s]")

    # Checkpoint after every step so partial results survive a crash
    pd.DataFrame(records).to_csv(out_csv, index=False)

total_h = (time.time() - t_start) / 3600
print(f"\nSweep complete in {total_h:.2f} h.  Data written to {out_csv}")

# ── Build result arrays ───────────────────────────────────────────────────────
df    = pd.DataFrame(records)
ar    = df['AR'].values
Q_L   = df['Q_L'].values
Q_b   = df['Q_b'].values
Q_VL  = df['Q_VL'].values
Q_Vw  = df['Q_Vw'].values
Q_tot = df['Q_tot'].values

# ── Plot ──────────────────────────────────────────────────────────────────────
# Colour palette chosen to match fig3a.png
c_L   = '#7b9cd4'   # blue   – Wall → Liquid
c_b   = '#c98fb5'   # pink   – Bottom
c_VL  = '#c87070'   # salmon – Vapour → Interface
c_Vw  = '#e8c058'   # gold   – Wall → Interface

fig, ax = plt.subplots(figsize=(5.5, 5), dpi=300)

# Stack absolute heat ingresses — total height == Q_tot at every AR
ax.stackplot(ar, Q_L, Q_b, Q_VL, Q_Vw,
             labels=['Wall → Liquid', 'Bottom',
                     'Vapour → Interface', 'Wall → Interface'],
             colors=[c_L, c_b, c_VL, c_Vw])

# Black line on top of the stack showing Q_tot — minimum is visually obvious
ax.plot(ar, Q_tot, color='black', linewidth=1.8, linestyle='-', zorder=5)

ax.set_xlim(ar[0], ar[-1])
ax.set_ylim(0)
ax.set_xlabel('Aspect Ratio H/D [-]', fontsize=11, fontweight='bold')
ax.set_ylabel('Heat Ingress / W', fontsize=11, fontweight='bold')
ax.legend(loc='upper right', fontsize=8, framealpha=0.85,
          edgecolor='gray', borderpad=0.6)
ax.grid(True, alpha=0.25, linestyle='--', color='gray')

plt.tight_layout()
fig.savefig(out_png, dpi=300, bbox_inches='tight')
fig.savefig(out_svg, bbox_inches='tight')
print(f"Figures saved:\n  {out_png}\n  {out_svg}")
