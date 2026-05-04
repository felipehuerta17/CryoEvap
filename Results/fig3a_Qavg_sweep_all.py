#!/usr/bin/env python3
"""
fig3a_Qavg_sweep_all.py
-----------------------
Parametric AR sweep for six cryogen/thermal-aspect-ratio scenarios:

  Cryogen   V (m³)   U_L=U_V (W/m²K)  q_b multiplier  Label
  --------  -------  ---------------  --------------  -----
  nitrogen    800       0.026           1.205           rq07
  nitrogen    800       0.026           6.5845          rq03
  hydrogen  50 000     3.73e-3          1.205           rq07
  hydrogen  50 000     3.73e-3          6.5845          rq03
  methane   50 000     0.038            1.205           rq07
  methane   50 000     0.038            6.5845          rq03

For each scenario the AR sweep and figure are produced independently so
partial results survive a crash (CSV checkpointed after every AR step).
Completed cases (CSV row count matches expected) are skipped on re-run.

Outputs (per case):
  Results/Data/fig3a_Qavg_{cryogen}_{label}.csv
  Results/Figures/fig3a_Qavg_{cryogen}_{label}.png
  Results/Figures/fig3a_Qavg_{cryogen}_{label}.svg
"""

import os
import sys

_project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, _project_root)

import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

# ── Shared constants ─────────────────────────────────────────────────────────
T_air         = 298.15   # K
P             = 100000   # Pa
eta_w         = 0.90
thick         = 0.02     # d_o = d_i × (1 + thick)
dz            = 0.1      # m
evap_time     = 720 * 3600   # s
time_interval = 3600         # s

# ── Case table ───────────────────────────────────────────────────────────────
# q_b_mult: q_b_fixed = q_b_mult × U_L × (T_air − T_sat)
# rq03 (η≈0.30): bottom-dominated  → optimal AR shifts right
# rq07 (η≈0.70): side-dominated    → optimal AR shifts left
CASES = [
    dict(cryogen='nitrogen', V_tank=800.0,   U_L=0.026,    LF=0.95, q_b_mult=1.205,  label='rq07'),
    dict(cryogen='nitrogen', V_tank=800.0,   U_L=0.026,    LF=0.95, q_b_mult=6.5845, label='rq03'),
    dict(cryogen='hydrogen', V_tank=50000.0, U_L=3.73e-3,  LF=0.90, q_b_mult=1.205,  label='rq07'),
    dict(cryogen='hydrogen', V_tank=50000.0, U_L=3.73e-3,  LF=0.90, q_b_mult=6.5845, label='rq03'),
    dict(cryogen='methane',  V_tank=50000.0, U_L=0.038,    LF=0.90, q_b_mult=1.205,  label='rq07'),
    dict(cryogen='methane',  V_tank=50000.0, U_L=0.038,    LF=0.90, q_b_mult=6.5845, label='rq03'),
]

# ── AR grid (same for all cases) ─────────────────────────────────────────────
# Denser below AR=0.5 where Q_tot changes most rapidly; extends to 3.5 to
# capture the rq03 minimum which sits at AR ≈ q_b_mult/4 ≈ 1.65
ar_lo   = np.linspace(0.05, 0.50, 20)
ar_mid  = np.linspace(0.55, 1.50, 20)
ar_hi   = np.linspace(1.60, 3.50, 20)
ar_array = np.unique(np.concatenate([ar_lo, ar_mid, ar_hi]))

# ── Output directories ────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
dir_data   = os.path.join(script_dir, 'Data')
dir_fig    = os.path.join(script_dir, 'Figures')
os.makedirs(dir_data, exist_ok=True)
os.makedirs(dir_fig,  exist_ok=True)

# ── Colour palette (consistent across all figures) ───────────────────────────
c_L   = '#7b9cd4'   # blue   – Wall → Liquid
c_b   = '#c98fb5'   # pink   – Bottom
c_VL  = '#c87070'   # salmon – Vapour → Interface
c_Vw  = '#e8c058'   # gold   – Wall → Interface

# ── Helper: run one sweep ────────────────────────────────────────────────────
def run_sweep(case):
    cryogen_name = case['cryogen']
    V_tank       = case['V_tank']
    U_L          = case['U_L']
    LF           = case['LF']
    q_b_mult     = case['q_b_mult']
    label        = case['label']

    tag      = f"{cryogen_name}_{label}"
    out_csv  = os.path.join(dir_data, f'fig3a_Qavg_{tag}.csv')
    out_png  = os.path.join(dir_fig,  f'fig3a_Qavg_{tag}.png')
    out_svg  = os.path.join(dir_fig,  f'fig3a_Qavg_{tag}.svg')

    n_ar = len(ar_array)

    # Skip if CSV is already complete
    if os.path.exists(out_csv):
        df_existing = pd.read_csv(out_csv)
        if len(df_existing) >= n_ar:
            print(f"\n[{tag}] CSV already complete ({len(df_existing)} rows) — skipping sweep.")
            return df_existing

    # Resume from checkpoint if partial CSV exists
    if os.path.exists(out_csv):
        df_existing = pd.read_csv(out_csv)
        done_ars = set(df_existing['AR'].values)
        records  = df_existing.to_dict('records')
        print(f"\n[{tag}] Resuming from {len(records)} existing rows.")
    else:
        done_ars = set()
        records  = []

    print(f"\n{'='*70}")
    print(f"  Case: {tag}  |  V={V_tank:.0f} m³  U_L={U_L:.4g}  q_b_mult={q_b_mult}  LF={LF}")
    print(f"  {n_ar} AR points, {evap_time/3600:.0f} h per run")
    print(f"{'='*70}")
    print(f"{'#':>4}  {'AR':>6}  {'Q_L':>10}  {'Q_b':>10}  {'Q_VL':>10}  "
          f"{'Q_Vw':>10}  {'Q_tot':>10}  {'BOR%/d':>8}  {'dt(s)':>6}")

    t_start = time.time()

    for i, a in enumerate(ar_array):
        # Round AR to avoid floating-point duplicate matches
        if round(a, 8) in {round(x, 8) for x in done_ars}:
            continue

        t0 = time.time()

        d_i  = ((4.0 * V_tank) / (np.pi * a)) ** (1.0 / 3.0)
        d_o  = d_i * (1.0 + thick)

        tank = Tank(d_i, d_o, V_tank, LF)

        cryogen = Cryogen(name=cryogen_name)
        cryogen.set_coolprops(P)
        tank.cryogen = cryogen

        q_b_fixed = q_b_mult * U_L * (T_air - cryogen.T_sat)
        tank.set_HeatTransProps(U_L, U_L, T_air,
                                q_b_fixed=q_b_fixed, Q_roof=0, eta_w=eta_w)
        tank.U_roof = 0   # Neumann BC — must be set after set_HeatTransProps

        n_z = max(4, 1 + int(np.round(tank.l_V / dz, 0)))
        tank.z_grid      = np.linspace(0, 1, n_z)
        tank.time_interval = time_interval

        tank.evaporate(evap_time)

        Q_L_mean  = float(np.mean(tank.data['Q_L']))
        Q_b_val   = float(tank.Q_b)
        Q_VL_mean = float(np.mean(tank.data['Q_VL']))
        Q_Vw_mean = float(np.mean(tank.data['Q_Vw']))
        Q_tot     = Q_L_mean + Q_b_val + Q_VL_mean + Q_Vw_mean
        bor       = tank.BOR() * 100.0

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

        # Checkpoint after every step
        pd.DataFrame(records).to_csv(out_csv, index=False)

    total_h = (time.time() - t_start) / 3600
    print(f"\n[{tag}] Sweep complete in {total_h:.2f} h.  CSV: {out_csv}")
    return pd.DataFrame(records)


# ── Helper: plot one case ─────────────────────────────────────────────────────
def plot_case(df, case):
    cryogen_name = case['cryogen']
    label        = case['label']
    q_b_mult     = case['q_b_mult']
    tag          = f"{cryogen_name}_{label}"
    out_png      = os.path.join(dir_fig, f'fig3a_Qavg_{tag}.png')
    out_svg      = os.path.join(dir_fig, f'fig3a_Qavg_{tag}.svg')

    # Sort by AR (in case of out-of-order checkpoints)
    df = df.sort_values('AR').reset_index(drop=True)

    ar    = df['AR'].values
    Q_L   = df['Q_L'].values
    Q_b   = df['Q_b'].values
    Q_VL  = df['Q_VL'].values
    Q_Vw  = df['Q_Vw'].values
    Q_tot = df['Q_tot'].values

    # Locate the minimum
    idx_min = np.argmin(Q_tot)
    ar_opt  = ar[idx_min]

    # Thermal aspect ratio label
    eta_label = '0.30' if label == 'rq03' else '0.70'
    title_str = (f"{cryogen_name.capitalize()}  –  "
                 f"η ≈ {eta_label}  (q_b / q_L = {q_b_mult:.4g})")

    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    ax.stackplot(ar, Q_L, Q_b, Q_VL, Q_Vw,
                 labels=['Wall → Liquid', 'Bottom',
                         'Vapour → Interface', 'Wall → Interface'],
                 colors=[c_L, c_b, c_VL, c_Vw])

    ax.plot(ar, Q_tot, color='black', linewidth=1.8, linestyle='-', zorder=5,
            label='Q$_{tot}$')

    # Mark minimum with a dashed vertical line
    ax.axvline(ar_opt, color='black', linewidth=0.9, linestyle='--', zorder=4)
    ar_span    = ar[-1] - ar[0]
    right_frac = (ar_opt - ar[0]) / ar_span
    if right_frac > 0.75:
        # Near right edge: place label to the left of the line
        ax.text(ar_opt - 0.06, ax.get_ylim()[1] * 0.97,
                f'AR* = {ar_opt:.2f}', fontsize=8, va='top', ha='right')
    else:
        ax.text(ar_opt + 0.06, ax.get_ylim()[1] * 0.97,
                f'AR* = {ar_opt:.2f}', fontsize=8, va='top', ha='left')

    ax.set_xlim(ar[0], ar[-1])
    ax.set_ylim(0)
    ax.set_xlabel('Aspect Ratio H/D  [–]', fontsize=11, fontweight='bold')
    ax.set_ylabel('Heat Ingress  /  W',    fontsize=11, fontweight='bold')
    ax.set_title(title_str, fontsize=10)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.85,
              edgecolor='gray', borderpad=0.6)
    ax.grid(True, alpha=0.25, linestyle='--', color='gray')

    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_svg, bbox_inches='tight')
    plt.close(fig)
    print(f"  Figures saved: {os.path.basename(out_png)}  {os.path.basename(out_svg)}")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    grand_start = time.time()

    for case in CASES:
        df = run_sweep(case)
        plot_case(df, case)

    grand_total = (time.time() - grand_start) / 3600
    print(f"\n{'='*70}")
    print(f"All cases complete in {grand_total:.2f} h.")
    print(f"Figures in:  {dir_fig}")
    print(f"Data in:     {dir_data}")
