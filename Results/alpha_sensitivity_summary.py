"""
Numerical summary of the thermal-diffusivity sensitivity analysis.

Transient period definitions (consistent with the paper):
  tau_wc = L^2 / (16 * alpha)         -- wall-cooling time scale
                                           L = r_o - r_i (radial wall thickness)
  tau_r  = first local minimum of Q_L  -- onset of periodic (diurnal) variation
  Second transient period = tau_r - tau_wc

Run after alpha_sensitivity.py.
"""

import pickle

import numpy as np
from scipy.signal import argrelmin

# Geometry
e      = 10.6972 * 2 - 10.24 * 2
V_tank = 88023.6952
a      = 0.5
d_i    = ((4 * V_tank) / (np.pi * a)) ** (1 / 3)
d_o    = d_i + e
L      = (d_o - d_i) / 2   # radial wall thickness (r_o - r_i)

K_W_BASE = 0.0411
CP_W     = 900.0
RHO_W    = 60.0
DELTA    = 0.088

cases  = ['alpha_low', 'alpha_base', 'alpha_high']
k_vals = {
    'alpha_low':  K_W_BASE * (1 - DELTA),
    'alpha_base': K_W_BASE,
    'alpha_high': K_W_BASE * (1 + DELTA),
}
labels = {
    'alpha_low':  'alpha -8.8%',
    'alpha_base': 'baseline   ',
    'alpha_high': 'alpha +8.8%',
}

print('=== Thermal diffusivity (alpha_w = k_w / rho_w / cp_w) sensitivity ===\n')
print(f'  Wall thickness  L  = {L:.4f} m  (= r_o - r_i)')
print(f'  Baseline: k_w = {K_W_BASE:.4f} W m^-1 K^-1, '
      f'alpha_w = {K_W_BASE/(RHO_W*CP_W)*1e6:.4f} x10^-6 m^2 s^-1')
for c in ['alpha_low', 'alpha_high']:
    kw = k_vals[c]
    sign = '-' if c == 'alpha_low' else '+'
    print(f'  {sign}8.8%  : k_w = {kw:.5f} W m^-1 K^-1, '
          f'alpha_w = {kw/(RHO_W*CP_W)*1e6:.4f} x10^-6 m^2 s^-1')

print()
print(f"{'Case':<14} {'alpha [1e-6]':>12} {'tau_wc [h]':>11} "
      f"{'tau_r [h]':>10} {'2nd period [h]':>15} "
      f"{'Q_peak [kW]':>13} {'Q_ss [kW]':>11} "
      f"{'BOG_ss [kg/h]':>14} {'<T_w>@72h [K]':>15}")
print('-' * 120)

ref_q_ss   = None
ref_bog_ss = None
results    = {}

for c in cases:
    kw    = k_vals[c]
    alpha = kw / (RHO_W * CP_W)
    tau_wc_h = (L ** 2 / (16 * alpha)) / 3600.0

    with open(f'Data/large_tank_data_LF95_3weeks_alpha_{c}.pkl', 'rb') as fh:
        d = pickle.load(fh)

    t_h      = d['Time'] / 3600.0
    q_kw     = d['Q_tot'] / 1000.0
    ql_kw    = d['Q_L'] / 1000.0
    bog_kg_h = d['BOG'] * 3600.0

    # tau_r: first local minimum of Q_L (onset of periodic diurnal oscillation)
    idx_min = argrelmin(ql_kw, order=3)[0]
    tau_r_h = t_h[idx_min[0]] if len(idx_min) else float('nan')
    second_period = tau_r_h - tau_wc_h

    # Peak Q_tot during first 70 h (exclude t=0 initial-condition spike)
    mask_tr = (t_h > 0) & (t_h <= 70)
    q_peak  = q_kw[mask_tr].max()

    # Quasi-steady values at t = 504 h
    idx_504  = int(np.argmin(np.abs(t_h - 504)))
    q_ss     = q_kw[idx_504]
    bog_ss   = bog_kg_h[idx_504]

    # Average wall temperature at ~72 h
    idx_72 = int(np.argmin(np.abs(t_h - 72)))
    Tw_72  = d['Tw_avg'][idx_72]

    print(f"{labels[c]:<14} {alpha*1e6:>12.4f} {tau_wc_h:>11.2f} "
          f"{tau_r_h:>10.1f} {second_period:>15.1f} "
          f"{q_peak:>13.2f} {q_ss:>11.2f} "
          f"{bog_ss:>14.2f} {Tw_72:>15.2f}")

    results[c] = dict(q_ss=q_ss, bog_ss=bog_ss, q_peak=q_peak,
                      tau_wc_h=tau_wc_h, tau_r_h=tau_r_h,
                      second_period=second_period)
    if c == 'alpha_base':
        ref_q_ss   = q_ss
        ref_bog_ss = bog_ss

print()
print('Relative change vs baseline (quasi-steady, t = 504 h):')
for c in ['alpha_low', 'alpha_high']:
    dq   = (results[c]['q_ss']   - ref_q_ss)   / ref_q_ss   * 100
    dbog = (results[c]['bog_ss'] - ref_bog_ss)  / ref_bog_ss * 100
    print(f"  {labels[c]}:  Delta Q_tot = {dq:+.2f} %   Delta BOG = {dbog:+.2f} %")

print()
print('Transient-period comparison vs baseline:')
base = results['alpha_base']
for c in ['alpha_low', 'alpha_high']:
    r = results[c]
    print(f"  {labels[c]}:  "
          f"tau_wc = {r['tau_wc_h']:.2f} h (Delta = {r['tau_wc_h']-base['tau_wc_h']:+.2f} h),  "
          f"tau_r = {r['tau_r_h']:.1f} h (Delta = {r['tau_r_h']-base['tau_r_h']:+.1f} h),  "
          f"2nd period = {r['second_period']:.1f} h (Delta = {r['second_period']-base['second_period']:+.1f} h)")
print(f"  Baseline:      "
      f"tau_wc = {base['tau_wc_h']:.2f} h,  "
      f"tau_r = {base['tau_r_h']:.1f} h,  "
      f"2nd period = {base['second_period']:.1f} h")

# ------------------------------------------------------------------
# Phase-lag analysis: time delay between Q_tot and Q_Li peaks
# in the quasi-steady (periodic) regime (t = 200-504 h).
# Cross-correlation gives the lag at which Q_Li best matches Q_tot.
# Positive lag = Q_Li peaks AFTER Q_tot (wall diffusion delay).
# ------------------------------------------------------------------
from scipy.signal import argrelmax
import numpy as np

print()
print('Phase-lag between Q_tot and Q_Li (quasi-steady, t=200-504 h, cross-correlation):')
print(f"  Positive lag = Q_Li peak lags Q_tot (wall diffusion delays liquid heat ingress)\n")
print(f"  {'Case':<14} {'phase lag [h]':>15} {'phase lag [min]':>16}")
print('  ' + '-' * 47)

for c in cases:
    with open(f'Data/large_tank_data_LF95_3weeks_alpha_{c}.pkl', 'rb') as fh:
        d = pickle.load(fh)
    t_h      = d['Time'] / 3600.0
    q_tot_kw = d['Q_tot'] / 1000.0
    ql_kw    = d['Q_L']   / 1000.0

    mask_qs = (t_h >= 200) & (t_h <= 504)
    qt_qs   = q_tot_kw[mask_qs] - q_tot_kw[mask_qs].mean()
    ql_qs   = ql_kw[mask_qs]    - ql_kw[mask_qs].mean()
    dt_h    = t_h[1] - t_h[0]

    # Full cross-correlation; lag > 0 means ql_qs lags qt_qs
    xcorr   = np.correlate(ql_qs, qt_qs, mode='full')
    lags    = np.arange(-(len(qt_qs)-1), len(qt_qs)) * dt_h
    lag_h   = lags[np.argmax(xcorr)]
    lag_min = lag_h * 60

    print(f"  {labels[c]:<14} {lag_h:>15.2f} {lag_min:>16.1f}")
