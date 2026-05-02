"""
Print a numerical summary of the h_env sensitivity analysis: short-term
peak heat ingress, steady-state heat ingress, BOR and average wall
temperature for the baseline and the +/- 50 % cases.

Run after h_sensitivity.py.
"""

import pickle

import numpy as np

cases = ['h_low', 'h_base', 'h_high']
labels = {'h_low': 'h_env x 0.5', 'h_base': 'baseline', 'h_high': 'h_env x 2.0'}
h_values = {'h_low': 7.4195, 'h_base': 14.8390, 'h_high': 29.6780}

# Geometry / fluid (same as paper)
e = 10.6972 * 2 - 10.24 * 2
V_tank = 88023.6952
a = 0.5
d_i = ((4 * V_tank) / (np.pi * a)) ** (1 / 3)
LF = 0.95
rho_L = 600.7  # rough — value not used; BOR is computed from BOG / V_L_init * rho_L

print('=== External convection coefficient (h_env) sensitivity ===\n')
print(f"{'case':<12} {'h_env [W/m2K]':>14} {'Q_tot peak [kW]':>16} "
      f"{'Q_tot @504h [kW]':>18} {'BOG @504h [kg/h]':>18} "
      f"{'<T_w>(72 h) [K]':>17}")

for c in cases:
    with open(f'Data/large_tank_data_LF95_3weeks_{c}.pkl', 'rb') as f:
        d = pickle.load(f)

    t_h = d['Time'] / 3600
    q_tot_kw = d['Q_tot'] / 1000.0
    bog_kg_h = d['BOG'] * 3600.0

    # Peak Q_tot during the first 70 h (transient peak)
    mask_tr = t_h <= 70
    q_peak = q_tot_kw[mask_tr].max()

    # Quasi-stationary values at t = 504 h (week 3)
    idx_504 = int(np.argmin(np.abs(t_h - 504)))
    q_ss = q_tot_kw[idx_504]
    bog_ss = bog_kg_h[idx_504]

    # Average wall temperature at the end of the transient (~72 h)
    idx_72 = int(np.argmin(np.abs(t_h - 72)))
    Tw_72 = d['Tw_avg'][idx_72]

    print(f"{labels[c]:<12} {h_values[c]:>14.4f} {q_peak:>16.2f} "
          f"{q_ss:>18.2f} {bog_ss:>18.2f} {Tw_72:>17.2f}")

# Relative differences vs baseline
print('\nRelative change vs baseline (3-week / steady)')
with open('Data/large_tank_data_LF95_3weeks_h_base.pkl', 'rb') as f:
    base = pickle.load(f)
idx_504_base = int(np.argmin(np.abs(base['Time']/3600 - 504)))
ref_q = base['Q_tot'][idx_504_base]
ref_bog = base['BOG'][idx_504_base]
for c in ['h_low', 'h_high']:
    with open(f'Data/large_tank_data_LF95_3weeks_{c}.pkl', 'rb') as f:
        d = pickle.load(f)
    idx_504 = int(np.argmin(np.abs(d['Time']/3600 - 504)))
    dq = (d['Q_tot'][idx_504] - ref_q) / ref_q * 100
    dbog = (d['BOG'][idx_504] - ref_bog) / ref_bog * 100
    print(f"  {labels[c]:<14}  Delta Q_tot = {dq:+.2f}%   Delta BOG = {dbog:+.2f}%")
