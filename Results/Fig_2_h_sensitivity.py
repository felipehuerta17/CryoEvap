"""
Fig. 2 with sensitivity bands on the external convection coefficient.

Overlays the wall temperature profiles obtained for h_env at the baseline
value (14.839 W m^-2 K^-1) and at +/- 50 % of that baseline.
"""

import pickle as pkl

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 12

# Geometry (same as Fig_2.py)
e       = 10.6972 * 2 - 10.24 * 2
V_tank  = 88023.6952
a       = 0.5
d_i     = ((4 * V_tank) / (np.pi * a)) ** (1 / 3)
d_o     = d_i + e
r_i     = d_i / 2
r_o     = d_o / 2

# Sensitivity cases
cases = ['h_low', 'h_base', 'h_high']
case_labels = {
    'h_low':  r'$h_{env} \times 0.5$',
    'h_base': r'$h_{env}$ (baseline)',
    'h_high': r'$h_{env} \times 2.0$',
}
linestyles = {'h_low': ':', 'h_base': '-', 'h_high': '--'}

# Load data
transient = {}
long_term = {}
for c in cases:
    with open(f'Data/large_tank_data_LF95_transient_{c}.pkl', 'rb') as f:
        transient[c] = pkl.load(f)
    with open(f'Data/large_tank_data_LF95_3weeks_{c}.pkl', 'rb') as f:
        long_term[c] = pkl.load(f)

# Time selections (kept identical to the published figure)
time_tau = [1, 2, 4, 6, 8, 10]
horas_objetivo1 = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]

cmap = plt.get_cmap('inferno', len(time_tau) + 1)

fig = plt.figure(figsize=(15, 6), dpi=200)

# ---------------- Subplot 1: short-term profiles ----------------
ax1 = plt.subplot(1, 2, 1)

# Reference grid (same node count for the three cases by construction)
T_w_ref = transient['h_base']['T_w_raw']
r_grid_short = np.linspace(0, 1, T_w_ref.shape[0]) * (d_o - d_i) / 2 + d_i / 2

# Plot baseline as solid colored lines (the published curves)
for j, i in enumerate(time_tau):
    ax1.plot(r_grid_short, transient['h_base']['T_w_raw'][:, i],
             color=cmap(j), lw=1.6,
             label=fr't$^*$ = {i/10:.1f}')

# Overlay -50 % and +50 % envelopes in light colors
for c in ['h_low', 'h_high']:
    for j, i in enumerate(time_tau):
        ax1.plot(r_grid_short, transient[c]['T_w_raw'][:, i],
                 linestyle=linestyles[c], color=cmap(j), lw=1.0, alpha=0.85)

# Custom legend that combines time labels and h_env style
from matplotlib.lines import Line2D
time_handles = [Line2D([0], [0], color=cmap(j), lw=1.6,
                       label=fr't$^*$ = {i/10:.1f}')
                for j, i in enumerate(time_tau)]
style_handles = [
    Line2D([0], [0], color='k', lw=1.6, ls='-',  label=case_labels['h_base']),
    Line2D([0], [0], color='k', lw=1.0, ls=':',  label=case_labels['h_low']),
    Line2D([0], [0], color='k', lw=1.0, ls='--', label=case_labels['h_high']),
]
ax1.legend(handles=time_handles + style_handles, fontsize=fontsize_legend, ncol=2)
ax1.set_xlabel('Radius / m', fontsize=fontsize_label)
ax1.set_ylabel('Wall temperature / K', fontsize=fontsize_label)
ax1.tick_params(labelsize=fontsize_ticks)

# ---------------- Subplot 2: long-term profiles ----------------
ax2 = plt.subplot(1, 2, 2)

T_w_long_ref = long_term['h_base']['T_w_raw']
r_grid_long = np.linspace(0, 1, T_w_long_ref.shape[0]) * (r_o - r_i) + r_i
selected_index1 = [np.where(np.isclose(long_term['h_base']['Time'] / 3600, h))[0][0]
                   for h in horas_objetivo1]

for i in range(len(selected_index1)):
    ax2.plot(r_grid_long,
             pd.DataFrame(long_term['h_base']['T_w_raw']).iloc[:, selected_index1[i]],
             color=cmap(i), lw=1.6,
             label=f"t = {long_term['h_base']['Time'][selected_index1[i]]/3600:.0f} h")

for c in ['h_low', 'h_high']:
    for i, idx in enumerate(selected_index1):
        ax2.plot(r_grid_long,
                 pd.DataFrame(long_term[c]['T_w_raw']).iloc[:, idx],
                 linestyle=linestyles[c], color=cmap(i), lw=1.0, alpha=0.85)

ax2.legend(loc='lower right', fontsize=fontsize_legend)
ax2.set_xlabel('Radius / m', fontsize=fontsize_label)
ax2.set_xlim(r_i * 0.999, r_o * 1.001)
ax2.tick_params(labelsize=fontsize_ticks)

plt.tight_layout()
plt.savefig('Figures/Fig_2_h_sensitivity.svg', bbox_inches='tight', dpi=200)
print('Saved Figures/Fig_2_h_sensitivity.svg')
