"""
Sensitivity figures for thermal diffusivity alpha_w = k_w / (rho_w * cp_w).

Panel (a) — Fig-2 style: radial wall-temperature profiles at short-term
(transient, first 24 h) and long-term (3-week) periods.

Panel (b) / Panel (c) — Fig-3 style: heat-ingress decomposition (transient)
and BOG comparison (wall model vs simplified, with alpha sensitivity bands).

Run alpha_sensitivity.py first to generate the required .pkl files.
"""

import pickle as pkl

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.signal import argrelmin, argrelmax

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 12

# Geometry (same as Fig_2.py / Fig_3.py)
e      = 10.6972 * 2 - 10.24 * 2
V_tank = 88023.6952
a      = 0.5
d_i    = ((4 * V_tank) / (np.pi * a)) ** (1 / 3)
d_o    = d_i + e
r_i    = d_i / 2
r_o    = d_o / 2

cases = ['alpha_low', 'alpha_base', 'alpha_high']
K_W_BASE = 0.0411
DELTA    = 0.088
k_vals   = {
    'alpha_low':  K_W_BASE * (1 - DELTA),
    'alpha_base': K_W_BASE,
    'alpha_high': K_W_BASE * (1 + DELTA),
}
case_labels = {
    'alpha_low':  r'$\alpha_w\times(1-8.8\%)$',
    'alpha_base': r'$\alpha_w$ (baseline)',
    'alpha_high': r'$\alpha_w\times(1+8.8\%)$',
}
linestyles = {'alpha_low': ':', 'alpha_base': '-', 'alpha_high': '--'}

# Load data
transient = {}
long_term = {}
for c in cases:
    with open(f'Data/large_tank_data_LF95_transient_alpha_{c}.pkl', 'rb') as f:
        transient[c] = pkl.load(f)
    with open(f'Data/large_tank_data_LF95_3weeks_alpha_{c}.pkl', 'rb') as f:
        long_term[c] = pkl.load(f)

with open('Data/large_tank_data_LF95_3weeks_wo_wall.pkl', 'rb') as f:
    data_wo_wall = pkl.load(f)

# ================================================================
# Figure A: wall temperature profiles (Fig-2 style)
# ================================================================
time_tau       = [1, 2, 4, 6, 8, 10]
horas_objetivo = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
cmap = plt.get_cmap('inferno', len(time_tau) + 1)

fig_a, (ax_a1, ax_a2) = plt.subplots(1, 2, figsize=(15, 6), dpi=200)

# Short-term
T_w_ref     = transient['alpha_base']['T_w_raw']
r_grid_short = np.linspace(0, 1, T_w_ref.shape[0]) * (d_o - d_i) / 2 + d_i / 2

for j, i in enumerate(time_tau):
    ax_a1.plot(r_grid_short, transient['alpha_base']['T_w_raw'][:, i],
               color=cmap(j), lw=1.6, label=fr't$^*$ = {i/10:.1f}')

for c in ['alpha_low', 'alpha_high']:
    for j, i in enumerate(time_tau):
        ax_a1.plot(r_grid_short, transient[c]['T_w_raw'][:, i],
                   linestyle=linestyles[c], color=cmap(j), lw=1.0, alpha=0.85)

time_handles = [Line2D([0], [0], color=cmap(j), lw=1.6, label=fr't$^*$ = {i/10:.1f}')
                for j, i in enumerate(time_tau)]
style_handles = [
    Line2D([0], [0], color='k', lw=1.6, ls='-',  label=case_labels['alpha_base']),
    Line2D([0], [0], color='k', lw=1.0, ls=':',  label=case_labels['alpha_low']),
    Line2D([0], [0], color='k', lw=1.0, ls='--', label=case_labels['alpha_high']),
]
ax_a1.legend(handles=time_handles + style_handles, fontsize=fontsize_legend, ncol=2)
ax_a1.set_xlabel('Radius / m', fontsize=fontsize_label)
ax_a1.set_ylabel('Wall temperature / K', fontsize=fontsize_label)
ax_a1.tick_params(labelsize=fontsize_ticks)

# Long-term
T_w_long_ref = long_term['alpha_base']['T_w_raw']
r_grid_long  = np.linspace(0, 1, T_w_long_ref.shape[0]) * (r_o - r_i) + r_i
selected_idx = [np.where(np.isclose(long_term['alpha_base']['Time'] / 3600, h))[0][0]
                for h in horas_objetivo]

for i, idx in enumerate(selected_idx):
    ax_a2.plot(r_grid_long,
               pd.DataFrame(long_term['alpha_base']['T_w_raw']).iloc[:, idx],
               color=cmap(i), lw=1.6,
               label=f"t = {long_term['alpha_base']['Time'][idx]/3600:.0f} h")

for c in ['alpha_low', 'alpha_high']:
    for i, idx in enumerate(selected_idx):
        ax_a2.plot(r_grid_long,
                   pd.DataFrame(long_term[c]['T_w_raw']).iloc[:, idx],
                   linestyle=linestyles[c], color=cmap(i), lw=1.0, alpha=0.85)

ax_a2.legend(loc='lower right', fontsize=fontsize_legend)
ax_a2.set_xlabel('Radius / m', fontsize=fontsize_label)
ax_a2.set_xlim(r_i * 0.999, r_o * 1.001)
ax_a2.tick_params(labelsize=fontsize_ticks)

fig_a.tight_layout()
fig_a.savefig('Figures/Fig_2_alpha_sensitivity.svg', bbox_inches='tight', dpi=200)
print('Saved Figures/Fig_2_alpha_sensitivity.svg')

# ================================================================
# Figure B: heat ingress + BOG (Fig-3 style)
# ================================================================
cmap_q   = plt.get_cmap('inferno', 6)
cmap_bog = plt.get_cmap('inferno', 4)
time_tr     = 70   # h — heat ingress panel
time_tr_bog = 168  # h — BOG panel

fig_b, (ax_b1, ax_b2) = plt.subplots(1, 2, figsize=(15, 6), dpi=200)

# ---------- Compute tau_wc and tau_r per case ----------
L = (d_o - d_i) / 2   # radial wall thickness
tau_info = {}
for c in cases:
    kw     = k_vals[c]
    alpha  = kw / (60.0 * 900.0)
    tau_wc = L**2 / (16 * alpha) / 3600.0   # h
    ql_all = long_term[c]['Q_L'] / 1000.0
    t_all  = long_term[c]['Time'] / 3600.0
    idx_m  = argrelmin(ql_all, order=3)[0]
    tau_r  = t_all[idx_m[0]] if len(idx_m) else float('nan')
    tau_info[c] = dict(tau_wc=tau_wc, tau_r=tau_r)

# ---------- Panel (a): heat-ingress decomposition ----------
base    = long_term['alpha_base']
mask_tr = base['Time'] <= time_tr * 3600

ax_b1.plot(base['Time'][mask_tr]/3600, base['Q_tot'][mask_tr]/1000,
           label=r'$\dot{Q}_{\mathrm{tot}}$', color=cmap_q(0))
ax_b1.plot(base['Time'][mask_tr]/3600, base['Q_L'][mask_tr]/1000,
           label=r'$\dot{Q}_{\mathrm{Li}}$', color=cmap_q(4))
ax_b1.plot(base['Time'][mask_tr]/3600, base['Q_b'][mask_tr]/1000,
           label=r'$\dot{Q}_{\mathrm{b}}$', color=cmap_q(3))
ax_b1.plot(base['Time'][mask_tr]/3600, base['Q_Vw'][mask_tr]/1000,
           label=r'$\dot{Q}_{\mathrm{Wi}}$', color=cmap_q(2))

for c in ['alpha_low', 'alpha_high']:
    d = long_term[c]
    m = d['Time'] <= time_tr * 3600
    ax_b1.plot(d['Time'][m]/3600, d['Q_tot'][m]/1000,
               linestyles[c], color=cmap_q(0), lw=1.0, alpha=0.85)
    ax_b1.plot(d['Time'][m]/3600, d['Q_L'][m]/1000,
               linestyles[c], color=cmap_q(4), lw=1.0, alpha=0.85)

# Vertical markers for tau_wc and tau_r (baseline only, annotated)
vline_colors = {'alpha_low': 'steelblue', 'alpha_base': 'dimgrey', 'alpha_high': 'firebrick'}
for c in cases:
    twc = tau_info[c]['tau_wc']
    tr  = tau_info[c]['tau_r']
    col = vline_colors[c]
    ax_b1.axvline(twc, color=col, lw=0.9, ls=linestyles[c], alpha=0.7)
    ax_b1.axvline(tr,  color=col, lw=0.9, ls=linestyles[c], alpha=0.7)

# Annotate baseline markers only
twc_b = tau_info['alpha_base']['tau_wc']
tr_b  = tau_info['alpha_base']['tau_r']
ax_b1.text(twc_b + 0.5, 90, r'$\tau_{wc}$', fontsize=fontsize_legend, color='dimgrey')
ax_b1.text(tr_b  + 0.5, 90, r'$\tau_{r}$',  fontsize=fontsize_legend, color='dimgrey')

ax_b1.set_xlim(1, time_tr)
ax_b1.set_ylim(0, 100)
ax_b1.set_ylabel(r'Heat transfer rate / kW', fontsize=fontsize_label)
ax_b1.set_xlabel('Time / h', fontsize=fontsize_label)
ax_b1.tick_params(labelsize=fontsize_ticks)

q_handles = [
    Line2D([0], [0], color=cmap_q(0), label=r'$\dot{Q}_{\mathrm{tot}}$'),
    Line2D([0], [0], color=cmap_q(4), label=r'$\dot{Q}_{\mathrm{Li}}$'),
    Line2D([0], [0], color=cmap_q(3), label=r'$\dot{Q}_{\mathrm{b}}$'),
    Line2D([0], [0], color=cmap_q(2), label=r'$\dot{Q}_{\mathrm{Wi}}$'),
]
style_handles_b = [
    Line2D([0], [0], color='k', ls='-',  label=case_labels['alpha_base']),
    Line2D([0], [0], color='k', ls=':',  label=case_labels['alpha_low']),
    Line2D([0], [0], color='k', ls='--', label=case_labels['alpha_high']),
]
ax_b1.legend(handles=q_handles + style_handles_b, fontsize=fontsize_legend,
             loc='upper right', ncol=2)

# ---------- Panel (b): BOG rate ----------
m_bog = base['Time'] < time_tr_bog * 3600
ax_b2.plot(base['Time'][m_bog]/3600, base['BOG'][m_bog] * 3600,
           label=case_labels['alpha_base'], color=cmap_bog(2), lw=1.6)

for c in ['alpha_low', 'alpha_high']:
    d = long_term[c]
    m = d['Time'] < time_tr_bog * 3600
    ax_b2.plot(d['Time'][m]/3600, d['BOG'][m] * 3600,
               linestyles[c], color=cmap_bog(2), lw=1.0, alpha=0.85,
               label=case_labels[c])

m_wo = data_wo_wall['Time'] < time_tr_bog * 3600
ax_b2.plot(data_wo_wall['Time'][m_wo]/3600,
           data_wo_wall['BOG'][m_wo] * 3600,
           label='Simplified (no wall)', color=cmap_bog(1), lw=1.6)

ax_b2.set_xlabel('Time / h', fontsize=fontsize_label)
ax_b2.set_xlim(1, 70)
ax_b2.set_ylim(50, 200)
ax_b2.set_ylabel(r'Boil-off gas rate / kg h$^{-1}$', fontsize=fontsize_label)
ax_b2.tick_params(labelsize=fontsize_ticks)
ax_b2.legend(loc='upper right', fontsize=fontsize_legend)

fig_b.tight_layout()
fig_b.savefig('Figures/Fig_3_alpha_sensitivity.svg', dpi=300, bbox_inches='tight')
print('Saved Figures/Fig_3_alpha_sensitivity.svg')

# ================================================================
# Figure C: BOG only — standalone three-line comparison
# ================================================================
fig_c, ax_c = plt.subplots(figsize=(8, 6), dpi=200)

colors = {
    'alpha_low':  cmap_bog(1),
    'alpha_base': cmap_bog(2),
    'alpha_high': cmap_bog(3),
}

for c in cases:
    d = long_term[c]
    m = d['Time'] < time_tr_bog * 3600
    ax_c.plot(d['Time'][m]/3600, d['BOG'][m] * 3600,
              linestyles[c], color=colors[c], lw=1.6, label=case_labels[c])

m_wo = data_wo_wall['Time'] < time_tr_bog * 3600
ax_c.plot(data_wo_wall['Time'][m_wo]/3600,
          data_wo_wall['BOG'][m_wo] * 3600,
          '-.', color='grey', lw=1.4, label='Simplified (no wall)')

ax_c.set_xlabel('Time / h', fontsize=fontsize_label)
ax_c.set_xlim(1, 70)
ax_c.set_ylim(50, 200)
ax_c.set_ylabel(r'Boil-off gas rate / kg h$^{-1}$', fontsize=fontsize_label)
ax_c.tick_params(labelsize=fontsize_ticks)
ax_c.legend(loc='upper right', fontsize=fontsize_legend)

fig_c.tight_layout()
fig_c.savefig('Figures/Fig_BOG_alpha_sensitivity.svg', dpi=300, bbox_inches='tight')
print('Saved Figures/Fig_BOG_alpha_sensitivity.svg')
