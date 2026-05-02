"""
Fig. 3 with sensitivity bands on the external convection coefficient.

Reproduces the two panels of Fig. 3 (heat ingress decomposition during the
transient and BOG comparison) overlaying h_env x 0.5 and x 1.5 against the
baseline run. The simplified (no-wall) curve in panel (b) is the original
case from the paper because the simplified model does not depend on h_env
through the wall ODE — only through U_L, which is unchanged by < 0.5 % across
the swept range.
"""

import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 12

# Sensitivity datasets (wall model)
cases = ['h_low', 'h_base', 'h_high']
ls = {'h_low': ':', 'h_base': '-', 'h_high': '--'}
case_labels = {
    'h_low':  r'$h_{env}\times0.5$',
    'h_base': r'$h_{env}$ baseline',
    'h_high': r'$h_{env}\times2.0$',
}
data = {}
for c in cases:
    with open(f'Data/large_tank_data_LF95_3weeks_{c}.pkl', 'rb') as f:
        data[c] = pickle.load(f)

# Simplified (wo-wall) baseline retained from the original paper
with open('Data/large_tank_data_LF95_3weeks_wo_wall.pkl', 'rb') as f:
    data_wo_wall = pickle.load(f)

# ---------------- Panel (a): heat ingresses during transient ---------------
cmap_q = plt.get_cmap('inferno', 6)
time_tr = 70  # h
mask_tr = data['h_base']['Time'] <= time_tr * 3600

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=200)

base = data['h_base']
ax1.plot(base['Time'][mask_tr]/3600, base['Q_tot'][mask_tr]/1000,
         label=r'$\dot{Q}_{\mathrm{tot}}$', color=cmap_q(0))
ax1.plot(base['Time'][mask_tr]/3600, base['Q_L'][mask_tr]/1000,
         label=r'$\dot{Q}_{\mathrm{Li}}$', color=cmap_q(4))
ax1.plot(base['Time'][mask_tr]/3600, base['Q_b'][mask_tr]/1000,
         label=r'$\dot{Q}_{\mathrm{b}}$', color=cmap_q(3))
ax1.plot(base['Time'][mask_tr]/3600, base['Q_Vw'][mask_tr]/1000,
         label=r'$\dot{Q}_{\mathrm{Wi}}$', color=cmap_q(2))

# Overlay sensitivity bands for Q_tot and Q_L (the two strongly affected curves)
for c in ['h_low', 'h_high']:
    d = data[c]
    m = d['Time'] <= time_tr * 3600
    ax1.plot(d['Time'][m]/3600, d['Q_tot'][m]/1000, ls[c], color=cmap_q(0),
             lw=1.0, alpha=0.85)
    ax1.plot(d['Time'][m]/3600, d['Q_L'][m]/1000, ls[c], color=cmap_q(4),
             lw=1.0, alpha=0.85)

ax1.set_xlim(1, time_tr)
ax1.set_ylim(0, 100)
ax1.set_ylabel(r'Heat transfer rate / kW', fontsize=fontsize_label)
ax1.set_xlabel('Time / h', fontsize=fontsize_label)
ax1.tick_params(labelsize=fontsize_ticks)

q_handles = [
    Line2D([0], [0], color=cmap_q(0), label=r'$\dot{Q}_{\mathrm{tot}}$'),
    Line2D([0], [0], color=cmap_q(4), label=r'$\dot{Q}_{\mathrm{Li}}$'),
    Line2D([0], [0], color=cmap_q(3), label=r'$\dot{Q}_{\mathrm{b}}$'),
    Line2D([0], [0], color=cmap_q(2), label=r'$\dot{Q}_{\mathrm{Wi}}$'),
]
style_handles = [
    Line2D([0], [0], color='k', ls='-',  label=case_labels['h_base']),
    Line2D([0], [0], color='k', ls=':',  label=case_labels['h_low']),
    Line2D([0], [0], color='k', ls='--', label=case_labels['h_high']),
]
ax1.legend(handles=q_handles + style_handles, fontsize=fontsize_legend,
           loc='upper right', ncol=2)

# ---------------- Panel (b): BOG transient ----------------
cmap2 = plt.get_cmap('inferno', 4)
time_tr_bog = 168

# Baseline wall-model and simplified curves (kept as in Fig_3.py)
m_bog = base['Time'] < time_tr_bog * 3600
ax2.plot(base['Time'][m_bog]/3600, base['BOG'][m_bog] * 3600,
         label='Wall model', color=cmap2(2), lw=1.6)

m_wo = data_wo_wall['Time'] < time_tr_bog * 3600
ax2.plot(data_wo_wall['Time'][m_wo]/3600,
         data_wo_wall['BOG'][m_wo] * 3600,
         label='Simplified', color=cmap2(1), lw=1.6)

# Overlay h_env sensitivity on the wall-model BOG
for c in ['h_low', 'h_high']:
    d = data[c]
    m = d['Time'] < time_tr_bog * 3600
    ax2.plot(d['Time'][m]/3600, d['BOG'][m] * 3600, ls[c], color=cmap2(2),
             lw=1.0, alpha=0.85,
             label=f"Wall, {case_labels[c]}")

ax2.set_xlabel('Time / h', fontsize=fontsize_label)
ax2.set_xlim(1, 70)
ax2.set_ylim(50, 200)
ax2.set_ylabel(r'Boil-off gas rate / kg h$^{-1}$', fontsize=fontsize_label)
ax2.tick_params(labelsize=fontsize_ticks)
ax2.legend(loc='upper right', fontsize=fontsize_legend, ncol=2)

plt.tight_layout()
plt.savefig('Figures/Fig_3_h_sensitivity.svg', dpi=300, bbox_inches='tight')
print('Saved Figures/Fig_3_h_sensitivity.svg')
