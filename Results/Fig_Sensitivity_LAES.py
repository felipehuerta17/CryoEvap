import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Params
params = [r'$U_L$', r'$U_V$', r'$\eta_w$', r'$\eta$']
x      = np.arange(len(params))
width  = 0.35

# Configuration
paleta = plt.get_cmap('inferno', 4)
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth       = 2.5

# Import data
LF95_data = pd.read_csv('Data/Parameter_sensitivity_LAES_95.csv')
LF05_data = pd.read_csv('Data/Parameter_sensitivity_LAES_05.csv')

# Datos LF 0.95
bor_sens_95 = LF95_data.Local_sensitivity.to_numpy()
des_sens_95 = LF95_data.Design_elasticity.to_numpy()

# Datos LF 0.05 
bor_sens_05 = LF05_data.Local_sensitivity.to_numpy()
des_sens_05 = LF05_data.Design_elasticity.to_numpy()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

rects1 = ax1.bar(x + width/2, bor_sens_95, width, label='LF=0.95', color=paleta(2), alpha=0.9)
rects2 = ax1.bar(x - width/2, bor_sens_05, width, label='LF=0.05', color=paleta(1), alpha=0.9)

ax1.set_ylabel('Local Sensitivity ($S_p$)', fontsize=fontsize_label)
ax1.text(0.03, 0.97, 'a)', transform=ax1.transAxes, fontsize=18, fontweight='bold', va='top')
ax1.set_xticks(x)
ax1.set_xticklabels(params, fontsize=fontsize_ticks)
ax1.tick_params(axis='both', labelsize=fontsize_ticks)
ax1.legend(prop={'size': fontsize_legend}, loc='lower left')
ax1.grid(axis='y', linestyle='--', alpha=0.3)
# ax1.set_ylim(0, 0.8)

rects3 = ax2.bar(x + width/2, des_sens_95, width, label='LF=0.95', color=paleta(2), alpha=0.9)
rects4 = ax2.bar(x - width/2, des_sens_05, width, label='LF=0.05', color=paleta(1), alpha=0.9)

ax2.set_ylabel(r'Design Elasticity ($\epsilon_{a^*,p}$)', fontsize=fontsize_label)
ax2.text(0.03, 0.97, 'b)', transform=ax2.transAxes, fontsize=18, fontweight='bold', va='top')
ax2.set_xticks(x)
ax2.set_xticklabels(params, fontsize=fontsize_ticks)
ax2.tick_params(axis='both', labelsize=fontsize_ticks)
ax2.legend(prop={'size': fontsize_legend})
ax2.grid(axis='y', linestyle='--', alpha=0.3)
ax2.axhline(0, color='black', linewidth=0.8)

plt.tight_layout()
plt.savefig("Figures/Fig_LAES_Sensitivity.svg", bbox_inches='tight', dpi=300)

