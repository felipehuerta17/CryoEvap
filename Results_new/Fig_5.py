import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# ---------------------------------------------------------
# 1. PARAMETERS & CONFIGURATION
# ---------------------------------------------------------
params = [r'$U_L$', r'$U_V$', r'$U_b$', r'$\eta_w$'] # Includes all parameters
x      = np.arange(len(params))
width  = 0.35

paleta = plt.get_cmap('inferno', 4)
fontsize_label  = 14
fontsize_ticks  = 16
fontsize_legend = 14
linewidth       = 2.5

# ---------------------------------------------------------
# 2. DATA IMPORT & AVERAGING
# ---------------------------------------------------------
folder = "Data/"

# Read CSVs using 'Param' as index to allow direct dataframe arithmetic
LF05_ru025 = pd.read_csv(os.path.join(folder, 'Sensitivity_LF05_ru_025.csv'), index_col='Param')
LF05_ru4   = pd.read_csv(os.path.join(folder, 'Sensitivity_LF05_ru_4.csv'), index_col='Param')

LF95_ru025 = pd.read_csv(os.path.join(folder, 'Sensitivity_LF95_ru_025.csv'), index_col='Param')
LF95_ru4   = pd.read_csv(os.path.join(folder, 'Sensitivity_LF95_ru_4.csv'), index_col='Param')

# Average the dataframes
LF05_avg = (LF05_ru025 + LF05_ru4) / 2.0
LF95_avg = (LF95_ru025 + LF95_ru4) / 2.0

# Extract values for LF 0.95
bor_sens_95 = LF95_avg['Local_sensitivity'].to_numpy()
des_sens_95 = LF95_avg['Design_elasticity'].to_numpy()

# Extract values for LF 0.05
bor_sens_05 = LF05_avg['Local_sensitivity'].to_numpy()
des_sens_05 = LF05_avg['Design_elasticity'].to_numpy()

# ---------------------------------------------------------
# 3. PLOTTING
# ---------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

# Left Subplot: Local Sensitivity
rects1 = ax1.bar(x + width/2, bor_sens_95, width, label=r'LF$_0$=0.95', color=paleta(2), alpha=0.9)
rects2 = ax1.bar(x - width/2, bor_sens_05, width, label=r'LF$_0$=0.05', color=paleta(1), alpha=0.9)

ax1.set_ylabel('Local Sensitivity ($S_p$)', fontsize=fontsize_label)
ax1.text(0.03, 0.97, 'a)', transform=ax1.transAxes, fontsize=18, fontweight='bold', va='top')
ax1.set_xticks(x)
ax1.set_xticklabels(params, fontsize=fontsize_ticks)
ax1.tick_params(axis='both', labelsize=fontsize_ticks)
ax1.legend(prop={'size': fontsize_legend}, loc='upper right')
ax1.grid(axis='y', linestyle='--', alpha=0.3)
ax1.set_ylim(0, 0.8)

# Right Subplot: Design Elasticity
rects3 = ax2.bar(x + width/2, des_sens_95, width, label=r'LF$_0$=0.95', color=paleta(2), alpha=0.9)
rects4 = ax2.bar(x - width/2, des_sens_05, width, label=r'LF$_0$=0.05', color=paleta(1), alpha=0.9)

ax2.set_ylabel(r'Design Elasticity ($\epsilon_{a^*,p}$)', fontsize=fontsize_label)
ax2.text(0.03, 0.97, 'b)', transform=ax2.transAxes, fontsize=18, fontweight='bold', va='top')
ax2.set_xticks(x)
ax2.set_xticklabels(params, fontsize=fontsize_ticks)
ax2.tick_params(axis='both', labelsize=fontsize_ticks)
ax2.legend(prop={'size': fontsize_legend})
ax2.grid(axis='y', linestyle='--', alpha=0.3)
ax2.axhline(0, color='black', linewidth=0.8)
# ax2.set_ylim(-1.2, 0.0)

# ---------------------------------------------------------
# 4. EXPORT FIGURE
# ---------------------------------------------------------
plt.tight_layout()
plt.savefig("Figures/Fig_5.svg", bbox_inches='tight', dpi=300)
plt.close()

# Print values for reporting
print("Local Sensitivity (LF=0.95):", bor_sens_95)
print("Local Sensitivity (LF=0.05):", bor_sens_05)
print("Design Elasticity (LF=0.95):", des_sens_95)
print("Design Elasticity (LF=0.05):", des_sens_05)
