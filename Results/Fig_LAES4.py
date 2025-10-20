import pandas as pd
import matplotlib.pyplot as plt

# Import data
data24h = pd.read_csv('Data/LAES_opti_168h_LFs.csv')
ratio = data24h['Geometric_AR']

# Optimal values of BOR
idx_min = data24h.iloc[:, 1:8].idxmin()
optimal_values = data24h.iloc[:, 1:8].min()

# Labels
labels = ["0.05", "0.20", "0.35", "0.50", "0.65", "0.80", "0.95"]

# Configuration
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth       = 2.5
paleta = plt.get_cmap('inferno', 9)
fig, axes = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

# Left subplot: BOR columns
ax = axes[0]
for i in range(7):
    ax.plot(ratio, data24h.iloc[:, 1 + i], label=fr'$LF ={labels[i]}$', color=paleta(i+1), linewidth=linewidth)
# Optimal points (one per LF)
ax.plot(ratio.iloc[idx_min.values], optimal_values.values, 'ko--', label='Optimal Values', markersize=5)
ax.set_xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
ax.set_ylabel('Boil-Off Rate (BOR) / %/day', fontsize=fontsize_label)
ax.tick_params(axis='both', labelsize=fontsize_ticks)
ax.set_ylim(0, 0.07)


# Right subplot: BOR columns
ax = axes[1]
for i in range(7):
    ax.plot(ratio, data24h.iloc[:, 8 + i], label=fr'$LF ={labels[i]}$', color=paleta(i+1), linewidth=linewidth)
ax.plot(ratio.iloc[idx_min.values], [data24h.iloc[idx_min.values[i], 8 + i] for i in range(7)], 'ko--', label='Optimal Values', markersize=5)
ax.set_xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
ax.set_ylabel('Thermal Aspect Ratio', fontsize=fontsize_label)
ax.tick_params(axis='both', labelsize=fontsize_ticks)
plt.legend(fontsize=fontsize_legend, loc=(-0.85, -0.3), ncol=4)

# plt.tight_layout()

plt.savefig("Figures/Fig_LAES4.svg", bbox_inches='tight', dpi=300)
# plt.show()