import pandas as pd
import matplotlib.pyplot as plt

folder = "Data/"
data_30d_rq03 = pd.read_csv(folder+'LAES_opti_30d_LFs_rq03.csv')
optimos_30d_rq03 = pd.read_csv(folder+'LAES_opti_30d_LFs_rq03_opts.csv')
data_30d_rq07 = pd.read_csv(folder+'LAES_opti_30d_LFs_rq07.csv')
optimos_30d_rq07 = pd.read_csv(folder+'LAES_opti_30d_LFs_rq07_opts.csv')


# Labels
label = ["LF 0.05", "LF 0.50", "LF 0.95"]
labels = ["BOR_LF_0.05", "BOR_LF_0.50", "BOR_LF_0.95"]
n = len(labels)
# Configuration
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth       = 2.5
paleta = plt.get_cmap('inferno', n+2)
fig, axes = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

# Left subplot: BOR columns
ax = axes[0]
for i in range(n):
    ax.plot(data_30d_rq03['Geometric_AR'], data_30d_rq03[labels[i]], label=label[i], linewidth=linewidth, color=paleta(i+1))
# Optimal points (one per LF)
ax.plot(optimos_30d_rq03['Optimal_AR'], optimos_30d_rq03['Min_BOR'], label='Optimal Values', linewidth=linewidth, color=paleta(0))
ax.set_xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
ax.set_ylabel('Boil-Off Rate (BOR) / %/day', fontsize=fontsize_label)
ax.tick_params(axis='both', labelsize=fontsize_ticks)
ax.set_ylim(0, 0.07)
ax.set_xlim([0.2, 3])
ax.text(0.03, 0.97, 'a)', transform=ax.transAxes, fontsize=18, fontweight='bold', va='top', font='Arial')


# Right subplot: BOR columns
ax = axes[1]
for i in range(n):
    ax.plot(data_30d_rq07['Geometric_AR'], data_30d_rq07[labels[i]], label=label[i], linewidth=linewidth, color=paleta(i+1))
ax.plot(optimos_30d_rq07['Optimal_AR'], optimos_30d_rq07['Min_BOR'], label='Optimal Values', linewidth=linewidth, color=paleta(0))
ax.set_xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
# ax.set_ylabel('Thermal Aspect Ratio', fontsize=fontsize_label)
ax.tick_params(axis='both', labelsize=fontsize_ticks)
ax.set_ylim(0, 0.2)
ax.set_xlim([0.2, 3])
plt.legend(fontsize=fontsize_legend, loc=(-0.85, -0.25), ncol=4)
ax.text(0.03, 0.97, 'b)', transform=ax.transAxes, fontsize=18, fontweight='bold', va='top', font='Arial')

plt.savefig("Figures/Fig_LAES_N1.svg", bbox_inches='tight', dpi=300)