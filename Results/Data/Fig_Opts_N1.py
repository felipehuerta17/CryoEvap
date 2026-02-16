import pandas as pd
import matplotlib.pyplot as plt

optimos = pd.read_csv('all_opts.csv')

# Labels
label = ["LAES", "LH2", "LNG"]
n = len(label)
# Configuration
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth       = 2.5
paleta = plt.get_cmap('inferno', n+2)
fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

# Left subplot: BOR columns rq03
ax = axes[0]
for i in range(n):
    data = optimos[optimos["Cryogen"] == label[i]+"_rq03"]
    ax.plot(data['Optimal_AR'], data['LF'], label=label[i], linewidth=linewidth, color=paleta(i+1))
# Optimal points (one per LF)
ax.set_xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
ax.set_ylabel('Liquid Filling', fontsize=fontsize_label)
ax.tick_params(axis='both', labelsize=fontsize_ticks)
ax.set_ylim(0.05, 0.95)
ax.set_xlim([0.4, 0.5])
ax.text(0.03, 0.97, 'a)', transform=ax.transAxes, fontsize=18, fontweight='bold', va='top', font='Arial')


# Right subplot: BOR columns rq07
ax = axes[1]
for i in range(n):
    data = optimos[optimos["Cryogen"] == label[i]+"_rq07"]
    ax.plot(data['Optimal_AR'], data['LF'], label=label[i], linewidth=linewidth, color=paleta(i+1))
ax.set_xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
# ax.set_ylabel('Thermal Aspect Ratio', fontsize=fontsize_label)
ax.tick_params(axis='both', labelsize=fontsize_ticks)
ax.set_ylim(0.05, 0.95)
ax.set_xlim([2.2, 2.8])
plt.legend(fontsize=fontsize_legend, loc=(-0.5, -0.25), ncol=4)
ax.text(0.03, 0.97, 'b)', transform=ax.transAxes, fontsize=18, fontweight='bold', va='top', font='Arial')

plt.savefig("Figures/Fig_opts_N1.svg", bbox_inches='tight', dpi=300)