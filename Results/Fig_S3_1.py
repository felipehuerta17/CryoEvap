import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv('Data/simulation_times_big_o.csv')

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

n = data['r_grid'] + data['z_grid']
t = data['mean_time']

p, log10_C = np.polyfit(np.log10(n), np.log10(t), 1)
C = 10**log10_C

n_fit = np.geomspace(n.min(), n.max(), 200)
t_fit = C * n_fit**p

cmap = plt.get_cmap('inferno')
color_data = cmap(0.25)
color_fit = cmap(0.80)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(n, t, 'o-', color=color_data, label='Empirical data')
ax.plot(n_fit, t_fit, '--', color=color_fit, label=fr'Fit: $O(N^{{{p:.2f}}})$')

ax.set_xscale('log', base=2)
ax.set_yscale('log')

k_min = int(np.floor(np.log2(n.min())))
k_max = int(np.ceil(np.log2(n.max())))
x_ticks = 2 ** np.arange(k_min, k_max + 1)
ax.set_xticks(x_ticks)
ax.set_xticklabels([fr'$2^{{{k}}}$' for k in range(k_min, k_max + 1)])

ax.set_xlabel(r'Total number of nodes ($N = r_{grid} + z_{grid}$)', fontsize=fontsize_label)
ax.set_ylabel('Simulation time (s)', fontsize=fontsize_label)

ax.tick_params(axis='both', which='both', labelsize=fontsize_ticks)
ax.legend(fontsize=fontsize_legend)
ax.grid(True, which='both', ls=':', alpha=0.7)

plt.tight_layout()
plt.savefig('Figures/Fig_S3_1.svg', format='svg')
