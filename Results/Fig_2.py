import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

fontsize_label = 14
fontsize_ticks = 14
fontsize_legend = 14

with open('Data/large_tank_LF95_transient.pkl', 'rb') as f:
    large_tank = pkl.load(f)

d_i    = 10.24*2     # Internal diameter / m
d_o    = 10.6972*2   # External diameter / m

r = np.linspace(0, 1, large_tank['T_w_raw'].shape[0]) * (d_o - d_i)/2 + d_i/2
time_tau = [1, 2, 4, 6, 8, 10]

cmap = plt.get_cmap('inferno', len(time_tau)+1)

plt.figure(figsize= (8,6), dpi=200)

for j, i in enumerate(time_tau):
    plt.plot(r, large_tank['T_w_raw'][:, i], label=fr't$^*$ = {i/10:.1f}', color = cmap(j))

plt.legend(fontsize = fontsize_legend)
plt.xlabel('Radius / m', fontsize=fontsize_label)
plt.ylabel('Temperature / K', fontsize=fontsize_label)
plt.xticks(fontsize=fontsize_ticks)
plt.yticks(fontsize=fontsize_ticks)

plt.savefig('Figures/Fig_2.svg', bbox_inches='tight', dpi=200)
