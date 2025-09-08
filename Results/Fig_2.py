
import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

with open('Data/large_tank_data_LF95_transient.pkl', 'rb') as f:
    large_tank = pkl.load(f)

with open('Data/large_tank_data_LF95_3weeks.pkl', 'rb') as f:
    data_LF95 = pkl.load(f)

# Dimensions
e       = 10.6972*2 - 10.24*2 # Wall thickness m
V_tank  = 88023.6952 #m^3
a       = 0.5 # Aspect ratio
d_i     = ((4 * V_tank)/(np.pi * a))**(1/3) # internal diameter / m
d_o     = d_i + e # external diameter / m
r_i     = d_i/2
r_o     = d_o/2
A_T     = np.pi * (0.5*d_i)**2 # m^2

r = np.linspace(0, 1, large_tank['T_w_raw'].shape[0]) * (d_o - d_i)/2 + d_i/2
z_adim = np.linspace(0, 1, large_tank['T_V_raw'].shape[0])
time_tau = [1, 2, 4, 6, 8, 10]

cmap = plt.get_cmap('inferno', len(time_tau)+1)

plt.figure(figsize= (15,6), dpi=200)

plt.subplot(1,2,1)
for j, i in enumerate(time_tau):
    plt.plot(r, large_tank['T_w_raw'][:, i], label=fr't$^*$ = {i/10:.1f}', color = cmap(j))

plt.legend(fontsize = fontsize_legend)
plt.xlabel('Radius / m', fontsize=fontsize_label)
plt.ylabel('Wall temperature / K', fontsize=fontsize_label)
plt.xticks(fontsize=fontsize_ticks)
plt.yticks(fontsize=fontsize_ticks)


# First subplot - Early hours temperature profiles
T_w_csv = pd.DataFrame(data_LF95['T_w_raw'])
r_adim = np.linspace(0, 1, T_w_csv.shape[0])
r_grid = r_adim * (r_o - r_i) + r_i

horas_objetivo1 = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
selected_index1 = [np.where(np.isclose(data_LF95['Time']/3600, h))[0][0] for h in horas_objetivo1]

plt.subplot(1,2,2)
for i in range(len(selected_index1)):
    plt.plot(r_grid, T_w_csv.iloc[:, selected_index1[i]], color=cmap(i), 
             label=f't = {data_LF95["Time"][selected_index1[i]]/3600:.0f} h')
plt.legend(loc='lower right', fontsize=fontsize_legend)
# plt.ylabel('Temperature / K', fontsize=fontsize_label)
plt.xlabel('Radius / m', fontsize=fontsize_label)
plt.xlim(r_i*0.999, r_o*1.001)
plt.tick_params(labelsize=fontsize_ticks)


plt.savefig('Figures/Fig_2.svg', bbox_inches='tight', dpi=200)