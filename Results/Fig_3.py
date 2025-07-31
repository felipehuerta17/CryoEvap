import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

with open('Data/large_tank_data_LF95_3weeks.pkl', 'rb') as f:
    data_LF95 = pickle.load(f)

# Create figure with subplots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(17, 5), dpi=300)

# First subplot - Early hours temperature profiles
T_w_csv = pd.DataFrame(data_LF95['T_w_raw'])
r_i = 10.24    
r_o = 10.6972  
r_adim = np.linspace(0, 1, T_w_csv.shape[0])
r_grid = r_adim * (r_o - r_i) + r_i

horas_objetivo1 = [5.0, 10.0, 15.0, 20.0, 25.0]
selected_index1 = [np.where(np.isclose(data_LF95['Time']/3600, h))[0][0] for h in horas_objetivo1]
cmap1 = plt.get_cmap('inferno', len(selected_index1)+1)

for i in range(len(selected_index1)):
    ax1.plot(r_grid, T_w_csv.iloc[:, selected_index1[i]], color=cmap1(i), 
             label=f't = {data_LF95["Time"][selected_index1[i]]/3600:.0f} h')
ax1.legend(loc='lower right', fontsize=fontsize_legend)
ax1.set_ylabel('Temperature / K', fontsize=fontsize_label)
ax1.set_xlabel('Radius / m', fontsize=fontsize_label)
ax1.set_xlim(r_i*0.999, r_o*1.001)
ax1.tick_params(labelsize=fontsize_ticks)

# Second subplot - Later hours temperature profiles
horas_objetivo2 = [50.0, 55.0, 60.0, 65.0, 70.0]
selected_index2 = [np.where(np.isclose(data_LF95['Time']/3600, h))[0][0] for h in horas_objetivo2]
cmap2 = plt.get_cmap('inferno', len(selected_index2)+1)

for i in range(len(selected_index2)):
    ax2.plot(r_grid, T_w_csv.iloc[:, selected_index2[i]], color=cmap2(i), 
             label=f't = {data_LF95["Time"][selected_index2[i]]/3600:.0f} h')
ax2.legend(loc='lower right', fontsize=fontsize_legend)
ax2.set_ylabel('Temperature / K', fontsize=fontsize_label)
ax2.set_xlabel('Radius / m', fontsize=fontsize_label)
ax2.set_xlim(r_i*0.999, r_o*1.001)
ax2.tick_params(labelsize=fontsize_ticks)

# Third subplot - Heat transfer rates
cmap3 = plt.get_cmap('inferno', 4)
ax3.plot(data_LF95['Time']/3600, -data_LF95['Q_w_L']/1000, label=r'$\dot{Q}_{\text{w,L}}$', color=cmap3(1))
ax3.plot(data_LF95['Time']/3600, -data_LF95['Q_env_w']/1000, label=r'$\dot{Q}_{\text{env,w}}$', color=cmap3(2))
ax3.set_xlabel('Time / h', fontsize=fontsize_label)
ax3.set_xlim(5, 84)
ax3.set_ylim(-1, 8)
ax3.set_ylabel('Heat transfer rate / kW', fontsize=fontsize_label)
ax3.legend(loc = 'lower right',fontsize=fontsize_legend)
ax3.tick_params(labelsize=fontsize_ticks)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig('Figures/Fig_3.svg', dpi=300, bbox_inches='tight')
