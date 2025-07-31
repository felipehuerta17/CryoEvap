
import pickle
import matplotlib.pyplot as plt
import numpy as np

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14


# Read pickle file
with open('Data/large_tank_data_LF95_3weeks.pkl', 'rb') as f:
    data_LF95 = pickle.load(f)

with open('Data/large_tank_data_LF95_3weeks_wo_wall.pkl', 'rb') as f:
    data_LF95_wo_wall = pickle.load(f)

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

# First subplot - Total, Liquid and Bottom heat transfer rates
cmap = plt.get_cmap('inferno', 5)
Q_b_no_wall_model = np.float64(1107.5268769106317)
Q_b_no_wall_vec = np.full_like(data_LF95['Time'], Q_b_no_wall_model)

ax1.plot(data_LF95['Time']/3600, data_LF95['Q_tot']/1000, label=r'$\dot{Q}_{\text{tot}}$', color=cmap(3))
ax1.plot(data_LF95_wo_wall['Time']/3600, data_LF95_wo_wall['Q_tot']/1000,'--', 
         label=r'$\dot{Q}_{\text{tot}}$ (no wall model)', color=cmap(3))

ax1.plot(data_LF95['Time']/3600, data_LF95['Q_L']/1000, label=r'$\dot{Q}_{\text{L}}$', color=cmap(2))
ax1.plot(data_LF95_wo_wall['Time']/3600, data_LF95_wo_wall['Q_L']/1000,'--', 
         label=r'$\dot{Q}_{\text{L}}$ (no wall model)', color=cmap(2))

ax1.plot(data_LF95['Time']/3600, data_LF95['Q_b']/1000, label=r'$\dot{Q}_{\text{b}}$', color=cmap(1))
ax1.plot(data_LF95_wo_wall['Time']/3600, Q_b_no_wall_vec/1000,'--', 
         label=r'$\dot{Q}_{\text{b}}$ (no wall model)', color=cmap(1))

ax1.set_xlim(0, 500)
ax1.set_ylim(0, 12)
ax1.set_ylabel(r'Heat Transfer Rate / kW', fontsize=fontsize_label)
ax1.set_xlabel('Time / h', fontsize=fontsize_label)
ax1.tick_params(labelsize=fontsize_ticks)
ax1.legend(fontsize=fontsize_legend)

# Second subplot - Wall and Vapor-Liquid Interface heat transfer rates
ax2.plot(data_LF95['Time']/3600, data_LF95['Q_Vw']/1000, label=r'$\dot{Q}_{\text{Wi}}$', color=cmap(3))
ax2.plot(data_LF95_wo_wall['Time']/3600, data_LF95_wo_wall['Q_Vw']/1000,'--', 
         label=r'$\dot{Q}_{\text{Wi}}$ (no wall model)', color=cmap(3))

ax2.plot(data_LF95['Time']/3600, data_LF95['Q_VL']/1000, label=r'$\dot{Q}_{\text{VL}}$', color=cmap(1))
ax2.plot(data_LF95_wo_wall['Time']/3600, data_LF95_wo_wall['Q_VL']/1000,'--', 
         label=r'$\dot{Q}_{\text{VL}}$ (no wall model)', color=cmap(1))

ax2.set_xlim(0, 500)
ax2.set_ylim(0, 0.25)
ax2.set_ylabel(r'Heat Transfer Rate / kW', fontsize=fontsize_label)
ax2.set_xlabel('Time / h', fontsize=fontsize_label)
ax2.tick_params(labelsize=fontsize_ticks)
ax2.legend(fontsize=fontsize_legend)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig('Figures/Fig_4.svg', dpi=300, bbox_inches='tight')
