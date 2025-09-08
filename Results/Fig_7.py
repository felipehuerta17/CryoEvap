import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d


fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

# Read pickle file
with open('Data/large_tank_data_LF95_3weeks.pkl', 'rb') as f:
    data_LF95 = pickle.load(f)
# Read pickle file
with open('Data/large_tank_data_LF50_3weeks.pkl', 'rb') as f:
    data_LF50 = pickle.load(f)

with open('Data/large_tank_data_LF50_3weeks_wo_wall.pkl', 'rb') as f:
    data_LF50_wo_wall = pickle.load(f)

with open('Data/large_tank_data_LF95_3weeks_wo_wall.pkl', 'rb') as f:
    data_LF95_wo_wall = pickle.load(f)
# Create figure with subplots

# First subplot - Total, Liquid and Bottom heat transfer rates (normalized)
cmap = plt.get_cmap('inferno', 6)

# Mask to select data starting from a specific time
time_tr = 240
time_to = 500
mask1 = data_LF95['Time'] <= time_tr*3600
mask2 = (data_LF95['Time'] >= time_tr*3600) & (data_LF95['Time'] < (time_to+time_tr)*3600)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), dpi=200)

# Crear figura y eje
# fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
ax1.plot(data_LF50['Time'][mask2]/3600, data_LF50['Q_tot'][mask2]/1000, label=r'$\dot{Q}_{\text{tot}}$', color=cmap(0))
ax1.plot(data_LF50['Time'][mask2]/3600, data_LF50['Q_L'][mask2]/1000, label=r'$\dot{Q}_{\text{Li}}$', color=cmap(4))
ax1.plot(data_LF50['Time'][mask2]/3600, data_LF50['Q_b'][mask2]/1000, label=r'$\dot{Q}_{\text{b}}$', color=cmap(3))
ax1.plot(data_LF50['Time'][mask2]/3600, data_LF50['Q_Vw'][mask2]/1000, label=r'$\dot{Q}_{\text{Wi}}$', color=cmap(2))

ax1.set_xlim(time_tr,time_to)
ax1.set_ylabel(r'Heat transfer rate / kW', fontsize=fontsize_label)
ax1.set_ylim(0, 50)
ax1.set_xlabel('Time / h', fontsize=fontsize_label)
ax1.tick_params(labelsize=fontsize_ticks)
ax1.legend(fontsize=fontsize_legend, loc='upper right', ncol=1, framealpha=1)

# Second subplot - BOG Temperature
cmap = plt.get_cmap('inferno', 4)

time_tr = 168
init_stationary = 239
end_stationary = 500 


mask_transient_LF50  = data_LF50['Time'] < time_tr *3600
mask_stationary_LF50 = (data_LF50['Time'] >= init_stationary *3600) & (data_LF50['Time'] < (end_stationary) * 3600)
mask_transient_LF95  = data_LF95['Time'] < time_tr *3600
mask_stationary_LF95 = (data_LF95['Time'] >= init_stationary * 3600) & (data_LF95['Time'] < (end_stationary) * 3600)
mask_transient_LF50_wo_wall  = data_LF50_wo_wall['Time'] < time_tr *3600
mask_stationary_LF50_wo_wall = (data_LF50_wo_wall['Time'] >= init_stationary *3600) & (data_LF50_wo_wall['Time'] < (end_stationary) * 3600)
mask_transient_LF95_wo_wall  = data_LF95_wo_wall['Time'] < time_tr *3600
mask_stationary_LF95_wo_wall = (data_LF95_wo_wall['Time'] >= init_stationary *3600) & (data_LF95_wo_wall['Time'] < (end_stationary) * 3600)

# Second subplot - stationary period
ax2.plot(data_LF50['Time'][mask_stationary_LF50]/3600, data_LF50['BOG'][mask_stationary_LF50]*3600, label='Wall model', color=cmap(2))

# Smooth the BOG data because of high frequency noise
bog_smooth = uniform_filter1d(
    data_LF50_wo_wall['BOG'][mask_stationary_LF50_wo_wall]*3600, 
    size=30
)
ax2.plot(
    data_LF50_wo_wall['Time'][mask_stationary_LF50_wo_wall]/3600,
    bog_smooth,
    label='Simplified (smoothed)', color=cmap(1)
)


ax2.set_xlabel('Time / h', fontsize=14)
ax2.set_xlim(init_stationary+1,end_stationary)
ax2.set_ylabel(r'Boil-off gas rate / $kg\ h^{-1}$', fontsize=fontsize_label)
# ax2.set_ylim(7, 12)
ax2.set_ylabel(r'Boil-off gas rate / $kg\ h^{-1}$', fontsize=fontsize_label)
ax2.tick_params(labelsize=fontsize_ticks)
ax2.legend(loc='upper right', fontsize=fontsize_legend, ncol=2, columnspacing=1.0)
plt.tight_layout()

plt.savefig('Figures/Fig_7.svg', dpi=300, bbox_inches='tight')

