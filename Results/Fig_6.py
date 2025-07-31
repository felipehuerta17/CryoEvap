import pickle
import matplotlib.pyplot as plt

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

# Read pickle file
with open('Data/large_tank_data_LF50_3weeks.pkl', 'rb') as f:
    data_LF50 = pickle.load(f)

with open('Data/large_tank_data_LF95_3weeks.pkl', 'rb') as f:
    data_LF95 = pickle.load(f)

with open('Data/large_tank_data_LF50_3weeks_wo_wall.pkl', 'rb') as f:
    data_LF50_wo_wall = pickle.load(f)

with open('Data/large_tank_data_LF95_3weeks_wo_wall.pkl', 'rb') as f:
    data_LF95_wo_wall = pickle.load(f)

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), dpi=200)

# First subplot - BOG Temperature
cmap = plt.get_cmap('inferno', 4)
ax1.plot(data_LF50['Time']/3600, data_LF50['T_BOG'], label='LF=0.50', color=cmap(1))
ax1.plot(data_LF50_wo_wall['Time']/3600, data_LF50_wo_wall['T_BOG'],'--', 
         label='LF=0.50 (no wall model)', color=cmap(1))
ax1.plot(data_LF95['Time']/3600, data_LF95['T_BOG'], label='LF=0.95', color=cmap(2))
ax1.plot(data_LF95_wo_wall['Time']/3600, data_LF95_wo_wall['T_BOG'],'--', 
         label='LF=0.95 (no wall model)', color=cmap(2))
ax1.set_xlabel('Time / h', fontsize=14)
ax1.set_xlim(1, 500)
ax1.set_ylabel('Boil-off Gas Temperature / K', fontsize=fontsize_label)
ax1.tick_params(labelsize=fontsize_ticks)

# Second subplot - BOG Rate
ax2.plot(data_LF50['Time']/3600, data_LF50['BOG'], label='LF=0.50', color=cmap(1))
ax2.plot(data_LF50_wo_wall['Time']/3600, data_LF50_wo_wall['BOG'],'--', 
         label='LF=0.50 (no wall model)', color=cmap(1))
ax2.plot(data_LF95['Time']/3600, data_LF95['BOG'], label='LF=0.95', color=cmap(2))
ax2.plot(data_LF95_wo_wall['Time']/3600, data_LF95_wo_wall['BOG'],'--', 
         label='LF=0.95 (no wall model)', color=cmap(2))
ax2.set_xlabel('Time / h', fontsize=14)
ax2.set_xlim(1, 500)
ax2.set_ylim(0, 0.01)
ax2.set_ylabel(r'Boil-off Rate / $kg\ h^{-1}$', fontsize=fontsize_label)
ax2.tick_params(labelsize=fontsize_ticks)

ax2.legend(loc=(1.05, 0.71), fontsize=fontsize_legend)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig('Figures/Fig_6.svg', dpi=300, bbox_inches='tight')  