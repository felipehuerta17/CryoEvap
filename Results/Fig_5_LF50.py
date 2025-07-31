import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

# Read pickle file
with open('Data/large_tank_data_LF50_3weeks.pkl', 'rb') as f:
    data_LF50 = pickle.load(f)

cmap = plt.get_cmap('inferno', 6)
T_env_avg = 5.3286+273.15     # K
T_range   = 15                # K

T_env = lambda t: T_env_avg + 0.5*T_range * np.sin(t * 2*np.pi/(24*3600))
plt.figure(figsize=(8, 6))
plt.plot(data_LF50['Time']/3600, T_env(data_LF50['Time']), label=r'$T_{\text{air}}$', color=cmap(0))
plt.plot(data_LF50['Time']/3600, data_LF50['Tw_avg'], label=r'$T_{\text{W,avg}}$', color=cmap(2))
plt.plot(data_LF50['Time']/3600, data_LF50['T_BOG'], label=r'$T_{\text{BOG}}$', color=cmap(3))
plt.plot(data_LF50['Time']/3600, data_LF50['Tv_avg'], label=r'$T_{\text{V,avg}}$', color=cmap(4))

plt.xlabel('Time / h', fontsize = fontsize_label)
plt.xticks(fontsize = fontsize_ticks)
plt.xlim(-1, 360)
plt.ylabel(r'Temperature / K', fontsize = fontsize_label)
plt.yticks(fontsize = fontsize_ticks)
plt.legend(loc =(1.05, 0.735), fontsize = fontsize_legend)
plt.savefig('Figures/Fig_5_LF50.svg', dpi=300, bbox_inches='tight')
