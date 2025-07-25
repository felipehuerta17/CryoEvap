import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Read pickle file
with open('Data/large_tank_data_LF95_3weeks.pkl', 'rb') as f:
    data_LF95 = pickle.load(f)

cmap = plt.get_cmap('inferno', 6)
T_env_avg = 5.3286+273.15     # K
T_range   = 15                # K

T_env = lambda t: T_env_avg + 0.5*T_range * np.sin(t * 2*np.pi/(24*3600))
plt.figure(figsize=(8, 6))
plt.plot(data_LF95['Time']/3600, T_env(data_LF95['Time']), label=r'$T_{\text{air}}$', color=cmap(0))
plt.plot(data_LF95['Time']/3600, data_LF95['Tw_avg'], label=r'$T_{\text{W,avg}}$', color=cmap(2))
plt.plot(data_LF95['Time']/3600, data_LF95['T_BOG'], label=r'$T_{\text{BOG}}$', color=cmap(3))
plt.plot(data_LF95['Time']/3600, data_LF95['Tv_avg'], label=r'$T_{\text{V,avg}}$', color=cmap(4))

plt.xlabel('Time / h', fontsize = 14)
plt.xticks(fontsize = 12)
plt.xlim(-1, 180)
plt.ylabel(r'Temperature / K', fontsize = 14)
plt.yticks(fontsize = 12)
plt.legend(loc =(1.05, 0.77), fontsize = 12)
plt.savefig('Figures/Fig_5.svg', dpi=300, bbox_inches='tight')
