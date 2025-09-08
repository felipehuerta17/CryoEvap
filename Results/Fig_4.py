
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

with open('Data/large_tank_data_LF95_3weeks.pkl', 'rb') as f:
    data_LF95 = pickle.load(f)

T_w_csv = pd.DataFrame(data_LF95['T_w_raw'])

# Create figure with subplots
# First subplot - Early hours temperature profiles
e      = 10.6972*2 - 10.24*2 # Wall thickness m
V_tank = 88023.6952 #m^3
a      = 0.5 # Aspect ratio
d_i    = ((4 * V_tank)/(np.pi * a))**(1/3) # internal diameter / m
d_o    = d_i + e # external diameter / m
r_i    = d_i/2
r_o    = d_o/2
r_adim = np.linspace(0, 1, T_w_csv.shape[0])
r_grid = r_adim * (r_o - r_i) + r_i


# Second subplot - Later hours temperature profiles
horas_objetivo2 = [50.0, 55.0, 60.0, 65.0, 70.0]
selected_index2 = [np.where(np.isclose(data_LF95['Time']/3600, h))[0][0] for h in horas_objetivo2]
cmap2 = plt.get_cmap('inferno', len(selected_index2)+1)
plt.figure(figsize=(7, 6), dpi=300)

for i in range(len(selected_index2)):
    plt.plot(r_grid, T_w_csv.iloc[:, selected_index2[i]], color=cmap2(i), 
             label=f't = {data_LF95["Time"][selected_index2[i]]/3600:.0f} h')
plt.legend(loc='lower right', fontsize=fontsize_legend)
plt.ylabel('Wall temperature / K', fontsize=fontsize_label)
plt.xlabel('Radius / m', fontsize=fontsize_label)
plt.xlim(r_i*0.999, r_o*1.001)
plt.tick_params(labelsize=fontsize_ticks)


# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig('Figures/Fig_4.svg', dpi=300, bbox_inches='tight')