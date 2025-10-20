import pandas as pd
import matplotlib.pyplot as plt

# Importar datos
data        = pd.read_csv('Data/LAES_opti_12h_LFs.csv')
BOR         = data['Boil-Off Rate (BOR)']
geometrical = data['Aspect Ratio']
thermal     = data['Thermal Aspect Ratio']

# Optimal value
min_val   = BOR.min()              
min_index = BOR.idxmin()      

# Plot configuration
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth = 2
paleta = plt.get_cmap('inferno', 7)

# Create figure
plt.figure(figsize=(8, 6), dpi=300)
plt.plot(data['Aspect Ratio'], data['Boil-Off Rate (BOR)'], label='Geometrical', color=paleta(2), linewidth=linewidth)
plt.plot(data['Thermal Aspect Ratio'], data['Boil-Off Rate (BOR)'], label='Thermal', color=paleta(4), linewidth=linewidth)
plt.plot([-1, 3], [min_val, min_val], '--', color='gray', label = "Minimum BOR", linewidth=linewidth)
plt.scatter(geometrical[min_index], BOR[min_index], color=paleta(2), zorder=5)
plt.scatter(thermal[min_index], BOR[min_index], color=paleta(4), zorder=5)
plt.legend(loc='upper right', fontsize=fontsize_legend)
plt.ylabel('Boil-Off Rate / %/day', fontsize=fontsize_label)
plt.xlabel('Aspect Ratio', fontsize=fontsize_label)
plt.xlim(0, 2)
plt.tick_params(labelsize=fontsize_ticks)
# plt.grid(True, alpha = 0.5)

# Adjust layout to prevent overlap
plt.tight_layout()

# Save figure
plt.savefig('Figures/Fig_LAES1.svg', dpi=300, bbox_inches='tight')