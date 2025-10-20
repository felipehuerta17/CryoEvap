import pandas as pd
import matplotlib.pyplot as plt

# Importar datos
data12h = pd.read_csv('Data/LAES_opti_12h_LFs.csv')
data24h = pd.read_csv('Data/LAES_opti_24h_LFs.csv')
data168h = pd.read_csv('Data/LAES_opti_168h_LFs.csv')

# Extraer variables relevantes
BOR_24h = data24h["BOR_LF_0.95"]
geom_AR_24h = data24h['Geometric_AR']

BOR_168h = data168h["BOR_LF_0.95"]
geom_AR_168h = data168h['Geometric_AR']

# Configuration
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth = 2.5
paleta = plt.get_cmap('inferno', 7)

# Create plot
plt.figure(figsize=(8, 6), dpi=300)
plt.plot(geom_AR_24h, BOR_24h, label=r'$t_{\text{final}} = $ 24 h', color=paleta(2), linewidth=linewidth)
plt.plot(geom_AR_168h, BOR_168h, label=r'$t_{\text{final}} = $ 168 h', color=paleta(5), linewidth=linewidth, linestyle=(0, (4, 4)))
plt.xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
plt.ylabel('Boil-Off Rate (BOR) / %/day', fontsize=fontsize_label)
plt.legend(fontsize=fontsize_legend, loc = 'upper right')
plt.tick_params(axis='both', labelsize=fontsize_ticks)
plt.xlim(0.2, 1)
plt.ylim(0.0029, 0.0035)
plt.tight_layout()

# Save Figure
plt.savefig("Figures/Fig_LAES2.svg", bbox_inches='tight', dpi=300)
