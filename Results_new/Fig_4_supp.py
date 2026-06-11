import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# Configuration
folder = "Data/"
data = pd.read_csv(folder + 'transient_period.csv')

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 11  
linewidth       = 2.5

a = data['a']
a_relation = 1 / (a**(1/3))

# 1. Búsqueda dinámica de columnas y valores r_U
col_map_05 = {}
col_map_95 = {}

for col in data.columns:
    # Busca el patrón 'tau_ru_NUMERO_LF'
    match = re.search(r'tau_ru_([0-9.]+)_LF', col)
    if match:
        ru_val = float(match.group(1))
        if 'LF_0.05' in col:
            col_map_05[ru_val] = col
        elif 'LF_0.95' in col:
            col_map_95[ru_val] = col

# Obtener todos los valores únicos de r_U ordenados de menor a mayor
ru_values = sorted(list(set(list(col_map_05.keys()) + list(col_map_95.keys()))))
n_curves = len(ru_values)

# Identificar el menor y mayor r_U
min_ru = ru_values[0]
max_ru = ru_values[-1]

# Paleta plasma ajustada a la cantidad de curvas exactas
paleta = plt.get_cmap('plasma', n_curves)

fig, axs = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

# 2. Bucle principal para graficar y calcular regresiones
for i, ru in enumerate(ru_values):
    color = paleta(i)
    
    # --- Subplot 1: LF = 0.05 ---
    if ru in col_map_05:
        col_name = col_map_05[ru]
        tau = data[col_name] / 3600
        
        coef = np.polyfit(a_relation, tau, 1)
        r2 = np.corrcoef(a_relation, tau)[0, 1]**2
        print(f"R^2 (r_U={ru:<5}, LF=0.05): {r2:.4f}")
        
        # Etiqueta solo si es el min o max r_U, de lo contrario None
        if ru == min_ru or ru == max_ru:
            label = rf'$r_U$ = {ru} ($\tau \sim {coef[0]:.1f}\,a^{{-1/3}}$)'
        else:
            label = None
            
        axs[0].plot(a, tau, label=label, color=color, linewidth=linewidth)

    # --- Subplot 2: LF = 0.95 ---
    if ru in col_map_95:
        col_name = col_map_95[ru]
        tau = data[col_name] / 3600
        
        coef = np.polyfit(a_relation, tau, 1)
        r2 = np.corrcoef(a_relation, tau)[0, 1]**2
        print(f"R^2 (r_U={ru:<5}, LF=0.95): {r2:.4f}")
        
        # Etiqueta solo si es el min o max r_U, de lo contrario None
        if ru == min_ru or ru == max_ru:
            label = rf'$r_U$ = {ru} ($\tau \sim {coef[0]:.1f}\,a^{{-1/3}}$)'
        else:
            label = None
            
        axs[1].plot(a, tau, label=label, color=color, linewidth=linewidth)

# 3. Formato del Subplot a) LF = 0.05
axs[0].set_xlabel(r'Geometric aspect ratio ($a$)', fontsize=fontsize_label)
axs[0].set_ylabel(r'Transient period ($\tau$) / h', fontsize=fontsize_label)
axs[0].tick_params(axis='both', labelsize=fontsize_ticks)
axs[0].legend(prop={'size': fontsize_legend})
axs[0].text(0.03, 0.97, 'a)', transform=axs[0].transAxes, fontsize=18, fontweight='bold', va='top')
axs[0].text(0.1, 0.97, r"LF$_0 = 0.05$", transform=axs[0].transAxes, 
            fontsize=16, va='top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
axs[0].set_ylim(80,500)  # Ajustar el límite inferior a 0
# 4. Formato del Subplot b) LF = 0.95
axs[1].set_xlabel(r'Geometric aspect ratio ($a$)', fontsize=fontsize_label)
axs[1].tick_params(axis='both', labelsize=fontsize_ticks)
axs[1].legend(prop={'size': fontsize_legend})
axs[1].text(0.03, 0.97, 'b)', transform=axs[1].transAxes, fontsize=18, fontweight='bold', va='top')
axs[1].text(0.1, 0.97, r"LF$_0 = 0.95$", transform=axs[1].transAxes, 
            fontsize=16, va='top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
axs[1].set_ylim(3,25)  # Ajustar el límite inferior a 0
plt.tight_layout()
plt.savefig("Figures/Fig_4_supp.svg", bbox_inches='tight', dpi=300)
plt.close()