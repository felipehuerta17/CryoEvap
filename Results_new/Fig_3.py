import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# ---------------------------------------------------------
# 1. CONFIGURACIÓN GENERAL
# ---------------------------------------------------------
folder = "Data/"

# Diccionario de configuración para automatizar la lectura y el ploteo
# Cada elemento de la lista corresponde a un subplot
configs = [
    {
        'ax_idx': 0,
        'title_letter': 'a)',
        'files': [
            {'path': 'T_V_LF5_30days_ru_025.csv', 'ru_label': '0.25', 'color_idx': 1},
            {'path': 'T_V_LF5_30days_ru_4.csv',   'ru_label': '4',    'color_idx': 2}
        ],
        'targets': [12, 72], # Horas objetivo para LF05
        'ylim': (86, 101), # Límite del eje Y para el primer subplot
        "text": r"LF$_0 = 0.05$" # Texto para el primer subplot
    },
    {
        'ax_idx': 1,
        'title_letter': 'b)',
        'files': [
            {'path': 'T_V_LF95_30days_ru_025.csv', 'ru_label': '0.25', 'color_idx': 1},
            {'path': 'T_V_LF95_30days_ru_4.csv',   'ru_label': '4',    'color_idx': 2}
        ],
        'targets': [1, 72], # Horas objetivo para LF95
        'ylim': (88, 88.65), # Límite del eje Y para el segundo subplot
        "text": r"LF$_0 = 0.95$" # Texto para el segundo subplot
    }
]

# Configuración visual
paleta = plt.get_cmap('inferno', 4)
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

# ---------------------------------------------------------
# 2. CREACIÓN DE FIGURA Y PLOTEO
# ---------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

for config in configs:
    ax = axs[config['ax_idx']]
    
    for file_info in config['files']:
        # Leer datos
        filepath = os.path.join(folder, file_info['path'])
        df = pd.read_csv(filepath)
        
        time_hours = df["Time (s)"] / 3600.0
        
        # Grid adimensional (todas las columnas menos la de Time)
        x = np.linspace(0, 1, len(df.columns) - 1)
        
        for i, target_time in enumerate(config['targets']):
            # Encuentra el índice con el tiempo más cercano al objetivo (más seguro que >=)
            idx = (time_hours - target_time).abs().idxmin()
            
            # Extraer el perfil de temperatura (todas las columnas desde la 1 en adelante)
            y = df.iloc[idx, 1:].to_numpy()
            
            # Estilo: línea continua para el primer target, punteada para el segundo
            linestyle = '-' if i == 0 else '--'
            label = rf'$r_U$ = {file_info["ru_label"]}, {target_time} h'
            
            ax.plot(x, y, linestyle=linestyle, label=label, 
                    color=paleta(file_info['color_idx']), linewidth=2)

    # Formato de cada subplot
    ax.set_xlabel(r'Dimensionless height / $\xi$', fontsize=fontsize_label)
    ax.tick_params(axis='both', labelsize=fontsize_ticks)
    ax.legend(fontsize=fontsize_legend)
    ax.text(0.03, 0.97, config['title_letter'], transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', fontname='Arial')
    ax.set_xlim(0, 1)
    ax.set_ylim(config['ylim']) # Ajusta el límite superior automáticamente
    ax.text(0.1, 0.97, config['text'], transform=ax.transAxes, 
            fontsize=16, va='top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'),)

# Añadir el label del eje Y solo al primer gráfico para mantenerlo limpio
axs[0].set_ylabel('Temperature / K', fontsize=fontsize_label)

plt.tight_layout()
os.makedirs("Figures", exist_ok=True)
plt.savefig("Figures/Fig_3.svg", bbox_inches='tight', dpi=300)
plt.close()