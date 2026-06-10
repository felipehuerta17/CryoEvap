import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# ---------------------------------------------------------
# 1. CARGA DE DATOS
# ---------------------------------------------------------
folder = "Data/"

data_files = {
    'ru_025': pd.read_csv(folder + 'LAES_opti_12h_LFs_ru_025.csv'),
    'ru_4': pd.read_csv(folder + 'LAES_opti_12h_LFs_ru_4.csv')
}

optimos_files = {
    'ru_025': pd.read_csv(folder + 'LAES_opti_12h_LFs_ru_025_opts.csv'),
    'ru_4': pd.read_csv(folder + 'LAES_opti_12h_LFs_ru_4_opts.csv')
}

# ---------------------------------------------------------
# 2. CONFIGURACIÓN VISUAL
# ---------------------------------------------------------
target_lfs = [0.05, 0.50, 0.95]

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth_curve = 3
linewidth_opt   = 3.5

colors = ['#4A1259', '#C54358', '#F99E1C'] 

# ---------------------------------------------------------
# 3. CREACIÓN DE FIGURA Y EJES QUEBRADOS
# ---------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(14, 6), dpi=300, 
                         gridspec_kw={'height_ratios': [1.2, 1], 'hspace': 0.08, 'wspace': 0.25})

# Ajuste de rangos: Nota cómo los yticks ahora están ESTRICTAMENTE dentro de los ylim
subplot_configs = [
    {
        'key': 'ru_025', 'xlim': [0.01, 0.5], 'letter': 'a)', 
        'ylim_bottom': [0.00, 0.008], 'yticks_bottom': [0.00, 0.005],
        'ylim_top': [0.025, 0.08],     'yticks_top': [0.03, 0.05, 0.07]
    },
    {
        'key': 'ru_4', 'xlim': [0.25, 3.0], 'letter': 'b)', 
        'ylim_bottom': [0.00, 0.015], 'yticks_bottom': [0.000, 0.005, 0.010, 0.015],
        'ylim_top': [0.06, 0.15],     'yticks_top': [0.08, 0.10, 0.12, 0.14]
    }
]

def plot_curves(ax, df_data, df_opts):
    for i, lf in enumerate(target_lfs):
        subset = df_data[np.isclose(df_data['Liquid_Filling'], lf, atol=1e-3)]
        ax.plot(subset['Geometric_AR'], subset['BOR'], 
                label=f"LF {lf:.2f}", linewidth=linewidth_curve, color=colors[i])
        
    ax.plot(df_opts['Optimal_Geometric_AR'], df_opts['Min_BOR'], 
            label='Optimal Values', linewidth=linewidth_opt, color='black')

for col_idx, config in enumerate(subplot_configs):
    ax_top = axes[0, col_idx]
    ax_bottom = axes[1, col_idx]
    
    df_data = data_files[config['key']]
    df_opts = optimos_files[config['key']]
    
    # 1. Graficar datos
    plot_curves(ax_top, df_data, df_opts)
    plot_curves(ax_bottom, df_data, df_opts)

    # 2. Aplicar límites exactos
    ax_top.set_xlim(config['xlim'])
    ax_bottom.set_xlim(config['xlim'])
    ax_top.set_ylim(config['ylim_top'])
    ax_bottom.set_ylim(config['ylim_bottom'])

    # 3. Aplicar Ticks explícitos
    ax_top.set_yticks(config['yticks_top'])
    ax_bottom.set_yticks(config['yticks_bottom'])
    
    # Formateo nativo de decimales (Reemplaza el antiguo set_yticklabels)
    if col_idx == 0:
        ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax_bottom.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    else:
        ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax_bottom.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    # 4. Magia del eje quebrado (Ocultar líneas centrales)
    ax_top.spines['bottom'].set_visible(False)
    ax_bottom.spines['top'].set_visible(False)
    
    # Quitar los ticks del eje X en el gráfico superior
    ax_top.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_bottom.tick_params(axis='x', which='both', bottom=True, top=False, labelsize=fontsize_ticks)
    ax_top.tick_params(axis='y', labelsize=fontsize_ticks)
    ax_bottom.tick_params(axis='y', labelsize=fontsize_ticks)

    # 5. Dibujar las líneas diagonales de corte (//)
    d = .015  # Inclinación/Largo de las marcas
    kwargs = dict(transform=ax_top.transAxes, color='black', clip_on=False, linewidth=1.5)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)        # Arriba-Izq
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # Arriba-Der
    kwargs.update(transform=ax_bottom.transAxes)  
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)        # Abajo-Izq
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # Abajo-Der

    # 6. Textos
    ax_bottom.set_xlabel('Geometrical Aspect Ratio', fontsize=fontsize_label)
    ax_top.text(0.04, 0.92, config['letter'], transform=ax_top.transAxes, 
                fontsize=16, fontweight='bold', va='top', fontname='Arial')

# Etiqueta global del eje Y compartida
fig.supylabel('Boil-Off Rate (BOR) / %/day', fontsize=fontsize_label, x=0.05, fontweight='medium')

# 7. Leyenda unificada al fondo
handles, labels = axes[1, 1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(), loc='lower center', 
           bbox_to_anchor=(0.5, -0.08), ncol=4, fontsize=fontsize_legend, frameon=True)

plt.savefig("Figures/Fig_2_V2.svg", bbox_inches='tight', dpi=300)
plt.close()