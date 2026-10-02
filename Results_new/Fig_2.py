import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
plt.rcParams['mathtext.fontset'] = 'cm'

# ---------------------------------------------------------
# Data loading
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

analytical_files = {
    'ru_025': pd.read_csv(folder + 'analytical_optimal_r_U_0.25.csv'),
    'ru_4': pd.read_csv(folder + 'analytical_optimal_r_U_4.0.csv')
}

# ---------------------------------------------------------
# Visual configuration
# ---------------------------------------------------------
target_lfs = [0.05, 0.50, 0.95]

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth_curve = 3
linewidth_opt   = 3.5

colors = ['#4A1259', '#C54358', '#F99E1C'] 

# ---------------------------------------------------------
# Figure and broken axes setup
# ---------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(14, 6), dpi=300, 
                         gridspec_kw={'height_ratios': [1.2, 1], 'hspace': 0.08, 'wspace': 0.2})

# Subplot configurations
subplot_configs = [
    {
        'key': 'ru_025', 'xlim': [0.01, 0.5], 'letter': 'a)', "text": r"$r_U = 0.25$",
        'ylim_bottom': [0.00, 0.8], 'yticks_bottom': [0.00, 0.5],
        'ylim_top': [2.5, 8.0],     'yticks_top': [2.5, 5.0, 7.5]
    },
    {
        'key': 'ru_4', 'xlim': [0.25, 3.0], 'letter': 'b)', "text": r"$r_U = 4.00$",
        'ylim_bottom': [0.00, 1.5], 'yticks_bottom': [0.000, 0.5, 1.0, 1.5],
        'ylim_top': [6.0, 15.0],     'yticks_top': [8, 10, 12, 14]
    }
]

def plot_curves(ax, df_data, df_opts, df_analytical, subplot_idx):
    for i, lf in enumerate(target_lfs):
        subset = df_data[np.isclose(df_data['LF'], lf, atol=1e-3)]
        ax.plot(subset['Geometric_AR'], subset['BOR']*100, 
                label=rf"LF$_0$ = {lf:.2f}", linewidth=linewidth_curve, color=colors[i])
        
    ax.plot(df_opts['Optimal_Geometric_AR'], df_opts['Min_BOR']*100, 
            label='Optimal Values', linewidth=linewidth_opt, color='black')
    
    # Geometric sampling
    n_points = 14
    sample_idx = np.unique(np.geomspace(1, len(df_analytical) - 1, n_points).astype(int))
    sample_idx = np.insert(sample_idx, 0, 0)
    
    if subplot_idx == 0:
        sample_idx = np.delete(sample_idx, 2)
    ax.plot(df_analytical['Optimal_AR'].iloc[sample_idx],
        df_analytical['BOR'].iloc[sample_idx] * 100, "x", markersize=10,
        markeredgewidth=2, label='Analytical Optimal Values', color="red")
            
    
for col_idx, config in enumerate(subplot_configs):
    ax_top = axes[0, col_idx]
    ax_bottom = axes[1, col_idx]
    
    df_data = data_files[config['key']]
    df_opts = optimos_files[config['key']]
    df_analytical = analytical_files[config['key']]

    # Plot curves
    plot_curves(ax_top, df_data, df_opts, df_analytical, 0)
    plot_curves(ax_bottom, df_data, df_opts, df_analytical, 1)

    # Axis limits
    ax_top.set_xlim(config['xlim'])
    ax_bottom.set_xlim(config['xlim'])
    ax_top.set_ylim(config['ylim_top'])
    ax_bottom.set_ylim(config['ylim_bottom'])

    # Axis ticks
    ax_top.set_yticks(config['yticks_top'])
    ax_bottom.set_yticks(config['yticks_bottom'])
    
    # Tick formatting
    if col_idx == 0:
        ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax_bottom.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    else:
        ax_top.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax_bottom.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    # Broken axis styling
    ax_top.spines['bottom'].set_visible(False)
    ax_bottom.spines['top'].set_visible(False)
    
    # Remove x-axis ticks on top subplot
    ax_top.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_bottom.tick_params(axis='x', which='both', bottom=True, top=False, labelsize=fontsize_ticks)
    ax_top.tick_params(axis='y', labelsize=fontsize_ticks)
    ax_bottom.tick_params(axis='y', labelsize=fontsize_ticks)

    # Diagonal cut marks
    d = .015  # Marker size
    kwargs = dict(transform=ax_top.transAxes, color='black', clip_on=False, linewidth=1.5)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)        # Top-left
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # Top-right
    kwargs.update(transform=ax_bottom.transAxes)  
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)        # Bottom-left
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # Bottom-right

    # Labels and annotations
    ax_bottom.set_xlabel(r"Geometrical Aspect Ratio (   )", fontsize=fontsize_label)
    ax_bottom.xaxis.label.set_position([0.5, -0.12]) 

    ax_bottom.text(0.7558, -0.155, r"$a$", transform=ax_bottom.transAxes,
                fontsize=fontsize_label + 3,  
                ha='center', va='top')

    ax_top.text(0.04, 0.92, config['letter'], transform=ax_top.transAxes, 
                fontsize=16, fontweight='bold', va='top', fontname='Arial')

    ax_top.text(0.75, 0.92, config['text'], transform=ax_top.transAxes, 
            fontsize=16, fontweight='bold', va='top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'),)
# Shared y-axis label
fig.supylabel('Boil-Off Rate (BOR) / %/day', fontsize=fontsize_label, x=0.065, fontweight='medium')

# Legend
handles, labels = axes[1, 1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(), loc='lower center', 
           bbox_to_anchor=(0.5, -0.08), ncol=5, fontsize=fontsize_legend, frameon=True)

plt.savefig("Figures/Fig_2.svg", bbox_inches='tight', dpi=300)
plt.close()