import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# ---------------------------------------------------------
# General configuration
# ---------------------------------------------------------
folder = "Data/"

# Subplot configurations
configs = [
    {
        'ax_idx': 0,
        'title_letter': 'a)',
        'files': [
            {'path': 'T_V_LF5_30days_ru_025.csv', 'ru_label': '0.25', 'color_idx': 1},
            {'path': 'T_V_LF5_30days_ru_4.csv',   'ru_label': '4',    'color_idx': 2}
        ],
        'targets': [12, 72],  # Target hours
        'ylim': (86, 101),
        "text": r"LF$_0 = 0.05$"
    },
    {
        'ax_idx': 1,
        'title_letter': 'b)',
        'files': [
            {'path': 'T_V_LF95_30days_ru_025.csv', 'ru_label': '0.25', 'color_idx': 1},
            {'path': 'T_V_LF95_30days_ru_4.csv',   'ru_label': '4',    'color_idx': 2}
        ],
        'targets': [1, 72],  # Target hours
        'ylim': (88, 88.65),
        "text": r"LF$_0 = 0.95$"
    }
]

# Visual configuration
paleta = plt.get_cmap('inferno', 4)
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

# ---------------------------------------------------------
# Figure setup and plotting
# ---------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

for config in configs:
    ax = axs[config['ax_idx']]
    
    for file_info in config['files']:
        # Load data
        filepath = os.path.join(folder, file_info['path'])
        df = pd.read_csv(filepath)
        
        time_hours = df["Time (s)"] / 3600.0
        
        # Dimensionless grid
        x = np.linspace(0, 1, len(df.columns) - 1)
        
        for i, target_time in enumerate(config['targets']):
            # Closest time index
            idx = (time_hours - target_time).abs().idxmin()
            
            # Temperature profile
            y = df.iloc[idx, 1:].to_numpy()
            
            # Line style
            linestyle = '-' if i == 0 else '--'
            label = rf'$r_U$ = {file_info["ru_label"]}, {target_time} h'
            
            ax.plot(x, y, linestyle=linestyle, label=label, 
                    color=paleta(file_info['color_idx']), linewidth=2)

    # Subplot formatting
    ax.set_xlabel(r'Dimensionless height / $\zeta$', fontsize=fontsize_label)
    ax.tick_params(axis='both', labelsize=fontsize_ticks)
    ax.legend(fontsize=fontsize_legend)
    ax.text(0.03, 0.97, config['title_letter'], transform=ax.transAxes, 
            fontsize=18, fontweight='bold', va='top', fontname='Arial')
    ax.set_xlim(0, 1)
    ax.set_ylim(config['ylim'])
    ax.text(0.1, 0.97, config['text'], transform=ax.transAxes, 
            fontsize=16, va='top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'),)

# Shared y-axis label
axs[0].set_ylabel('Temperature / K', fontsize=fontsize_label)

plt.tight_layout()
os.makedirs("Figures", exist_ok=True)
plt.savefig("Figures/Fig_3.svg", bbox_inches='tight', dpi=300)
plt.close()