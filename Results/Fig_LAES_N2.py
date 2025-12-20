import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

folder = "Data/"
data_LF05_rq03 = pd.read_csv(folder + "T_V_LF05_30days_rq03.csv")
data_LF05_rq07 = pd.read_csv(folder + "T_V_LF05_30days_rq07.csv")
data_LF95_rq03 = pd.read_csv(folder + "T_V_LF95_30days_rq03.csv")
data_LF95_rq07 = pd.read_csv(folder + "T_V_LF95_30days_rq07.csv")
target_1_LF05 = 12 # hours
target_2_LF05 = 48 # hours

target_1_LF95 = 1 # hours
target_2_LF95 = 48 # hours

time_hours_LF05_rq03 = data_LF05_rq03["Time (s)"]/3600
time_hours_LF05_rq07 = data_LF05_rq07["Time (s)"]/3600
time_hours_LF95_rq07 = data_LF95_rq07["Time (s)"]/3600
time_hours_LF95_rq03 = data_LF95_rq03["Time (s)"]/3600

# Extract time at target points
# LF05
time_at_target_1_LF05_rq03 = time_hours_LF05_rq03[time_hours_LF05_rq03 >= target_1_LF05].idxmin()
time_at_target_2_LF05_rq03 = time_hours_LF05_rq03[time_hours_LF05_rq03 >= target_2_LF05].idxmin()
time_at_target_1_LF05_rq07 = time_hours_LF05_rq07[time_hours_LF05_rq07 >= target_1_LF05].idxmin()
time_at_target_2_LF05_rq07 = time_hours_LF05_rq07[time_hours_LF05_rq07 >= target_2_LF05].idxmin()
# LF95
time_at_target_1_LF95_rq03 = time_hours_LF95_rq03[time_hours_LF95_rq03 >= target_1_LF95].idxmin()
time_at_target_2_LF95_rq03 = time_hours_LF95_rq03[time_hours_LF95_rq03 >= target_2_LF95].idxmin()
time_at_target_1_LF95_rq07 = time_hours_LF95_rq07[time_hours_LF95_rq07 >= target_1_LF95].idxmin()
time_at_target_2_LF95_rq07 = time_hours_LF95_rq07[time_hours_LF95_rq07 >= target_2_LF95].idxmin()

# Adimensional grid
x = np.linspace(0, 1, len(data_LF05_rq03.columns)-1)

paleta = plt.get_cmap('inferno', 4)
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14

# Plot each target point
fig, axs = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

# first subplot: LF05
axs[0].plot(x, data_LF05_rq03.iloc[time_at_target_1_LF05_rq03, 1:], label=r'$r$ = 30%, 12 h', color=paleta(1), linewidth=2)
axs[0].plot(x, data_LF05_rq03.iloc[time_at_target_2_LF05_rq03, 1:], '--', label=r'$r$ = 30%, 48 h', color=paleta(1), linewidth=2)
axs[0].plot(x, data_LF05_rq07.iloc[time_at_target_1_LF05_rq07, 1:], label=r'$r$ = 70%, 12 h', color=paleta(2), linewidth=2)
axs[0].plot(x, data_LF05_rq07.iloc[time_at_target_2_LF05_rq07, 1:], '--', label=r'$r$ = 70%, 48 h', color=paleta(2), linewidth=2)
axs[0].set_xlabel(r'Dimensionless height / $\xi$', fontsize=fontsize_label)
axs[0].set_ylabel('Temperature / K', fontsize=fontsize_label)
axs[0].tick_params(axis='both', labelsize=fontsize_ticks)
axs[0].legend(fontsize=fontsize_legend)
axs[0].text(0.03, 0.97, 'a)', transform=axs[0].transAxes, fontsize=18, fontweight='bold', va='top', font='Arial')
axs[0].set_xlim(0, 1)
# second subplot: LF95
axs[1].plot(x, data_LF95_rq03.iloc[time_at_target_1_LF95_rq03, 1:], label=r'$r$ = 30%, 1 h', color=paleta(1), linewidth=2)
axs[1].plot(x, data_LF95_rq03.iloc[time_at_target_2_LF95_rq03, 1:], '--', label=r'$r$ = 30%, 48 h', color=paleta(1), linewidth=2)
axs[1].plot(x, data_LF95_rq07.iloc[time_at_target_1_LF95_rq07, 1:], label=r'$r$ = 70%, 1 h', color=paleta(2), linewidth=2)
axs[1].plot(x, data_LF95_rq07.iloc[time_at_target_2_LF95_rq07, 1:], '--', label=r'$r$ = 70%, 48 h', color=paleta(2), linewidth=2)
axs[1].set_xlabel(r'Dimensionless height / $\xi$', fontsize=fontsize_label)
axs[1].tick_params(axis='both', labelsize=fontsize_ticks)
axs[1].legend(fontsize=fontsize_legend)
axs[1].text(0.03, 0.97, 'b)', transform=axs[1].transAxes, fontsize=18, fontweight='bold', va='top', font='Arial')
axs[1].set_xlim(0, 1)
plt.tight_layout()
plt.savefig("Figures/Fig_LAES_N2.svg", bbox_inches='tight', dpi=300)
