import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
files = ['results_4x4.pkl', 'results_8x8.pkl', 'results_16x16.pkl', 'results_32x32.pkl', 'results_64x64.pkl']
data = []
for i in range(len(files)):
    with open('Data/'+files[i], 'rb') as f:
        data.append(pkl.load(f))

time_days = data[0]['Time'] / (3600 * 24)
cmap = plt.get_cmap("inferno", len(files)+1)
colors = cmap(np.array([0, 1, 2, 3, 4]))

plt.figure(figsize=(12, 5), dpi=300)

ax1 = plt.subplot(1, 2, 1)
ax2 = plt.subplot(1, 2, 2)

for i, (file, color) in enumerate(zip(files, colors)):
    label = file.replace("results_", "").replace(".pkl", "")
    ax1.plot(time_days, data[i]["Tv_avg"], label=label, color=color)
    ax2.plot(time_days, data[i]["Tw_avg"], label=label, color=color)

    if i == 0:
        ref_tv = data[-1]["Tv_avg"]
        ref_tw = data[-1]["Tw_avg"]
        diff_profiles = {}

    print('MSE Tv:', np.mean((data[i]["Tv_avg"] - ref_tv) ** 2))
    print('MSE Tw:', np.mean((data[i]["Tw_avg"] - ref_tw) ** 2))

ax1.set_xlabel("Time / days")
ax1.set_ylabel("Vapour Average temperature / K")
ax1.grid(True, alpha=0.3)

ax2.set_xlabel("Time / days")
ax2.set_ylabel("Wall Average temperature / K")
ax2.grid(True, alpha=0.3)
handles, labels = ax1.get_legend_handles_labels()
plt.gcf().legend(
    handles,
    labels,
    title="Number of nodes (r x z)",
    loc="lower center",
    bbox_to_anchor=(0.5, 0.04),
    ncol=5,
    fontsize="small"
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.25)
plt.savefig("Figures/grid_comparison.svg", dpi=300)


def calculate_gci(f1, f2, f3, r=2.0):
    """
    Based on https://www.dynamore.eu/en/downloads/papers/dynamore/de/download/papers/forum08/dokumente/I-I-03.pdf
    """
    diff_32 = f3 - f2
    diff_21 = f2 - f1
    
    if diff_21 == 0:
        return np.inf, 0.0  
        
    p = np.log(abs(diff_32 / diff_21)) / np.log(r)
    
    f_exact = f2 - diff_32/(r**p - 1)

    error = abs((f_exact - f2)/f2)
    
    return p, error

# 1. Vapour Phase
tv_f1 = data[-1]["Tv_avg"][-1] 
tv_f2 = data[-2]["Tv_avg"][-1]
tv_f3 = data[-3]["Tv_avg"][-1]

p_vap, gci_vap = calculate_gci(tv_f1, tv_f2, tv_f3, r=2.0)

# 2. Wall Phase
tw_f1 = data[-1]["Tw_avg"][-1]
tw_f2 = data[-2]["Tw_avg"][-1]
tw_f3 = data[-3]["Tw_avg"][-1]

p_wall, gci_wall = calculate_gci(tw_f1, tw_f2, tw_f3, r=2.0)

print(f"Vapour Phase. Order of Convergence (p): {p_vap:.3f} | GCI: {gci_vap:.5f}%")
print(f"Wall Phase. Order of Convergence (p): {p_wall:.3f} | GCI: {gci_wall:.5f}%")
