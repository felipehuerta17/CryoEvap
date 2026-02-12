import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


folder = "Data/"
data = pd.read_csv(folder+'transient_period_LAES.csv')
# Configuration
fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth       = 2.5
paleta = plt.get_cmap('inferno', 4)

a = data['a']
a_relation = 1/(a**(1/3))
tau_rq_03_LF95 = data['tau_rq_03_LF95']/3600
tau_rq_07_LF95 = data['tau_rq_07_LF95']/3600
tau_rq_03_LF05 = data['tau_rq_03_LF05']/3600
tau_rq_07_LF05 = data['tau_rq_07_LF05']/3600

coef07_95 = np.polyfit(a_relation, tau_rq_07_LF95, 1)
coef03_95 = np.polyfit(a_relation, tau_rq_03_LF95, 1)
coef07_05 = np.polyfit(a_relation, tau_rq_07_LF05, 1)
coef03_05 = np.polyfit(a_relation, tau_rq_03_LF05, 1)

slope07_95 = coef07_95[0]
slope03_95 = coef03_95[0]
slope07_05 = coef07_05[0]
slope03_05 = coef03_05[0]


fig, axs = plt.subplots(1, 2, figsize=(15, 6), dpi=300)
# first subplot: LF05
axs[0].plot(a, tau_rq_03_LF05, label=r'$r$ = 30%', color=paleta(1), linewidth=linewidth)
axs[0].plot(a, tau_rq_07_LF05, label=r'$r$ = 70%', color=paleta(2), linewidth=linewidth)
axs[0].set_xlabel(r'Geometric aspect ratio ($a$)', fontsize=fontsize_label)
axs[0].set_ylabel(r'Transient period ($\tau$) / h', fontsize=fontsize_label)
axs[0].tick_params(axis='both', labelsize=fontsize_ticks)
axs[0].legend(prop={'size': fontsize_legend})
axs[0].text(0.03, 0.97, 'a)', transform=axs[0].transAxes, fontsize=18, fontweight='bold', va='top')
axs[0].set_ylim(0, 570)

# second subplot: LF95
axs[1].plot(a, tau_rq_03_LF95, label=r'$r$ = 30%', color=paleta(1), linewidth=linewidth)
axs[1].plot(a, tau_rq_07_LF95, label=r'$r$ = 70%', color=paleta(2), linewidth=linewidth)
axs[1].set_xlabel(r'Geometric aspect ratio ($a$)', fontsize=fontsize_label)
axs[1].tick_params(axis='both', labelsize=fontsize_ticks)
axs[1].legend(prop={'size': fontsize_legend})
axs[1].text(0.03, 0.97, 'b)', transform=axs[1].transAxes, fontsize=18, fontweight='bold', va='top')
axs[1].set_ylim(0, 28)

axs[0].text(0.83, 0.375, rf'$\tau = {slope03_05:.3f}\,a^{{-1/3}}$', transform=axs[0].transAxes,
            fontsize=14, ha='center', va='center', color=paleta(2))

axs[0].text(0.83, 0.18, rf'$\tau = {slope07_05:.3f}\,a^{{-1/3}}$', transform=axs[0].transAxes,
            fontsize=14, ha='center', va='center', color=paleta(1))

axs[1].text(0.83, 0.375, rf'$\tau = {slope03_95:.3f}\,a^{{-1/3}}$', transform=axs[1].transAxes,
            fontsize=14, ha='center', va='center', color=paleta(2))

axs[1].text(0.83, 0.18, rf'$\tau = {slope07_95:.3f}\,a^{{-1/3}}$', transform=axs[1].transAxes,
            fontsize=14, ha='center', va='center', color=paleta(1))


plt.tight_layout()
plt.savefig("Figures/Fig_LAES_N3.svg", bbox_inches='tight', dpi=300)
# plt.show()

