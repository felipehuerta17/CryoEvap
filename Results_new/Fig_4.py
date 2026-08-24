import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Configuration
folder = "Data/"
data = pd.read_csv(folder + 'transient_period.csv')

fontsize_label  = 14
fontsize_ticks  = 14
fontsize_legend = 14
linewidth       = 2.5
paleta = plt.get_cmap('inferno', 4)

a = data['a']
a_relation = 1 / (a**(1/3))

tau_rq_03_LF95 = data['tau_ru_0.25_LF_0.95'] / 3600
tau_rq_07_LF95 = data['tau_ru_4.00_LF_0.95'] / 3600
tau_rq_03_LF05 = data['tau_ru_0.25_LF_0.05'] / 3600
tau_rq_07_LF05 = data['tau_ru_4.00_LF_0.05'] / 3600

# 1. Ajuste lineal (Regresión)
coef03_95 = np.polyfit(a_relation, tau_rq_03_LF95, 1)
coef07_95 = np.polyfit(a_relation, tau_rq_07_LF95, 1)
coef03_05 = np.polyfit(a_relation, tau_rq_03_LF05, 1)
coef07_05 = np.polyfit(a_relation, tau_rq_07_LF05, 1)

# 2. Evaluación de las curvas ajustadas
fit_03_95 = np.polyval(coef03_95, a_relation)
fit_07_95 = np.polyval(coef07_95, a_relation)
fit_03_05 = np.polyval(coef03_05, a_relation)
fit_07_05 = np.polyval(coef07_05, a_relation)

# Calculate R^2 values
r2_03_95 = np.corrcoef(a_relation, tau_rq_03_LF95)[0, 1]**2
r2_07_95 = np.corrcoef(a_relation, tau_rq_07_LF95)[0, 1]**2
r2_03_05 = np.corrcoef(a_relation, tau_rq_03_LF05)[0, 1]**2
r2_07_05 = np.corrcoef(a_relation, tau_rq_07_LF05)[0, 1]**2

print(f"R^2 (r_U=0.25, LF=0.95): {r2_03_95:.4f}")
print(f"R^2 (r_U=4.00, LF=0.95): {r2_07_95:.4f}")
print(f"R^2 (r_U=0.25, LF=0.05): {r2_03_05:.4f}")
print(f"R^2 (r_U=4.00, LF=0.05): {r2_07_05:.4f}")

fig, axs = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

axs[0].plot(a, tau_rq_03_LF05, label=r'$r_U$ = 0.25', color=paleta(1), linewidth=linewidth)
axs[0].plot(a, tau_rq_07_LF05, label=r'$r_U$ = 4.00', color=paleta(2), linewidth=linewidth)
axs[0].set_xlabel(r'Geometric aspect ratio ($a$)', fontsize=fontsize_label)
axs[0].set_ylabel(r'Transient period ($\tau$) / h', fontsize=fontsize_label)
axs[0].tick_params(axis='both', labelsize=fontsize_ticks)
axs[0].legend(prop={'size': fontsize_legend})
axs[0].text(0.03, 0.97, 'a)', transform=axs[0].transAxes, fontsize=18, fontweight='bold', va='top')
axs[0].set_ylim(0, 570)

eq_03_05 = rf'$\tau = {coef03_05[0]:.2f}\,a^{{-1/3}} + {coef03_05[1]:.2f}$'
eq_07_05 = rf'$\tau = {coef07_05[0]:.2f}\,a^{{-1/3}} + {coef07_05[1]:.2f}$'

axs[0].text(0.6, 0.53, eq_03_05, transform=axs[0].transAxes, fontsize=14, color=paleta(1)) # Corrección de color y posición
axs[0].text(0.565, 0.34, eq_07_05, transform=axs[0].transAxes, fontsize=14, color=paleta(2))
axs[0].text(0.1, 0.97, r"LF$_0 = 0.05$", transform=axs[0].transAxes, 
            fontsize=16, va='top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'),)


axs[1].plot(a, tau_rq_03_LF95, label=r'$r_U$ = 0.25', color=paleta(1), linewidth=linewidth)
axs[1].plot(a, tau_rq_07_LF95, label=r'$r_U$ = 4.00', color=paleta(2), linewidth=linewidth)
axs[1].set_xlabel(r'Geometric aspect ratio ($a$)', fontsize=fontsize_label)
axs[1].tick_params(axis='both', labelsize=fontsize_ticks)
axs[1].legend(prop={'size': fontsize_legend})
axs[1].text(0.03, 0.97, 'b)', transform=axs[1].transAxes, fontsize=18, fontweight='bold', va='top')
axs[1].set_ylim(0, 28)
axs[1].text(0.1, 0.97, r"LF$_0 = 0.95$", transform=axs[1].transAxes, 
            fontsize=16, va='top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'),)

# Textos de ecuaciones con colores corregidos
eq_03_95 = rf'$\tau = {coef03_95[0]:.2f}\,a^{{-1/3}} + {coef03_95[1]:.2f}$'
eq_07_95 = rf'$\tau = {coef07_95[0]:.2f}\,a^{{-1/3}} + {coef07_95[1]:.2f}$'

axs[1].text(0.635, 0.53, eq_03_95, transform=axs[1].transAxes, fontsize=14, color=paleta(1))
axs[1].text(0.6, 0.34, eq_07_95, transform=axs[1].transAxes, fontsize=14, color=paleta(2))

plt.tight_layout()
plt.savefig("Figures/Fig_4.svg", bbox_inches='tight', dpi=300)
plt.close()