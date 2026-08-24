import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

a  = np.linspace(0.05, 2.5, 100)
LF = np.linspace(0.05, 0.95, 10)
opti_a = 1/(2*LF)
V_T = 1
cte = (np.pi/4)**(1/3)
print(cte)
# Set the seaborn palette
sns.set_palette("inferno", n_colors=len(LF))

plt.figure(figsize=(7, 5),)
# Plot each curve for different LF values
for lf in LF:
    ratio = cte*(lf*4*a**(1/3) + a**(-2/3))
    plt.plot(a, ratio, label=fr'LF$_0$={lf:.2f}')

# Update xlabel and ylabel with fontsize 12
ratios_opti = cte*(LF*4*opti_a**(1/3) + opti_a**(-2/3))
plt.plot(opti_a, ratios_opti, 'o--', markersize=8, label=fr'Optimal $a$')
plt.xlabel("Geometrical Aspect Ratio", fontsize=12)
plt.ylabel(r"Ratio $(A_L + A_b)/V_T$", fontsize=12)
plt.xlim(0, 2.5)
plt.legend(loc=(1.05, 0.41), fontsize=10)
plt.savefig('Figures/Fig_S2.svg', bbox_inches='tight')
plt.close()