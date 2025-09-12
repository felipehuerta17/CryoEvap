# Libraries
import CoolProp.CoolProp as CP
from scipy.integrate import simpson
import numpy as np
import pandas as pd

# Calculate properties for each fluid
fluids   = ['ammonia', 'hydrogen', 'methane', 'nitrogen']
pressure = 100000 # Pa
ranges_temperature = [[250, 400], [100, 400], [120, 400], [100, 400]]

# Initialize lists to store coefficients
coeffs_kV_list = []
coeffs_rhoV_list = []
coeffs_cpV_list = []

for i in range(len(fluids)):
    T_vec = np.linspace(ranges_temperature[i][0], ranges_temperature[i][1], 200)
    rho_V = [CP.PropsSI('D','P', pressure,'T',T, fluids[i]) for T in T_vec]
    cp_V  = [CP.PropsSI('C','P', pressure,'T', T, fluids[i]) for T in T_vec]
    k_V   = [CP.PropsSI('L','P', pressure,'T', T, fluids[i]) for T in T_vec]

    coeffs_kV_list.append(np.polyfit(T_vec, k_V, 4))
    coeffs_rhoV_list.append(np.polyfit(T_vec, rho_V, 4))
    coeffs_cpV_list.append(np.polyfit(T_vec, cp_V, 6))

# Convert lists to DataFrames with fluids as columns and polynomial weights as rows
df_coeffs_kV = pd.DataFrame(np.array(coeffs_kV_list).T, columns=fluids)
df_coeffs_rhoV = pd.DataFrame(np.array(coeffs_rhoV_list).T, columns=fluids)
df_coeffs_cpV = pd.DataFrame(np.array(coeffs_cpV_list).T, columns=fluids)

# Rename rows to x^n format (polynomial degree)
df_coeffs_kV.index = [f'x^{i}' for i in range(len(df_coeffs_kV.index)-1, -1, -1)]
df_coeffs_rhoV.index = [f'x^{i}' for i in range(len(df_coeffs_rhoV.index)-1, -1, -1)]
df_coeffs_cpV.index = [f'x^{i}' for i in range(len(df_coeffs_cpV.index)-1, -1, -1)]

# Save to CSV
df_coeffs_kV.to_csv('../cryoevap/cryogens/Coeffs/coeffs_kV.csv')
df_coeffs_rhoV.to_csv('../cryoevap/cryogens/Coeffs/coeffs_rhoV.csv')
df_coeffs_cpV.to_csv('../cryoevap/cryogens/Coeffs/coeffs_cpV.csv')

