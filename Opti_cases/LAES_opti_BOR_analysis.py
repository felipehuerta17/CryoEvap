import sys
import os
sys.path.append("..")

import numpy as np
import pandas as pd
import jax
import jax.numpy as jnp

from cryoevap.optimize import TankOptimizerJAX, load_coolprop_coeffs
from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

jax.config.update("jax_enable_x64", True)

# ---------------------------------------------------------
# TANK & CRYOGEN SETUP
# ---------------------------------------------------------
Q_roof = 0               
T_air  = 18.08 + 273.15  
V_tank = 4831            

U_L = 8.72e-2            
U_V = 8.72e-2            
P = 101325 * 3           

# Base tank initialization
a_init = 0.5
d_i_init = ((4 * V_tank) / (np.pi * a_init))**(1/3)

cryogen = Cryogen(name="nitrogen")
cryogen.set_coolprops(P)

large_tank = Tank(d_i_init, d_i_init * 1.02, V_tank, 0.95)
large_tank.set_HeatTransProps(U_L, U_V, T_air, q_b_fixed=None, Q_roof=0, eta_w=0.90)
large_tank.cryogen = cryogen

n_z = 50 
large_tank.z_grid = np.linspace(0, 1, n_z)
large_tank.U_roof = 0

evap_time = 3600 * 12
large_tank.time_interval = 1200
large_tank.plot_interval = evap_time / 6

# ---------------------------------------------------------
# OPTIMIZATION & DATA EXPORT FOR THE TABLE
# ---------------------------------------------------------
coeffs = load_coolprop_coeffs()
opti   = TankOptimizerJAX(large_tank)

folder = "../Results_new/Data/"
os.makedirs(folder, exist_ok=True)

# Parameters to fill the table
r_U_values = [0.25, 1.0, 4.0]
LF_values = [0.95, 0.50, 0.05] 

results_table = []

print("Calculating BOR for the r_U and LF combinations...")

for r_U in r_U_values:
    # Update the bottom heat transfer coefficient
    opti.tank.U_b = U_L * r_U
    
    for lf in LF_values:
        # Update the tank state
        opti.tank.LF = float(lf)
        opti.params = opti._build_params(opti.tank, coeffs)
        
        # Optimize the aspect ratio
        opt_ar, opt_tar, min_bor = opti.optimize(
            t_final=evap_time,
            coarse_samples=100,
            fine_samples=500,
            ar_min=0.05,
            ar_max=3.0
        )
        
        row_data = {
            "r_U": r_U,
            "LF": f"{lf*100:g}%",
            "Optimal a Value": float(opt_ar),
            "Optimal a BOR": float(min_bor)
        }
        
        # 4. Compute the BOR directly for a=0.5 and a=2.0
        bor_05 = opti._objective_bor(0.5, opti.params, evap_time)
        bor_20 = opti._objective_bor(2.0, opti.params, evap_time)
        
        row_data["BOR a=0.5"] = float(bor_05)
        row_data["BOR a=2.0"] = float(bor_20)
        
        results_table.append(row_data)
        print(f"Processed: r_U = {r_U}, LF = {lf*100:g}% | a_opt = {float(opt_ar):.4f}")

# ---------------------------------------------------------
# TABLE FORMATTING AND EXPORT
# ---------------------------------------------------------
df_results = pd.DataFrame(results_table)

# Save the raw results for backup
filename_raw = os.path.join(folder, 'BOR_table_raw.csv')
df_results.to_csv(filename_raw, index=False)

# Transpose and structure to match imagen.png
df_pivot = df_results.set_index(["r_U", "LF"]).T

filename_matrix = os.path.join(folder, 'BOR_table_matrix.csv')
df_pivot.to_csv(filename_matrix)

print(f"\nMatrix exported successfully to: {filename_matrix}")