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
# 1. TANK & CRYOGEN SETUP
# ---------------------------------------------------------
Q_roof = 0              # Roof heat ingress / W
T_air  = 18.08 + 273.15 # Temperature of the environment / K

V_tank = 4831           # Tank volume / m^3

a   = 0.5
d_i = ((4 * V_tank) / (np.pi * a))**(1/3) # Internal diameter / m

# Thickness in % of the internal diameter
ST  = 1.02
d_o = d_i * ST          # External diameter / m

# Set overall heat transfer coefficient through the walls for liquid and vapour
U_L = 8.72e-2           # W/m2/K
U_V = 8.72e-2           # W/m2/K

# Set overall bottom heat transfer coefficient respect the r_U
r_U = 4
U_b = r_U * U_L         # W/m2/K

# Initial liquid filling / Dimensionless
LF = 0.95

# Wall heat fraction
eta = 0.9

# Specify tank operating pressure
P = 101325 * 3          # Pa

# Initialise cryogen
cryogen = Cryogen(name="nitrogen")
cryogen.set_coolprops(P)

# Initialize large-scale tank
large_tank = Tank(d_i, d_o, V_tank, LF)
large_tank.set_HeatTransProps(U_L, U_V, T_air, q_b_fixed=None, Q_roof=0, eta_w=0.90)
large_tank.U_b = U_L * r_U

# Set cryogen
large_tank.cryogen = cryogen

print(f"The initial evaporation rate of {cryogen.name} is {large_tank.b_l_dot * 3600:.1f} kg/h")

# Define vertical spacing and computation grid
n_z = 100 
large_tank.z_grid = np.linspace(0, 1, n_z)

# Insulated roof
large_tank.U_roof = 0

# Time properties configuration
large_tank.time_interval = 1200
evap_time = 3600 * 24 * 30  # 30 days simulation execution runtime

# ---------------------------------------------------------
# 2. DATA LOADING AND SIMULATION
# ---------------------------------------------------------
folder_data = "../Results/Data/"         
folder_results = "../Results/Data/" 
os.makedirs(folder_results, exist_ok=True)

# 1. Read the CSV containing previously generated optimal configurations
df_opts = pd.read_csv(os.path.join(folder_data, 'LAES_opti_12h_LFs_ru_4_opts.csv'))

# 2. Extract specific optimal aspect ratios using np.isclose to handle float precision safely
opt_ar_95 = df_opts.loc[np.isclose(df_opts['LF'], 0.95, atol=1e-3), 'Optimal_Geometric_AR'].values[0]
opt_ar_05 = df_opts.loc[np.isclose(df_opts['LF'], 0.05, atol=1e-3), 'Optimal_Geometric_AR'].values[0]

print(f"Loaded optimal configurations -> LF 0.95: AR={opt_ar_95:.4f} | LF 0.05: AR={opt_ar_05:.4f}")

# 3. Initialize the JAX optimizer and simulation engine
coeffs = load_coolprop_coeffs('../cryoevap/cryogens/Coeffs/')
opti = TankOptimizerJAX(large_tank)

# 4. Iterate over both scenarios to simulate and export profiles
scenarios = [(0.95, opt_ar_95), (0.05, opt_ar_05)]

for lf_target, ar_target in scenarios:
    print(f"Simulating thermal profile for LF = {lf_target} over 30 days...")
    
    # Update tank properties and regenerate parameter mapping
    opti.tank.LF = lf_target
    opti.params = opti._build_params(opti.tank, coeffs)
    
    # Solve the system using the designated optimal aspect ratio
    sol = opti.simulate(aspect_ratio=ar_target, params=opti.params, t_final=evap_time)
    
    # Convert vapor temperature profile time series to DataFrame
    df_sol = pd.DataFrame(sol.ys[:, 1:])
    df_sol.insert(0, "Time (s)", sol.ts)
    
    # Export results to target folder
    lf_str = int(lf_target * 100) 
    filename = os.path.join(folder_results, f'T_V_LF{lf_str}_30days_ru_4.csv')
    df_sol.to_csv(filename, index=False)
    
    print(f"Successfully saved profile to: {filename}")