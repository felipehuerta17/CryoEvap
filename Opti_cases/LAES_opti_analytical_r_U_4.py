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

# Analytical optimal aspect ratio function
a_opt_analytical = lambda ru, lf: (ru) / (2 *(lf + (1-lf)*eta))

# ---------------------------------------------------------
# TANK & CRYOGEN SETUP
# ---------------------------------------------------------
Q_roof = 0               # Roof heat ingress / W
T_air  = 18.08 + 273.15  # Temperature of the environment / K

V_tank = 4831            # Tank volume / m^3

a   = 0.5
d_i = ((4 * V_tank) / (np.pi * a))**(1/3) # Internal diameter / m

# Thickness in % of the internal diameter
ST  = 1.02
d_o = d_i * ST           # External diameter / m

# Overall heat transfer coefficient through walls
U_L = 8.72e-2            # W/m2/K
U_V = 8.72e-2            # W/m2/K

# Overall bottom heat transfer coefficient respect to r_U
r_U = 4
U_b = r_U * U_L          # W/m2/K

# Initial liquid filling / Dimensionless
LF = 0.95

# Wall heat fraction
eta = 0.9

# Specify tank operating pressure
P = 101325 * 3           # Pa

# Initialise cryogen
cryogen = Cryogen(name="nitrogen")
cryogen.set_coolprops(P)

# Initialize large-scale tank
large_tank = Tank(d_i, d_o, V_tank, LF)
large_tank.set_HeatTransProps(U_L, U_V, T_air, q_b_fixed=None, Q_roof=0, eta_w=0.90)
large_tank.U_b = U_L * r_U

# Set cryogen
large_tank.cryogen = cryogen

# Define vertical spacing and computational grid
n_z = 100 
large_tank.z_grid = np.linspace(0, 1, n_z)

# Insulated roof
large_tank.U_roof = 0

# Time steps configuration
evap_time = 3600 * 12
large_tank.time_interval = 100
large_tank.plot_interval = evap_time / 6

# ---------------------------------------------------------
# OPTIMIZATION & DATA EXPORT
# ---------------------------------------------------------

# Load coefficients once globally to optimize performance
coeffs = load_coolprop_coeffs()
opti   = TankOptimizerJAX(large_tank)

folder = "../Results_new/Data/"
os.makedirs(folder, exist_ok=True)

LF_array_fine = jnp.linspace(0.05, 0.95, 100)
results_opt = []

# Store original LF to prevent accidental global state mutation
original_lf = opti.tank.LF

for lf_eval in LF_array_fine:
    lf_eval_float = float(lf_eval)
    
    # Update internal tank state and rebuild parameters smoothly
    opti.tank.LF = lf_eval_float
    opti.params = opti._build_params(opti.tank, coeffs)
    a = a_opt_analytical(r_U, lf_eval_float)  # Analytical initial guess for aspect ratio

    # Optimize aspect ratio
    BOR = opti._objective_bor(a, opti.params, t_final=evap_time)
    
    results_opt.append({
        "LF": lf_eval_float,
        "Optimal_AR": float(a),
        "BOR": float(BOR)
    })
    
    # Console feedback for tracking progress
    if len(results_opt) % 20 == 0:
        print(f"Processed LF {lf_eval_float:.2f}: AR = {float(a):.4f}, BOR = {float(BOR):.4e}")

# Restore original tank state
opti.tank.LF = original_lf
opti.params = opti._build_params(opti.tank, coeffs)

# Export optimal results
df_optimos = pd.DataFrame(results_opt)
filename_opts = os.path.join(folder, 'LAES_opti_12h_LFs_ru_4_opts_analytical.csv')
df_optimos.to_csv(filename_opts, index=False)
print(f"High-resolution optimals saved to: {filename_opts}")