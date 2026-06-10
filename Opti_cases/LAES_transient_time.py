# Ensure that python finds the submodules
import sys
sys.path.append("..") # Adds higher directory to python modules path.

# Scientific computing
import numpy as np

#import optimization class
from cryoevap.optimize import Opti_jax

# Visualisation
import matplotlib.pyplot as plt

## Module imports
# Import the storage tank Class
from cryoevap.storage_tanks import Tank

# Import Cryogen class
from cryoevap.cryogens import Cryogen

# Import JAX library
import jax
import jax.numpy as jnp

# Import pandas for data handling
import pandas as pd

import itertools

# Set JAX to use 64-bit precision
jax.config.update("jax_enable_x64", True)

# LNG tank properties
Q_roof = 0              # Roof heat ingress / W
T_air  = 18.08 + 273.15  # Temperature of the environment K

V_tank = 4831 # Tank volume / m^3

a   = 0.5
d_i = ((4 * V_tank)/(np.pi * a))**(1/3) # internal diameter / m

# Thickness of the in % of the internal diameter
ST  = 1.02
d_o = d_i * ST # external diameter / m

# Set overall heat transfer coefficient through the walls for liquid and vapour
U_L = 8.72e-2 # W/m2/K
U_V = 8.72e-2 # W/m2/K

# Set overall bottom heat transfer coefficient respect the r_U
r_U = 0.25
U_b = r_U * U_L # W/m2/K

# Initial liquid filling / Dimensionless
LF = 0.95

# Wall heat fraction
eta = 0.9
# 
# Specify tank operating pressure
P = 101325*3 # Pa

# # Initialise cryogen
cryogen = Cryogen(name = "nitrogen")
cryogen.set_coolprops(P)

def estimate_transient_period(a, r_U, LF):
    eta_w = 0.9

    d_i = ((4 * V_tank)/(np.pi * a))**(1/3) # internal diameter / m
    ST  = 1.02
    d_o = d_i * ST # external diameter / m

    # Initialize large-scale tank
    large_tank = Tank(d_i, d_o, V_tank, LF)
    large_tank.set_HeatTransProps(U_L, U_V, T_air, q_b_fixed = None, Q_roof = 0, eta_w = 0.90)
    large_tank.cryogen = cryogen
    large_tank.U_b = U_L*r_U

    # Minimum number of hours to achieve steady state 
    return large_tank.tau

a_vals   = np.linspace(0.2, 2.5, 500)
r_U_vals = np.linspace(0.25, 4.0, 50) 
LF_vals  = [0.05, 0.95]

results = {'a': a_vals}

for r_U, LF in itertools.product(r_U_vals, LF_vals):
    col_name = f'tau_ru_{r_U:.2f}_LF_{LF:.2f}'
    results[col_name] = estimate_transient_period(a_vals, r_U, LF)

# 4. Guardar en CSV
df = pd.DataFrame(results)
df.to_csv('../Results/Data/transient_period.csv', index=False)
print("CSV guardado exitosamente.")
