import sys
import os
import time
sys.path.append("..")

import numpy as np
import jax
import jax.numpy as jnp
import pandas as pd

from cryoevap.optimize import TankOptimizerJAX, load_coolprop_coeffs
from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

jax.config.update("jax_enable_x64", True)

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

# Wall heat fraction
eta = 0.9

# Specify tank operating pressure
P = 101325 * 3           # Pa

# Initialise cryogen
cryogen = Cryogen(name="nitrogen")
cryogen.set_coolprops(P)

# Initialize large-scale tank
LF_initial = 0.05
large_tank = Tank(d_i, d_o, V_tank, LF_initial)
large_tank.set_HeatTransProps(U_L, U_V, T_air, q_b_fixed=None, Q_roof=0, eta_w=0.90)
large_tank.U_b = U_L * r_U

# Set cryogen
large_tank.cryogen = cryogen

# Define vertical spacing and computational grid
n_z = 50 
large_tank.z_grid = np.linspace(0, 1, n_z)

# Insulated roof
large_tank.U_roof = 0

# Time steps configuration
evap_time = 3600 * 24 * 7
large_tank.time_interval = 1200
large_tank.plot_interval = evap_time / 6

# ---------------------------------------------------------
# PERFORMANCE TESTING
# ---------------------------------------------------------

coeffs = load_coolprop_coeffs()
opti   = TankOptimizerJAX(large_tank)

test_LFs = [0.05, 0.95]
num_runs = 3

print(f"Performance Test for r_U = {r_U}")
print("-" * 40)

# warm-up to exclude JAX JIT compilation time from the average
print("Running warm-up (JIT compilation)...")
opti.tank.LF = float(0.05)
opti.params = opti._build_params(opti.tank, coeffs)
_ = opti.optimize(t_final=evap_time, coarse_samples=100, fine_samples=500, ar_min=0.05, ar_max=3.0)
print("Warm-up complete.\n")

results = []

for lf in test_LFs:
    opti.tank.LF = float(lf)
    opti.params = opti._build_params(opti.tank, coeffs)
    
    times = []
    print(f"Testing LF = {lf}")
    for i in range(num_runs):
        start_time = time.time()
        
        opt_ar, opt_tar, min_bor = opti.optimize(
            t_final=evap_time,
            coarse_samples=100,
            fine_samples=500,
            ar_min=0.05,
            ar_max=3.0
        )
        
        # Ensure asynchronous JAX dispatch completes by fetching a value
        _ = float(opt_ar)
        
        end_time = time.time()
        elapsed = end_time - start_time
        times.append(elapsed)
        results.append({
            "LF": lf,
            "run": i + 1,
            "cpu_time_s": elapsed,
            "optimal_value": opt_ar,
        })
        print(f"  Run {i+1}: {elapsed:.4f} seconds")
        
    avg_time = sum(times) / num_runs
    print(f"  Average CPU time for LF={lf}: {avg_time:.4f} seconds\n")

df_results = pd.DataFrame(results)
folder     = "../Results_new/Data/"
file_name  = "cpu_time_JAX.csv"
df_results.to_csv(folder+file_name, index=False)
