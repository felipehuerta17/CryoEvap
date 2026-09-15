import sys
import os
import time
sys.path.append("..")

import numpy as np
import pandas as pd

# Import the old scipy version instead of JAX
from cryoevap.optimize.opti import Opti
from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

# ---------------------------------------------------------
# TANK & CRYOGEN SETUP
# ---------------------------------------------------------
def create_tank(LF_initial):
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
    return large_tank



# ---------------------------------------------------------
# PERFORMANCE TESTING (SCIPY)
# ---------------------------------------------------------
evap_time_hours = 24 * 7
test_LFs = [0.05, 0.95]
num_runs = 3
results  = []

for lf in test_LFs:

    
    times = []
    print(f"Testing LF = {lf} with SciPy")
    for i in range(num_runs):
        
        # Update LF on the tank
        large_tank = create_tank(lf)
        large_tank.time_interval = 1200
        large_tank.plot_interval = evap_time_hours * 3600 / 6

        # Instantiate the old Opti class for the current tank setup
        opti_scipy = Opti(
            Tank_obj=large_tank, 
            time=evap_time_hours, 
            dz=0.1, 
            thickness=0.02, 
            x0=1.0, 
            bounds=[0.05, 3.0],
        )

        start_time = time.time()
        
        # Run optimization using the SciPy implementation
        opt_ar = opti_scipy.aspect(verbose=3, tol=1e-12)
        
        end_time = time.time()
        elapsed = end_time - start_time
        times.append(elapsed)
        results.append({
            "LF": lf,
            "run": i + 1,
            "cpu_time_s": elapsed,
            "optimal_value": opt_ar,
        })
        print(f"  Run {i+1}: {elapsed:.4f} seconds (Opt AR: {opt_ar:.4f})")
        
    avg_time = sum(times) / num_runs
    print(f"  Average CPU time for LF={lf}: {avg_time:.4f} seconds\n")

df_results = pd.DataFrame(results)
folder     = "../Results_new/Data/"
file_name  = "cpu_time_SCIPY.csv"
df_results.to_csv(os.path.join(folder, file_name), index=False)
print(f"Saved results to {os.path.join(folder, file_name)}")
