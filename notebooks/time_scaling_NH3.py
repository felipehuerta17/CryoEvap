# Ensure that python finds the submodules
import sys
sys.path.append("..") # Adds higher directory to python modules path.

# Scientific computing
import numpy as np

# Visualisation
import matplotlib.pyplot as plt

## Module imports
# Import the storage tank Class
from cryoevap.storage_tanks import Tank 

# Import Cryogen class
from cryoevap.cryogens import Cryogen

# Import Pandas for data handling
import pandas as pd

# Import time for tracking simulation time
import time

# Dimensions
# Wall thickness 
e = 10.6972*2 - 10.24*2 #m

# Vertically orientated cylindrical tank volume
V_tank = 200000.0 #m^3

# Aspect ratio
a = 0.5

# Internal diameter
d_i  = ((4 * V_tank)/(np.pi * a))**(1/3) # internal diameter / m
d_o  = d_i + e # external diameter / m

# Initial liquid filling / Dimensionless
LF = 0.95

# LNG tank properties
Q_roof = 0           # Roof heat ingress / W
T_air  = 5.3286+273.15 # Temperature of the environment K

# Set overall heat transfer coefficient through the walls for liquid and vapour
U_L = 8.86344e-02   # W/m2/K
U_V = 8.86344e-02   # W/m2/K
U_b = 8.80129e-2      # W/m2/K
h_L = 135.08  # W/m2/K

# Set wall tank properties
Init_T = True   # True: T_air, False: T_L
k_w    = 0.0411 # W/mK
cp_w   = 900    # J/kg K
rho_w  = 60     # kg/m3

# Specify heat transfer rate at the bottom to prevent ground heating
Q_b = None # W, 

# Specify tank operating pressure
P = 101325 # Pa

T_env_avg = 5.3286+273.15     # K
T_range   = 15                # K
h_env     = 14.839               # W/m2K   

p_anual   = np.array([-0.00022, 0.08395, 0.989792])

evap_time = 3600 * 24 * 60 # 2 months in seconds 

def run_simulation(n_nodes_r, n_nodes_z):
    # Initialise cryogen
    ammonia = Cryogen(name = "ammonia")
    ammonia.set_coolprops(P)

    # Initialize large-scale tank
    large_tank = Tank(d_i, d_o, V_tank, LF)

    # Set heat transfer properties
    large_tank.set_HeatTransProps(U_L, U_V, T_air, Q_b_fixed = None, Q_roof = 0, eta_w = 0.70, k_w = k_w, rho_w = rho_w, cp_w = cp_w, h_L = h_L, T_init = Init_T)
   
    # Bottom heat ingress
    large_tank.U_b = U_b

    # Insulated roof
    large_tank.U_roof = 0

    # Set cryogen
    large_tank.cryogen = ammonia

    # Set environmental props
    large_tank.set_EnvironmentalProps(T_avg_day = T_env_avg, T_range_day = T_range, h_env = h_env,
                                  p_anual=None)

    # Define dimensionless computational grid
    large_tank.z_grid = np.linspace(0, 1, n_nodes_z)
    large_tank.r_grid = np.linspace(0, 1, n_nodes_r)

    # Time step to record data,
    large_tank.time_interval = 3600 * 3

    # Time step to plot each vapour temperature profile
    large_tank.plot_interval = evap_time/6

    # Simulate the evaporation
    init_time = time.time()
    large_tank.evaporate(evap_time)
    final_time = time.time()

    return final_time - init_time

# Run the simulation for different grid sizes
r_grid = np.array([4, 8, 16, 32, 64, 128])
z_grid = np.array([4, 8, 16, 32, 64, 128])
simulation_times = []
desviation_times = []

for n_r, n_z in zip(r_grid, z_grid):
    print(f"Running simulation for grid size {n_r}x{n_z}...")
    sim_time = []
    for i in range(3): # Run each simulation 3 times to get an average time
        sim_time.append(run_simulation(n_r, n_z))
    mean_time  = np.mean(sim_time)
    desv_time  = np.std(sim_time)
    simulation_times.append(mean_time)
    desviation_times.append(desv_time)
    print(f"Grid size {n_r}x{n_z}: {mean_time:.2f} seconds")

# Export results to a DataFrame
results_df = pd.DataFrame({
    'r_grid': r_grid,
    'z_grid': z_grid,
    'mean_time': simulation_times,
    'std_time': desviation_times
})
# Save results to a CSV file
results_df.to_csv('simulation_times_big_o.csv', index=False)