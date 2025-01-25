# Step 1: Third-party module imports
import numpy as np # Scientific computing
import matplotlib.pyplot as plt # Visualisation
# CryoEvap module imports
from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

# Step 2: Initialise tank object
Q_roof = 0   		# Roof heat ingress / W
d_i = 8 	    	# Internal diameter / m
d_o = 8.4   		# External diameter / m
T_air = 293.15 		# Temperature of the environment K
U_L = 3.73e-3 		# Liquid overall heat transfer coefficient W/m^2/K
U_V = 3.73e-3 		# Vapour overall heat transfer coefficient W/m^2/K
Q_b = 100 		    # Heat transfer rate at the bottom / W
V_tank = 2033   	# Tank volume / m^3
LF = 0.50     		# Initial liquid filling / -
P = 101325  		# Tank operating pressure / Pa
mid_tank = Tank(d_i, d_o, V_tank, LF) # Initialize mid-scale tank
mid_tank.set_HeatTransProps(U_L, U_V, T_air, Q_b, Q_roof, eta_w= 0.8)

# Step 3: Initialise cryogen
hydrogen = Cryogen(name = "hydrogen")
hydrogen.set_coolprops(P)
mid_tank.cryogen = hydrogen

# Step 4: Simulation setup
# Calculate initial evaporation rate
print("The initial evaporation rate of " + hydrogen.name + " is %.1f kg/h" % (mid_tank.b_l_dot * 3600)) 
# Estimate transient period duration
print("Transient period = %.3f s " % mid_tank.tau)
# Minimum number of hours to achieve steady state 
tau_h = (np.floor(mid_tank.tau / 3600) + 1)

dz = 0.1 # grid spacing / m
n_z = 1 + int(np.round(mid_tank.l_V/dz, 0)) # Number of nodes
mid_tank.z_grid = np.linspace(0, 1, n_z) # Set dimensionless grid
mid_tank.U_roof = 0 # Roof overall heat transfer coefficient W/m^2/K
evap_time = 3600 * tau_h # Define evaporation time / s
mid_tank.time_interval = 60 # Time-step to record data
mid_tank.plot_interval = evap_time/6 # Interval to plot vapour temperature profiles
mid_tank.evaporate(evap_time) # Simulate the evaporation

# Step 5: Visualisation
mid_tank.plot_tv(t_unit="h") # Vapour temperature
plt.savefig("temp_vap.png", bbox_inches = 'tight')
mid_tank.plot_Q( unit="W", t_unit="h") # Heat transfer rates
plt.savefig("Q_plots.png", bbox_inches = 'tight')
mid_tank.plot_V_L(unit="L", t_unit="min") # Liquid volume
plt.savefig("V_L.png", bbox_inches = 'tight')
mid_tank.plot_BOG(unit='g/h', t_unit="min") # Boil-off gas and evaporation 
plt.savefig("BOG.png", bbox_inches = 'tight')