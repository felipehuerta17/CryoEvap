# Ensure that python finds the submodules
import sys
sys.path.append("..") # Adds higher directory to python modules path.
import numpy as np # Scientific computing
import matplotlib.pyplot as plt # Visualisation
from scipy.integrate import simpson

# Import the storage tank and cryogen classes
from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen
from scipy.optimize import minimize, dual_annealing

# Simulation set up
# LNG tank properties

d_i = 2.106         # Internal diameter / m
d_o = 2.106+2.08e-3 # External diameter / m
V_tank = 4.89       # Spherical tank volume / m^3
LF = 0.95           # Initial liquid filling / -
P = 117e3           # Pa
Q_b = 0 # W, 
Q_roof = 0          # Roof heat ingress / W
T_air = 350         # Temperature of the environment K

# Set overall heat transfer coefficient through the walls for liquid and vapour
U_L = 0.01063 # W/m2/K
U_V = 0.01 # W/m2/K

# Vapour and liquid geometry
Geo_v1 = "spherical"
Geo_l1 = "spherical"

# Initialize large-scale tank
small_tank1 = Tank(d_i, d_o, V_tank, Geo_v1, Geo_l1,LF)
small_tank1.set_HeatTransProps(U_L, U_V, T_air, Q_roof, Q_b, eta_w = 0)

# Initialise and set cryogen
hydrogen1 = Cryogen(name = "hydrogen")
hydrogen1.set_coolprops(P)
small_tank1.cryogen = hydrogen1

# Define vertical spacing
dz = 0.0075/2 # n_z = 78 for this configuration

# Calculate number of nodes
n_z = 1 + int(np.round(small_tank1.l_V/dz, 0))

# Define dimensionless computational grid
small_tank1.z_grid = np.linspace(0, 1, n_z)

# Insulated roof
small_tank1.U_roof = 0

evap_time = 3600 * 4 # Simulation time / s

# Time step to record data, relevant for plotting integrated quantities
# such as the vapour to liquid heat transfer rate, Q_VL
small_tank1.time_interval = 60

# Time step to plot each vapour temperature profile
small_tank1.plot_interval = evap_time/6

def fobj(U_V):
    # U_L, U_V = p
    small_tank1 = Tank(d_i, d_o, V_tank, Geo_v1, Geo_l1, LF)
    small_tank1.set_HeatTransProps(U_L, U_V, T_air, Q_roof, Q_b, eta_w = 0)
    small_tank1.cryogen = hydrogen1
    small_tank1.z_grid = np.linspace(0, 1, n_z)
    small_tank1.U_roof = 0 # Insulated roof
    small_tank1.time_interval = 60 # Save timesteps every 60 s
    small_tank1.plot_interval = evap_time/6 # Vapour temperature plot intervals
    small_tank1.evaporate(evap_time)

    # Calculate objective functions
    T_V_end = small_tank1.sol.y[1:, -1] # Get the vapour temperature at the last time-step
    # Compute vertical vapour temperature gradient at the end of the simulation
    Tvz_end = (T_V_end[-1] - T_V_end[0])/(small_tank1.l_V * 100)
    print("dTv_dz|t=end = %.3f K/cm" % Tvz_end)
    J2 = (Tvz_end- 0.636)**2
    print("T_roof = %.3f K" % T_V_end[-1])
    print("T_sat = %.3f K" % T_V_end[0])
    print("J2 = %.3f" % J2)
    print("U_V = %.3f W m^-2 K^-1" % U_V)
    return np.sqrt(J2)

#U_V_list = []
#fobj_list = []
#for U_V in np.linspace(0.056-5e-4, 0.056+5e-4, 20):
#    J = fobj(U_V)
#    U_V_list.append(U_V)
#    fobj_list.append(J)

#U_V_list = np.array(U_V_list)
#fobj_list = np.array(fobj_list)
#plt.plot(U_V_list, fobj_list, marker='o')
#plt.savefig("fobj_vs_UV_H2.png", bbox_inches='tight')
#  There may be multiple minima so we will try with DA
#results = dict()
#results['DA'] = dual_annealing(fobj, bounds=[(0.1, 0.4)], maxiter=20)
res = minimize(fobj, [U_V], bounds=[(0.0556, 0.056)], method='L-BFGS-B')
print("Optimized U_V = %.5f W/m2/K" % res.x[0])

