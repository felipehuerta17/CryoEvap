# CryoEvap
A Python software for the simulation of the evaporation of cryogenic liquids in storage tanks

### Requirements

* CoolProp >= 6.4.1
* SciPy >= 1.11.3
* Matplotlib >= 3.8.0
* Jupyter >= 1.0.0

The software requires Python 3 and it has been tested on Windows and Linux. The most straightforward way to install Python 3 is through Anaconda. The required modules are Matplotlib, NumPy, SciPy and CoolProp. They can be installed using a package manager such as conda or pip. In a system terminal with pip and python added to the path, the packages can be installed typing pip install matplotlib numpy scipy coolprop --upgrade. To take advantage of CryoEvap interactivity, it is suggested to use Jupyter Notebook as the integrated development environment (IDE). 

Once all the prerequisites are satisfied, the source code can be obtained from the GitHub repository by cloning the repository branch with the command 

`git clone https://github.com/felipehuerta17/CryoEvap.git`

To install CryoEvap, open Anaconda Prompt (on Windows) or Terminal (on Linux/Mac) and navigate to the path where the CryoEvap folder is located. If the repository was cloned, type cd /path/CryoEvap. If the compressed file was downloaded, the directory will be /path/CryoEvap-1.0.0. From this folder, install the package typing in the terminal 

`pip install .`

A successful installation will show in the terminal the message “Successfully installed cryoevap”.

### How to use CryoEvap

The workflow of CryoEvap is divided in five steps: module import, tank initialisation, cryogen initialisation, simulation setup and visualisation. The following code snippet illustrates the minimum code necessary to simulate the evaporation of LN2 in a lab-scale storage tank, and it can be used as a template for any scenario.  In step 1, NumPy and Matplotlib are imported to directly visualise the results in the notebook. Additionally, the CryoEvap Tank and Cryogen classes are imported, as they act as an interface for all simulation functionalities. In step 2, the geometrical and heat transfer properties of the tank are defined, as well as the operating pressure and initial liquid filling. In step 3, the cryogen is first initialised with its name and then their properties at the operating pressure are set. In step 4, the grid spacing is set to calculate the number of nodes in the vertical direction, z, and then this value is set as a tank property. The overall roof heat transfer coefficient is also set to illustrate that tank properties can be modified after the tank is constructed. The time_interval property establishes the time-step at which simulation data will be recorded for post-processing. The plot_interval property defines the interval at which vapour temperatures will be plotted. Step 4 ends with the function evaporate, which receives the simulation time and performs the simulation. Finally, Step 5 illustrates the syntax to produce the figures that summarise the results. The parameters t_unit and unit allow the user to control the units of time and the corresponding dependent variables to improve visualisations. Further examples that set up the scenarios and perform the simulations can be found in the Jupyter Notebooks located in the /notebooks folder of CryoEvap GitHub repository.

```python
# Step 1: Third-party module imports
import numpy as np # Scientific computing
import matplotlib.pyplot as plt # Visualisation
# CryoEvap module imports
from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

# Step 2: Initialise tank object
Q_roof = 0 		# Roof heat ingress / W
d_i = 0.201 		# Internal diameter / m
d_o = 0.204   		# External diameter / m
T_air = 298.15 		# Temperature of the environment K
U_L = 0.026 		# Liquid overall heat transfer coefficient W/m^2/K
U_V = 0.026 		# Vapour overall heat transfer coefficient W/m^2/K
Q_b = 0 		# Heat transfer rate at the bottom / W
V_tank = 6.75e-3 	# Tank volume / m^3
LF = 0.278  		# Initial liquid filling / -
P = 100000 		# Tank operating pressure / Pa
small_tank = Tank(d_i, d_o, V_tank, LF) # Initialize large-scale tank
small_tank.set_HeatTransProps(U_L, U_V, T_air, Q_roof, Q_b, eta_w = 0.963)

# Step 3: Initialise cryogen
nitrogen = Cryogen(name = "nitrogen")
nitrogen.set_coolprops(P)
small_tank.cryogen = nitrogen	# Set initialised cryogen as a tank property

# Step 4: Simulation setup
dz = 0.005 # grid spacing / m
n_z = 1 + int(np.round(small_tank.l_V/dz, 0)) # Number of nodes
small_tank.z_grid = np.linspace(0, 1, n_z) # Set dimensionless grid
small_tank.U_roof = 0 # Roof overall heat transfer coefficient W/m^2/K
evap_time = 3600*16 # Define evaporation time / s
small_tank.time_interval = 60 # Time-step to record data
small_tank.plot_interval = evap_time/6 # Interval to plot vapour temperature profiles
small_tank.evaporate(evap_time) # Simulate the evaporation

# Step 5: Visualisation
small_tank.plot_tv(t_unit="min") # Vapour temperature
small_tank.plot_Q( unit="W", t_unit="min") # Heat transfer rates
small_tank.plot_V_L(unit="L", t_unit="min") # Liquid volume
small_tank.plot_BOG(unit='g/h', t_unit="min") # Boil-off gas and evaporation rates
```
