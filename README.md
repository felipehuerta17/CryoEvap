# CryoEvap
A Python framework for the simulation and optimization of cryogenic liquid evaporation in storage tanks.

### Requirements

* CoolProp >= 6.4.1
* NumPy >= 1.26.0
* SciPy >= 1.11.3
* Matplotlib >= 3.8.0
* Pandas >= 2.0.0
* JAX >= 0.4.26
* Diffrax >= 0.5.0
* Equinox >= 0.11.0
* Jupyter >= 1.0.0

Dependencies can be installed via pip:
```bash
pip install -r requirements.txt
```

To install CryoEvap in editable mode from the repository root:
```bash
pip install -e .
```

### Quickstart / Simplified Base Case

The following example illustrates a simplified base case to initialize a storage tank, simulate liquid nitrogen (LN2) evaporation, and perform JAX-accelerated aspect ratio optimization to minimize the boil-off rate (BOR).

```python
import numpy as np
import matplotlib.pyplot as plt
import jax

from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen
from cryoevap.optimize import TankOptimizerJAX

# Enable 64-bit precision for JAX calculations
jax.config.update("jax_enable_x64", True)

# 1. Tank Initialization (Base Case)
V_tank = 100.0         # Tank volume [m^3]
aspect_ratio = 1.0     # Height-to-diameter ratio (H/D)
d_i = ((4 * V_tank) / (np.pi * aspect_ratio))**(1/3)  # Internal diameter [m]
d_o = d_i * 1.02       # External diameter with 2% wall thickness [m]
LF = 0.5               # Initial liquid filling fraction [-]

tank = Tank(d_i, d_o, V_tank, LF)
tank.set_HeatTransProps(U_L=0.05, U_V=0.05, T_air=298.15, q_b_fixed=None, Q_roof=0, eta_w=0.90)
tank.U_b = 0.05        # Bottom heat transfer coefficient [W/m^2/K]

# 2. Cryogen Setup
nitrogen = Cryogen(name="nitrogen")
nitrogen.set_coolprops(101325)  # Operating pressure [Pa]
tank.cryogen = nitrogen

# 3. Discretization Grid & Output Intervals
tank.z_grid = np.linspace(0, 1, 25)  # Dimensionless vertical grid
tank.time_interval = 60              # Record data every 60 s

# 4. Evaporation Simulation
evap_time = 3600 * 2                 # Simulation duration [s]
tank.plot_interval = evap_time / 4
tank.evaporate(evap_time)

# Visualise results
tank.plot_BOG(unit="kg/h", t_unit="h")
plt.show()

# 5. JAX-Accelerated Optimization
# Find the optimal aspect ratio minimizing the Boil-Off Rate (BOR)
opti = TankOptimizerJAX(tank)
opt_ar, opt_tar, min_bor = opti.optimize(
    t_final=evap_time,
    coarse_samples=20,
    fine_samples=50,
    ar_min=0.2,
    ar_max=3.0
)
print(f"Optimal Aspect Ratio: {opt_ar:.4f} | Min BOR: {min_bor:.4e}")
```

Further examples, parameter sensitivity analyses, and case studies can be found in the `/notebooks` and `/Opti_cases` directories.
