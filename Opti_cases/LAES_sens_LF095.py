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
# TANK & CRYOGEN SETUP
# ---------------------------------------------------------
Q_roof = 0               # Roof heat ingress / W
T_air  = 18.08 + 273.15  # Temperature of the environment / K

V_tank = 4831.0          # Tank volume / m^3

a   = 0.5
d_i = ((4 * V_tank) / (np.pi * a))**(1/3) # Internal diameter / m

# Thickness in % of the internal diameter
ST  = 1.02
d_o = d_i * ST           # External diameter / m

# Base wall heat transfer coefficients
U_L = 8.72e-2            # W/m2/K
U_V = 8.72e-2            # W/m2/K

# Initial liquid filling / Dimensionless
LF = 0.95

# Specify tank operating pressure
P = 101325 * 3           # Pa

# Initialise cryogen
cryogen = Cryogen(name="nitrogen")
cryogen.set_coolprops(P)

# Initialize base large-scale tank object
large_tank = Tank(d_i, d_o, V_tank, LF)
large_tank.set_HeatTransProps(U_L, U_V, T_air, q_b_fixed=None, Q_roof=0, eta_w=0.90)

# Set cryogen
large_tank.cryogen = cryogen

# Define vertical spacing and computational grid
dz = 0.1
n_z = 1 + int(np.round(large_tank.l_V / dz, 0))
large_tank.z_grid = np.linspace(0, 1, n_z)

# Insulated roof
large_tank.U_roof = 0

# Define simulation time properties
evap_time = 3600 * 12
large_tank.time_interval = 1200.0

# Define execution configurations for both r_U values
r_U_scenarios = [4.0, 0.25]
folder_data = "../Results_new/Data/"
os.makedirs(folder_data, exist_ok=True)

# Load polynomial thermophysical coefficients once globally
coeffs = load_coolprop_coeffs()
opti = TankOptimizerJAX(large_tank)

# Define JAX auto-diff functions for design elasticity evaluation
grad_fun = jax.grad(opti._objective_bor, argnums=0, allow_int=True) 
hess_fun = jax.jacfwd(grad_fun, argnums=0)
mixed_grad_fun = jax.jacfwd(grad_fun, argnums=1) 

params_to_analyze = ['U_L', 'U_V', 'U_b', 'eta_w']

# ---------------------------------------------------------
# MULTI-SCENARIO EXECUTION LOOP
# ---------------------------------------------------------
for r_U in r_U_scenarios:
    # Resolve file name convention string based on numerical value
    ru_str = "4" if r_U == 4.0 else "025"
    print(f"\n=========================================================")
    print(f"RUNNING ANALYSIS FOR: r_U = {r_U} (Output suffix: ru_{ru_str})")
    print(f"=========================================================")

    # Update specific tank instance bottom properties
    opti.tank.U_b = U_L * r_U
    opti.params = opti._build_params(opti.tank, coeffs)

    # Load corresponding optimization file
    opts_file = os.path.join(folder_data, f'LAES_opti_12h_LFs_ru_{ru_str}_opts.csv')
    if not os.path.exists(opts_file):
        raise FileNotFoundError(f"Required baseline optimization file missing: {opts_file}")

    df_opts = pd.read_csv(opts_file)
    opt_data = df_opts[np.isclose(df_opts['LF'], LF, atol=1e-3)].iloc[0]
    optimal_aspect_ratio = float(opt_data['Optimal_Geometric_AR'])
    
    # Calculate baseline BOR at local optimum
    bor_base = float(opti._objective_bor(optimal_aspect_ratio, opti.params, evap_time))
    print(f"Optimal AR loaded: {optimal_aspect_ratio:.6f} | Base BOR: {bor_base:.6f} %/day")

    # Compute Local Sensitivities
    local_sensitivities = {}
    for param in params_to_analyze:
        val_param = float(opti.params[param])
        grad = float(opti.calculate_sensibility_param(param, optimal_aspect_ratio, evap_time))
        local_sensitivities[param] = grad * (val_param / bor_base)

    # Compute Design Elasticities (Optimum Sensitivity via Implicit Function Theorem)
    H_xx = float(hess_fun(optimal_aspect_ratio, opti.params, evap_time))
    grads_xp = mixed_grad_fun(optimal_aspect_ratio, opti.params, evap_time)
    
    design_elasticities = {}
    for param in params_to_analyze:
        d2f_dxdp = float(grads_xp[param])
        dx_dp = - (d2f_dxdp / H_xx)
        val_param = float(opti.params[param])
        design_elasticities[param] = dx_dp * (val_param / optimal_aspect_ratio)

    # Export specific execution scenario results to independent CSV
    lf_str = int(LF * 100)
    export_filename = os.path.join(folder_data, f"Sensitivity_LF{lf_str:02d}_ru_{ru_str}.csv")
    
    df_export = pd.DataFrame({
        'Param': params_to_analyze,
        'Local_sensitivity': [local_sensitivities[p] for p in params_to_analyze],
        'Design_elasticity': [design_elasticities[p] for p in params_to_analyze]
    })
    
    df_export.to_csv(export_filename, index=False)
    print(f"Results successfully exported to: {export_filename}")