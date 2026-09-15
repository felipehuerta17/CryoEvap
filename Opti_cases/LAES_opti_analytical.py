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
n_z = 50 
large_tank.z_grid = np.linspace(0, 1, n_z)

# Insulated roof
large_tank.U_roof = 0

# Time steps configuration
evap_time = 3600 * 12
large_tank.time_interval = 1200
large_tank.plot_interval = evap_time / 6

# ---------------------------------------------------------
# OPTIMIZATION & DATA EXPORT
# ---------------------------------------------------------

# Load coefficients once globally to optimize performance
coeffs = load_coolprop_coeffs()
opti   = TankOptimizerJAX(large_tank)


def analytical_fixed_point(opti_params, r_U, evap_time, a_guess=1.0, tol=1e-8, max_iter=1000):

    V_tank = float(opti_params['V'])
    LF_0   = float(opti_params['LF'])
    
    # Extraemos el factor geométrico constante (ST = 1.02)
    xi     = float(opti_params['d_o'] / opti_params['d_i'])
    eta_w  = float(opti_params['eta_w'])
    
    T_L    = float(opti_params['T_sat'])
    T_air  = float(opti_params['T_air'])
    U_V    = float(opti_params['U_V'])
    U_L    = float(opti_params['U_L'])
    h_L    = float(opti_params['h_L'])
    h_V    = float(opti_params['h_V'])
    rho_L  = float(opti_params['rho_L'])
    
    k_V_poly   = opti_params['k_V_poly']
    cp_V_poly  = opti_params['cp_V_poly']
    rho_V_poly = opti_params['rho_V_poly']
    
    delta_T_sat = T_air - T_L
    a_current = float(a_guess)
    delta_T_V_bar = delta_T_sat / 2.0  
    
    # Picard iterations
    for _ in range(max_iter):
        d_i = ((4.0 * V_tank) / (jnp.pi * a_current))**(1/3)
        d_o = d_i * xi
        
        A_T = jnp.pi * d_i**2 / 4.0
        l_tank = V_tank / A_T
        l_V = l_tank * (1.0 - LF_0)
        
        A_L = jnp.pi * d_o * l_tank * LF_0
        A_V = jnp.pi * d_o * l_tank * (1.0 - LF_0)
        
        T_avg = T_air - delta_T_V_bar
        kV    = float(jnp.polyval(k_V_poly, T_avg))
        cp_V  = float(jnp.polyval(cp_V_poly, T_avg))
        rho_V = float(jnp.polyval(rho_V_poly, T_avg))
        
        Q_L  = U_L * A_L * delta_T_sat
        Q_b  = r_U * U_L * A_T * delta_T_sat 
        Q_wi = U_V * A_V * eta_w * delta_T_V_bar
        
        BL_0 = (Q_L + Q_b + Q_wi) / (h_V - h_L)
        v_z0 = (4.0 * BL_0) / (rho_V * jnp.pi * d_i**2)
        
        H = rho_V * cp_V * v_z0
        S = (4.0 * U_V * d_o / d_i**2) * (1.0 - eta_w)
        
        delta_H = jnp.sqrt(H**2 + 4.0 * kV * S)
        chi_plus = (H + delta_H) / (2.0 * kV)
        chi_minus = (H - delta_H) / (2.0 * kV)
        
        R = (chi_minus / chi_plus) * jnp.exp(-l_V * delta_H / kV)
        I1 = (1.0 / (1.0 - R)) * (jnp.exp(l_V * chi_minus) - 1.0) / chi_minus
        I2 = -chi_minus * jnp.exp(l_V * chi_minus) / (chi_plus**2 * (1.0 - R)) + R / (chi_plus * (1.0 - R))
        
        delta_T_V_bar_new = float(- (I1 + I2) * (T_L - T_air) / l_V)
        delta_T_V_bar_new = float(jnp.clip(delta_T_V_bar_new, 0.0, delta_T_sat))
        
        term_denominator = LF_0 + (U_V / U_L) * (1.0 - LF_0) * eta_w * (delta_T_V_bar_new / delta_T_sat)

        # Optimal aspect ratio based on the current delta_T_V_bar
        a_new = r_U / (2.0 * xi * term_denominator)
        
        if abs(a_new - a_current) < tol:
            a_current = a_new
            delta_T_V_bar = delta_T_V_bar_new
            break
            
        a_current = a_new
        delta_T_V_bar = delta_T_V_bar_new

    # Analytical BOR
    h_LV = h_V - h_L
    Q_slope = (4.0 * d_o / d_i**2) * (U_L * delta_T_sat - eta_w * U_V * delta_T_V_bar)
    C_wneq = - Q_slope / (rho_L * h_LV)
    
    Q_b_total = r_U * U_L * A_T * delta_T_sat
    Q_wall_total = eta_w * U_V * (jnp.pi * d_o * l_tank) * delta_T_V_bar
    D_wneq = - (Q_b_total + Q_wall_total) / (rho_L * h_LV)
    
    V_L_0 = LF_0 * V_tank
    V_L_tf = (D_wneq / C_wneq) * (jnp.exp(C_wneq * evap_time) - 1.0) + V_L_0 * jnp.exp(C_wneq * evap_time)
    bor_exact = ((V_L_0 - V_L_tf) / V_L_0) * (86400.0 / evap_time)
    
    return float(a_current), float(bor_exact)

LF_array_fine = jnp.linspace(0.05, 0.95, 100)
for r_U in [0.25, 4]:
    results_opt = []

    for lf_eval in LF_array_fine:

        lf_eval_float = float(lf_eval)
        opti.tank.LF  = lf_eval_float
        opti.params   = opti._build_params(opti.tank, coeffs)

        a_optimo, bor_minimo = analytical_fixed_point(
            opti_params=opti.params,
            r_U=r_U,
            evap_time=evap_time,
            a_guess=0.5 
        )

        results_opt.append({
            "LF": lf_eval_float,
            "Optimal_AR": a_optimo,
            "BOR": bor_minimo
        }) 
    
    pd.DataFrame(results_opt).to_csv(f"../Results_new/Data/analytical_optimal_r_U_{r_U}.csv", index=False)
    print(f"Results for r_U = {r_U} saved to ../Results_new/Data/analytical_optimal_r_U_{r_U}.csv")