import jax
import jax.numpy as jnp
from jax.example_libraries import optimizers
from jax import grad
from diffrax import diffeqsolve, ODETerm, Tsit5, SaveAt, PIDController, BacksolveAdjoint
import pandas as pd
import matplotlib.pyplot as plt
# Import polyfit function from CoolProp
folder   = '../Cryoevap/cryogens/Coeffs/'
cp_V_df  = pd.read_csv(folder + 'coeffs_cpV.csv')
k_V_df   = pd.read_csv(folder + 'coeffs_kV.csv')
rho_V_df = pd.read_csv(folder + 'coeffs_rhoV.csv')

class Opti_jax:
    def __init__(self, tank_obj):
        self.tank = tank_obj

        # Crear dict con parámetros JAX para usar en funciones
        self.params = self.make_params(self, self.tank)

    @staticmethod
    def make_params(self, tank):
        cryo = tank.cryogen
        return {
            "V":      tank.V,
            "d_i":    tank.d_i,
            "d_o":    tank.d_o,
            "LF":     tank.LF,
            "eta_w":  tank.eta_w,
            "T_air":  tank.T_air,
            "U_L":    tank.U_L,
            "U_V":    tank.U_V,
            "dz":     tank.z_grid[1] - tank.z_grid[0],
            "z_grid": tank.z_grid,
            "T_sat":  cryo.T_sat,
            "h_L":    cryo.h_L,
            "h_V":    cryo.h_V,
            "rho_L":  cryo.rho_L,
            'time_interval': tank.time_interval,
            "cp_V_poly": jnp.array(cp_V_df[cryo.name].values, dtype=jnp.float64),
            "k_V_poly": jnp.array(k_V_df[cryo.name].values, dtype=jnp.float64),
            "rho_V_poly": jnp.array(rho_V_df[cryo.name].values, dtype=jnp.float64)
        }
    # Funciones internas adaptadas de tu código
    @staticmethod
    @jax.jit
    def cp_V_fun(T, p):
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def k_V_fun(T, p):
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def rho_V_fun(T, p):
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def sys_isobaric_jax(t, y, args):
        aspect_ratio, p = args
        d_i = ((4 * p["V"])/(jnp.pi * aspect_ratio))**(1/3)
        d_o = d_i + 0.02 

        V_L = y[0]
        T_V = y[1:]

        A_T   = jnp.pi * d_i**2 / 4
        l     = p["V"] / A_T
        l_dry = l * (1 - p['LF'])
        dz    = p["dz"] * l_dry

        k_V   = jnp.mean(Opti_jax.k_V_fun(T_V, p["k_V_poly"]))
        cp_V  = jnp.mean(Opti_jax.cp_V_fun(T_V, p["cp_V_poly"]))
        rho_V = jnp.mean(Opti_jax.rho_V_fun(T_V, p["rho_V_poly"]))

        A_V = jnp.pi * d_o * l * (1 - p['LF'])
        A_L = jnp.pi * d_o * l * p['LF']

        Q_Lin = p["U_L"] * A_L * (p["T_air"] - p["T_sat"])
        Q_VL  = k_V * A_T * (-3 * T_V[0] + 4 * T_V[1] - T_V[2]) / (2 * dz)
        Q_b   = p["U_L"] * A_T * (p["T_air"] - p["T_sat"])
        Q_wi  = p["U_V"] * A_V * p["eta_w"] * (p["T_air"] - jnp.mean(T_V))

        BL_0  = (Q_Lin + Q_b + Q_wi) / (p["h_V"] - p["h_L"])
        v_z   = 4 * BL_0 / (rho_V * jnp.pi * d_i ** 2)
        v_int = v_z * (rho_V / p["rho_L"])

        alpha = k_V / (rho_V * cp_V)

        dT      = jnp.zeros_like(T_V)
        dT_dz   = (T_V[1:-1] - T_V[:-2]) / dz
        d2T_dz2 = (T_V[:-2] - 2 * T_V[1:-1] + T_V[2:]) / (dz ** 2)

        S_wall = (4 * p["U_V"] * d_o / (d_i ** 2)) * (p["T_air"] - T_V[1:-1]) * (1 - p["eta_w"])

        dT = dT.at[1:-1].set(alpha * d2T_dz2 - (v_z - v_int) * dT_dz + (alpha / k_V) * S_wall)

        dT = dT.at[0].set(0.0)
        dT = dT.at[-1].set((4 * dT[-2] - dT[-3]) / 3)

        dV = (-1 / p["rho_L"]) * (Q_Lin + Q_b + Q_wi + Q_VL) / (p["h_V"] - p["h_L"])

        return jnp.concatenate([jnp.array([dV], dtype=jnp.float64), dT])

    @staticmethod
    @jax.jit
    def BOR(V_L, t):
        return (1.0 - V_L[-1] / V_L[0]) * (86400.0 / t[-1])

    def evaporate(self, aspect_ratio):
        VL_0 = jnp.array(self.params['V'] * self.params['LF'], dtype=jnp.float64)
        Tv_0 = jnp.ones(len(self.params['z_grid']), dtype=jnp.float64) * self.params['T_sat']
        IC = jnp.concatenate([jnp.array([VL_0], dtype=jnp.float64), Tv_0])

        term = ODETerm(self.sys_isobaric_jax)

        sol = diffeqsolve(
            term,
            solver=Tsit5(),
            t0=0,
            t1=self.time,
            dt0=0.01,
            y0=IC,
            args=(aspect_ratio, self.params),
            saveat=SaveAt(ts=jnp.arange(0, self.time + 1, self.params['time_interval'])),
            max_steps=1000000,
            stepsize_controller=PIDController(rtol=1e-8, atol=1e-8),
            adjoint=BacksolveAdjoint()
        )
        return sol

    def objective_function(self, aspect_ratio):
        a_eff = 0.1 + jnp.exp(aspect_ratio)
        sol   = self.evaporate(a_eff)
        V_L   = sol.ys[:, 0]
        t     = sol.ts
        return self.BOR(V_L, t)

    def optimize(self, verbose=True, t_final = 3600*24, max_iter = 100, lr = 1e-1, x0 = 1.0):
        train_loss = []

        self.time = t_final

        # Inicializar optimizador para el parámetro aspect_ratio (escalares float64)
        opt_init, opt_update, get_params = optimizers.adam(lr)
        opt_state = opt_init(jnp.array(x0, dtype=jnp.float64))  # valor inicial aspect ratio = 1.0

        for i in range(max_iter):
            params          = get_params(opt_state)
            loss_val, grads = jax.value_and_grad(self.objective_function)(params)
            opt_state       = opt_update(i, grads, opt_state)
            train_loss.append(loss_val)

            if verbose:
                print(f"Iter {i}, Loss {loss_val:.6e}, Param {params:.6e}")

        self.opt_state            = opt_state
        self.optimal_aspect_ratio = 0.1 + jnp.exp(get_params(opt_state))
        self.train_loss           = train_loss
        return self.optimal_aspect_ratio

    def plot_loss_history(self):
        plt.plot(self.train_loss)
        plt.xlabel('Iteration')
        plt.ylabel('BOR')
        plt.title('Training Loss History')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)
        plt.show()

        pass

    def plot_surface_respone(self, a_array, t_final):
        self.time = t_final
        BOR_values = jnp.array([self.objective_function(jnp.log(a) - 0.1) for a in a_array])
        plt.plot(a_array, BOR_values)
        plt.xlabel('Aspect Ratio')
        plt.ylabel('Boil-Off Rate (BOR)')
        plt.title('Response Surface of Boil-Off Rate vs Aspect Ratio')
        plt.grid(True)
        plt.show()
