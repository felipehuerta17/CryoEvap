import jax
import jax.numpy as jnp
from jax.example_libraries import optimizers
from diffrax import diffeqsolve, ODETerm, Tsit5, SaveAt, PIDController, BacksolveAdjoint
import pandas as pd
import matplotlib.pyplot as plt
jax.config.update("jax_enable_x64", True)

# Import polyfit function from CoolProp
folder   = '../Cryoevap/cryogens/Coeffs/'
cp_V_df  = pd.read_csv(folder + 'coeffs_cpV.csv')
k_V_df   = pd.read_csv(folder + 'coeffs_kV.csv')
rho_V_df = pd.read_csv(folder + 'coeffs_rhoV.csv')

class Opti_jax:
    """
    Class Opti_jax
    ------------------
    The Opti_jax class is designed to perform JAX-accelerated optimization on a cryogenic storage tank system.
    It uses automatic differentiation for gradient-based optimization of the tank's aspect ratio to minimize
    boil-off rate (BOR). The class leverages JAX for both performance and differentiability of the entire
    simulation pipeline.

    Parameters
    ----------
    tank_obj : Tank
        An instance of the Tank class representing the cryogenic tank to be optimized.
        Must have all required heat transfer properties and a defined cryogen.

    Attributes
    ----------
    tank : Tank
        The provided Tank object to be optimized.
    params : dict
        Dictionary containing all the parameters needed for the simulation, extracted from the tank object.
    time : float
        Simulation time in seconds, set during optimization.
    opt_state : OptState
        JAX optimizer state after optimization.
    optimal_aspect_ratio : float
        The computed optimal aspect ratio after optimization.
    train_loss : list
        History of the objective function values (BOR) during optimization.

    Methods
    -------
    make_params(self, tank)
        Static method to extract and format all necessary parameters from the tank object.
    cp_V_fun(T, p), k_V_fun(T, p), rho_V_fun(T, p)
        JIT-compiled functions to calculate vapor specific heat capacity, thermal conductivity,
        and density using polynomial coefficients.
    sys_isobaric_jax(t, y, args)
        JIT-compiled ODE system representing the isobaric evaporation process.
    BOR(V_L, t)
        JIT-compiled function to calculate the boil-off rate from liquid volume and time.
    evaporate(aspect_ratio)
        Simulates the evaporation process for a given aspect ratio using diffrax solvers.
    objective_function(aspect_ratio)
        Calculates the boil-off rate for a given aspect ratio, used as the optimization objective.
    optimize(verbose=True, t_final=3600*24, max_iter=100, lr=1e-1, x0=1.0)
        Performs gradient-based optimization to find the aspect ratio that minimizes boil-off rate.
        Returns the optimal aspect ratio found.
    plot_loss_history()
        Plots the convergence history of the optimization process.
    plot_surface_respone(a_array, t_final)
        Plots the response surface of the boil-off rate as a function of the aspect ratio.
    """

    def __init__(self, tank_obj):
        self.tank   = tank_obj
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
            "cp_V_poly":  jnp.array(cp_V_df[cryo.name].values, dtype=jnp.float64),
            "k_V_poly":   jnp.array(k_V_df[cryo.name].values, dtype=jnp.float64),
            "rho_V_poly": jnp.array(rho_V_df[cryo.name].values, dtype=jnp.float64)
        }
    # Funciones internas adaptadas de tu código
    @staticmethod
    @jax.jit
    def cp_V_fun(T, p):
        """
        Calculate the specific heat capacity of the vapor at temperature T using polynomial coefficients p.
        
        Parameters
        ----------
        T : jnp.ndarray
            Temperature array.
        p : jnp.ndarray
            Polynomial coefficients for specific heat capacity.
        
        Returns
        -------
        jnp.ndarray
            Specific heat capacity evaluated at temperature T.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def k_V_fun(T, p):
        """
        Calculate the thermal conductivity of the vapor at temperature T using polynomial coefficients p.
        
        Parameters
        ----------
        T : jnp.ndarray
            Temperature array.
        p : jnp.ndarray
            Polynomial coefficients for thermal conductivity.
        
        Returns
        -------
        jnp.ndarray
            Thermal conductivity evaluated at temperature T.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def rho_V_fun(T, p):
        """
        Calculate the density of the vapor at temperature T using polynomial coefficients p.
        
        Parameters
        ----------
        T : jnp.ndarray
            Temperature array.
        p : jnp.ndarray
            Polynomial coefficients for density.
        
        Returns
        -------
        jnp.ndarray
            Density evaluated at temperature T.
        """
        return jnp.polyval(p, T)

    @staticmethod
    @jax.jit
    def sys_isobaric_jax(t, y, args):
        """
        ODE system representing the isobaric evaporation process.
        
        Parameters
        ----------
        t : float
            Current time (not used in this system).
        y : jnp.ndarray
            State vector containing liquid volume and temperatures.
        args : tuple
            Tuple containing aspect ratio and parameters dictionary.
        
        Returns
        -------
        jnp.ndarray
            Derivatives of the state vector.
        """

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
        """
        Calculate the Boil-Off Rate (BOR) per day from liquid volume and time.

        Parameters
        ----------
        V_L : jnp.ndarray
            Array of liquid volumes at different time steps.
        t : jnp.ndarray
            Array of time steps corresponding to the liquid volumes.
        Returns
        -------
        float
            Boil-Off Rate (BOR) per day, calculated as the fraction of initial liquid volume lost per day.
        """

        return (1.0 - V_L[-1] / V_L[0]) * (86400.0 / t[-1])

    def evaporate(self, aspect_ratio):
        """
        Simulates the evaporation process for a given aspect ratio using diffrax solvers.
        
        Parameters
        ----------
        aspect_ratio : float
            Aspect ratio to be used in the simulation, defined as the ratio of height to diameter.
        
        Returns
        -------
        sol : diffrax.Solution
            Solution object containing the results of the ODE integration.
        """

        VL_0 = jnp.array(self.params['V'] * self.params['LF'], dtype=jnp.float64)
        Tv_0 = jnp.ones(len(self.params['z_grid']), dtype=jnp.float64) * self.params['T_sat']
        IC   = jnp.concatenate([jnp.array([VL_0], dtype=jnp.float64), Tv_0])

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
        """
        Objective function to minimize the Boil-Off Rate (BOR) with respect to the aspect ratio.
        This function simulates the evaporation process for a given aspect ratio and calculates the BOR.
        
        Parameters
        ----------
        aspect_ratio : float
            Aspect ratio to be used in the simulation, defined as the ratio of height to diameter.
        
        Returns
        -------
        float
            The Boil-Off Rate (BOR) calculated from the simulation results.
        """

        a_eff = jnp.exp(aspect_ratio)
        sol   = self.evaporate(a_eff)
        V_L   = sol.ys[:, 0]
        t     = sol.ts
        return self.BOR(V_L, t)

    def optimize(self, verbose=True, t_final = 3600*24, max_iter = 100, lr = 1e-1, x0 = 1.0, lr_2 = 1e-3, iter_decay_lr = 100 ):
        """
        Performs gradient-based optimization to find the aspect ratio that minimizes boil-off rate (BOR).
        Save the loss history, optimal state and optimal aspect ratio in the class attributes.

        Parameters
        ----------
        verbose : bool, optional
            If True, prints optimization progress. Default is True.
        t_final : float, optional
            Final simulation time in seconds. Default is 3600*24 (one day).
        max_iter : int, optional
            Maximum number of optimization iterations. Default is 100.
        lr : float, optional
            Learning rate for the optimizer. Default is 1e-1.
        x0 : float, optional
            Initial guess for the aspect ratio, defined as a logarithmic scale. Default is 1.0.
        
        Returns
        -------
        float
            The optimal aspect ratio found during the optimization process.
        
        """
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
                print(f"Iter {i}, Loss {loss_val:.6e}, Param {float(jnp.exp(params)):.4g}")

            if i == iter_decay_lr:
                opt_init, opt_update, get_params = optimizers.adam(lr_2)
                opt_state = opt_init(get_params(opt_state))



        self.opt_state            = opt_state
        self.optimal_aspect_ratio = jnp.exp(get_params(opt_state))       
        self.train_loss           = train_loss
        return self.optimal_aspect_ratio

    def optimize_grid_with_refinement(self, verbose=True, t_final=3600*24, 
                                    coarse_samples=100, fine_samples=100,
                                    aspect_ratio_min=0.2, aspect_ratio_max=3.0, 
                                    refinement_window=0.2):
        """
        Two-phase optimization: coarse grid search followed by fine grid search around the best point.
        """
        self.time = t_final
        
        # Phase 1: Coarse grid search
        if verbose:
            print(f"Phase 1: Coarse grid search with {coarse_samples} samples...")
        
        aspect_ratios  = jnp.linspace(aspect_ratio_min, aspect_ratio_max, coarse_samples)
        bor_values     = jax.vmap(lambda a: self.objective_function(jnp.log(a)))(aspect_ratios)
        
        min_idx        = jnp.argmin(bor_values)
        coarse_optimal = aspect_ratios[min_idx]
        
        if verbose:
            print(f"Coarse search optimal: {coarse_optimal:.6f}, BOR: {bor_values[min_idx]:.6e}")
        
        # Phase 2: Fine grid search around the best point
        if verbose:
            print(f"Phase 2: Fine grid search with {fine_samples} samples...")
        
        # Define refined search range
        refined_min = max(aspect_ratio_min, coarse_optimal - refinement_window/2)
        refined_max = min(aspect_ratio_max, coarse_optimal + refinement_window/2)
        
        refined_aspect_ratios = jnp.linspace(refined_min, refined_max, fine_samples)
        refined_bor_values    = jax.vmap(lambda a: self.objective_function(jnp.log(a)))(refined_aspect_ratios)
        
        refined_min_idx = jnp.argmin(refined_bor_values)
        optimal_aspect_ratio = refined_aspect_ratios[refined_min_idx]
        min_bor = refined_bor_values[refined_min_idx]
        
        # Store results
        self.optimal_aspect_ratio = optimal_aspect_ratio
        # Combine coarse and fine evaluations for visualization
        self.train_loss = bor_values.tolist() + refined_bor_values.tolist()
        
        if verbose:
            print(f"Refined optimal aspect ratio: {optimal_aspect_ratio:.6f}")
            print(f"Minimum BOR: {min_bor:.6e}")
        
        return optimal_aspect_ratio


    def plot_loss_history(self):
        """
        Plots the convergence history of the optimization process.
        This method generates a log-log plot of the training loss (Boil-Off Rate) over iterations.
        """
        plt.plot(self.train_loss)
        plt.xlabel('Iteration')
        plt.ylabel('BOR')
        plt.title('Training Loss History')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True)

        pass

    def plot_surface_response(self, a_array, t_final):
        """
        Plots the response surface of the boil-off rate (BOR) as a function of the aspect ratio.
        
        Parameters
        ----------
        a_array : jnp.ndarray
            Array of aspect ratios for which to compute the boil-off rates.
        t_final : float
            Final simulation time in seconds, used to set the time for the evaporation simulation.
        """

        self.time = t_final
        BOR_values = jax.vmap(lambda a: self.objective_function(jnp.log(a)))(a_array)
        plt.plot(a_array, BOR_values, label=r't_final = ' + str(t_final/3600) + ' h')
        plt.xlabel('Aspect Ratio')
        plt.ylabel('Boil-Off Rate (BOR)')
        plt.title('Response Surface of Boil-Off Rate vs Aspect Ratio')
        plt.legend()
        plt.grid(True)

        pass