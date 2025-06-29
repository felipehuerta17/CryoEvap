import jax
import jax.numpy as jnp
from jax.example_libraries import optimizers
from diffrax import diffeqsolve, ODETerm, Tsit5, SaveAt, PIDController, BacksolveAdjoint
import pandas as pd
import matplotlib.pyplot as plt
jax.config.update("jax_enable_x64", True)

# Import polyfit function from CoolProp
folder   = '../cryoevap/cryogens/Coeffs/'
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
            "rho_V_poly": jnp.array(rho_V_df[cryo.name].values, dtype=jnp.float64),
            "q_b_fixed": tank.q_b_fixed
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
    def q_b_fun(p):
        if p["q_b_fixed"] is None:
            "If q_b_fixed is not set, calculate"
            return p["U_L"] * (p["T_air"] - p["T_sat"])
        else:
            return p["q_b_fixed"]
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
        Q_b   = Opti_jax.q_b_fun(p) * A_T
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

    def plot_surface_response_liquid_filling(self, a_array, lf_array, t_final):
        """
        Plots the response surface of the boil-off rate (BOR) as a function of the aspect ratio, for each
        liquid filling provided.
        
        Parameters
        ----------
        a_array : jnp.ndarray
            Array of aspect ratios for which to compute the boil-off rates.
        lf_array : jnp.ndarray
            Array of liquid filling for which to compute the tank.
        t_final : float
            Final simulation time in seconds, used to set the time for the evaporation simulation.
        """
        LF_og = self.params['LF']
        opt_a_values = jnp.array([])
        opt_bor_values = jnp.array([])
        plt.figure()
        for LF in lf_array:
            self.tank.LF = LF
            self.params = self.make_params(self, self.tank)
            self.time = t_final
            BOR_values = jax.vmap(lambda a: self.objective_function(jnp.log(a)))(a_array)
            plt.plot(a_array, BOR_values, label=r'LF = ' + str(self.tank.LF) )
            aspect_ratio = self.optimize_grid_with_refinement(verbose=False, t_final=self.time, coarse_samples=100, fine_samples=100,
                                   aspect_ratio_min=0.2, aspect_ratio_max=3, refinement_window=0.1)
            optimal_BOR = self.objective_function(jnp.log(aspect_ratio))
            opt_a_values = jnp.append(opt_a_values, aspect_ratio)
            opt_bor_values = jnp.append(opt_bor_values, optimal_BOR)
        print(f"optimal Aspect Ratio: {opt_a_values}")
        print(f"optimal BOR: {opt_bor_values}")
        plt.plot(opt_a_values, opt_bor_values, color='red', label='optimal values',linestyle='--')
        plt.xlabel('Aspect Ratio')
        plt.ylabel('Boil-Off Rate (BOR)')
        plt.title('Response Surface of Boil-Off Rate vs Aspect Ratio | t=' + str(t_final/3600) + ' h')
        plt.grid(True)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.axis('tight')
        self.tank.LF = LF_og
        self.params = self.make_params(self, self.tank)
        return plt.show()
    
    def plot_surface_response_bottom_heat(self, q_b_array, t_final):
        """
        Plots the response surface of the optimized aspect ratio as a function of the bottom heat flux.
        
        Parameters
        ----------
        q_b_array : jnp.ndarray
            Array of bottom heat fluxes for which to compute the boil-off rates.
        t_final : float
            Final simulation time in seconds, used to set the time for the evaporation simulation.
        """
        # Redefinir el tanke para cada q_b
        q_b_og = self.params['q_b_fixed']
        a_values = jnp.array([])
        self.time = t_final
        for q_b in q_b_array:
            self.tank.q_b_fixed = q_b
            self.params = self.make_params(self, self.tank)
        # calcular el optimo para cada q_b
        # plotear
            optimal_aspect_ratio = self.optimize_grid_with_refinement(verbose=False, t_final=720*3600, coarse_samples=100, fine_samples=100,
                                   aspect_ratio_min=0.2, aspect_ratio_max=3, refinement_window=0.1)
            a_values = jnp.append(a_values, optimal_aspect_ratio)
        plt.plot(q_b_array, a_values)
        plt.xlabel('Heat flux | w/m^2')
        plt.ylabel('Aspect Ratio')
        plt.title('Response Surface of optimized Aspect Ratio vs Heat Flux | t=' + str(t_final/3600) + ' h')
        plt.grid(True)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.axis('tight')
        self.tank.q_b_fixed = q_b_og
        self.params = self.make_params(self, self.tank)
        return plt.show()


    @staticmethod
    @jax.jit
    def thermal_sys_isobaric_jax(t, y, args):
        """
        ODE system representing the isobaric evaporation process.
        
        Parameters
        ----------
        t : float
            Current time (not used in this system).
        y : jnp.ndarray
            State vector containing liquid volume and temperatures.
        args : tuple
            Tuple containing thermal aspect ratio and parameters dictionary.
        
        Returns
        -------
        jnp.ndarray
            Derivatives of the state vector.
        """

        thermal_aspect_ratio, p = args

        V_L = y[0]
        T_V = y[1:]

        aspect_ratio = (thermal_aspect_ratio*Opti_jax.q_b_fun(p))/( 4*1.02*(1-thermal_aspect_ratio)*((p["U_L"]*p["LF"]*(p["T_air"] - p["T_sat"]))+(p["U_V"]*(1-p["LF"])*p["eta_w"] * (p["T_air"] - jnp.mean(T_V)))))

        d_i = ((4 * p["V"])/(jnp.pi * aspect_ratio))**(1/3)
        d_o = d_i * 1.02 

        

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
        Q_b   = Opti_jax.q_b_fun(p)*A_T
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

    def thermal_evaporate(self, thermal_aspect_ratio):
        """
        Simulates the evaporation process for a given thermal aspect ratio using diffrax solvers.
        
        Parameters
        ----------
        thermal_aspect_ratio : float
            Thermal Aspect ratio to be used in the simulation, defined as the ratio of height to diameter.
        
        Returns
        -------
        sol : diffrax.Solution
            Solution object containing the results of the ODE integration.
        """

        VL_0 = jnp.array(self.params['V'] * self.params['LF'], dtype=jnp.float64)
        Tv_0 = jnp.ones(len(self.params['z_grid']), dtype=jnp.float64) * self.params['T_sat']
        IC   = jnp.concatenate([jnp.array([VL_0], dtype=jnp.float64), Tv_0])

        term = ODETerm(self.thermal_sys_isobaric_jax)

        sol = diffeqsolve(
            term,
            solver=Tsit5(),
            t0=0,
            t1=self.time,
            dt0=0.01,
            y0=IC,
            args=(thermal_aspect_ratio, self.params),
            saveat=SaveAt(ts=jnp.arange(0, self.time + 1, self.params['time_interval'])),
            max_steps=1000000,
            stepsize_controller=PIDController(rtol=1e-8, atol=1e-8),
            adjoint=BacksolveAdjoint()
        )
        return sol

    def thermal_objective_function(self, thermal_aspect_ratio):
        """
        Objective function to minimize the Boil-Off Rate (BOR) with respect to the thermal aspect ratio.
        This function simulates the evaporation process for a given aspect ratio and calculates the BOR.
        
        Parameters
        ----------
        thermal_aspect_ratio : float
            Thermal Aspect ratio to be used in the simulation, defined as the ratio of height to diameter.
        
        Returns
        -------
        float
            The Boil-Off Rate (BOR) calculated from the simulation results.
        """

        a_eff = jnp.exp(thermal_aspect_ratio)
        sol   = self.thermal_evaporate(a_eff)
        V_L   = sol.ys[:, 0]
        t     = sol.ts
        return self.BOR(V_L, t)

    def thermal_optimize(self, verbose=True, t_final = 3600*24, max_iter = 100, lr = 1e-1, x0 = 1.0, lr_2 = 1e-3, iter_decay_lr = 100 ):
        """
        Performs gradient-based optimization to find the thermal aspect ratio that minimizes boil-off rate (BOR).
        Save the loss history, optimal state and optimal thermal aspect ratio in the class attributes.

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
            Initial guess for the thermal aspect ratio, defined as a logarithmic scale. Default is 1.0.
        
        Returns
        -------
        float
            The optimal thermal aspect ratio found during the optimization process.
        
        """
        train_loss = []

        self.time = t_final

        # Inicializar optimizador para el parámetro aspect_ratio (escalares float64)
        opt_init, opt_update, get_params = optimizers.adam(lr)
        opt_state = opt_init(jnp.array(x0, dtype=jnp.float64))  # valor inicial aspect ratio = 1.0

        for i in range(max_iter):
            params          = get_params(opt_state)
            loss_val, grads = jax.value_and_grad(self.thermal_objective_function)(params)
            opt_state       = opt_update(i, grads, opt_state)
            train_loss.append(loss_val)

            if verbose:
                print(f"Iter {i}, Loss {loss_val:.6e}, Param {float(jnp.exp(params)):.4g}")

            if i == iter_decay_lr:
                opt_init, opt_update, get_params = optimizers.adam(lr_2)
                opt_state = opt_init(get_params(opt_state))



        self.thermal_opt_state            = opt_state
        self.optimal_thermal_aspect_ratio = jnp.exp(get_params(opt_state))       
        self.thermal_train_loss           = train_loss
        return self.optimal_thermal_aspect_ratio

    def thermal_optimize_grid_with_refinement(self, verbose=True, t_final=3600*24, 
                                    coarse_samples=100, fine_samples=100,
                                    thermal_aspect_ratio_min=0.01, thermal_aspect_ratio_max=1, 
                                    refinement_window=0.2):
        """
        Two-phase optimization: coarse grid search followed by fine grid search around the best point.
        """
        self.time = t_final
        
        # Phase 1: Coarse grid search
        if verbose:
            print(f"Phase 1: Coarse grid search with {coarse_samples} samples...")
        
        thermal_aspect_ratios  = jnp.linspace(thermal_aspect_ratio_min, thermal_aspect_ratio_max, coarse_samples)
        bor_values     = jax.vmap(lambda a: self.thermal_objective_function(jnp.log(a)))(thermal_aspect_ratios)
        
        min_idx        = jnp.argmin(bor_values)
        coarse_optimal = thermal_aspect_ratios[min_idx]
        
        if verbose:
            print(f"Coarse search optimal: {coarse_optimal:.6f}, BOR: {bor_values[min_idx]:.6e}")
        
        # Phase 2: Fine grid search around the best point
        if verbose:
            print(f"Phase 2: Fine grid search with {fine_samples} samples...")
        
        # Define refined search range
        refined_min = max(thermal_aspect_ratio_min, coarse_optimal - refinement_window/2)
        refined_max = min(thermal_aspect_ratio_max, coarse_optimal + refinement_window/2)
        
        refined_thermal_aspect_ratios = jnp.linspace(refined_min, refined_max, fine_samples)
        refined_bor_values    = jax.vmap(lambda a: self.thermal_objective_function(jnp.log(a)))(refined_thermal_aspect_ratios)
        
        refined_min_idx = jnp.argmin(refined_bor_values)
        optimal_thermal_aspect_ratio = refined_thermal_aspect_ratios[refined_min_idx]
        min_bor = refined_bor_values[refined_min_idx]
        
        # Store results
        self.optimal_thermal_aspect_ratio = optimal_thermal_aspect_ratio
        # Combine coarse and fine evaluations for visualization
        self.thermal_train_loss = bor_values.tolist() + refined_bor_values.tolist()
        
        if verbose:
            print(f"Refined optimal thermal aspect ratio: {optimal_thermal_aspect_ratio:.6f}")
            print(f"Minimum BOR: {min_bor:.6e}")
        
        return optimal_thermal_aspect_ratio
    def plot_thermal_surface_response(self, thermal_a_array, t_final):
        """
        Plots the response surface of the boil-off rate (BOR) as a function of the thermal aspect ratio.
        
        Parameters
        ----------
        thermal_a_array : jnp.ndarray
            Array of thermal aspect ratios for which to compute the boil-off rates.
        t_final : float
            Final simulation time in seconds, used to set the time for the evaporation simulation.
        """

        self.time = t_final
        BOR_values = jax.vmap(lambda a: self.thermal_objective_function(jnp.log(a)))(thermal_a_array)
        plt.plot(thermal_a_array, BOR_values, label=r't_final = ' + str(t_final/3600) + ' h')
        plt.xlabel('Thermal Aspect Ratio')
        plt.ylabel('Boil-Off Rate (BOR)')
        plt.title('Response Surface of Boil-Off Rate vs Thermal Aspect Ratio')
        plt.legend()
        plt.grid(True)

        pass
    def plot_thermal_aspect_ratio(self, thermal_a_array, log = False):
        a_array = []
        p = self.params
        T_V = jnp.ones(len(self.params['z_grid']), dtype=jnp.float64) * self.params['T_sat']
        for thermal_aspect_ratio in thermal_a_array:
            aspect_ratio = (thermal_aspect_ratio*Opti_jax.q_b_fun(p))/( 4*1.02*(1-thermal_aspect_ratio)*((p["U_L"]*p["LF"]*(p["T_air"] - p["T_sat"]))+(p["U_V"]*(1-p["LF"])*p["eta_w"] * (p["T_air"] - jnp.mean(T_V)))))
            a_array.append(aspect_ratio)
        if log:
            plt.plot(thermal_a_array, jnp.log(jnp.array(a_array)))
        else:
            plt.plot(thermal_a_array, jnp.array(a_array))
        plt.xlabel('Thermal Aspect Ratio')
        plt.ylabel('Aspect Ratio')
        plt.title('Response Surface of Aspect Ratio vs Thermal Aspect Ratio')
        plt.legend()
        plt.grid(True)

    def plot_surface_response_thermal_liquid_filling(self, thermal_a_array, lf_array, t_final):
        """
        Plots the response surface of the boil-off rate (BOR) as a function of the thermal aspect ratio, for each
        liquid filling provided.
        
        Parameters
        ----------
        thermal_a_array : jnp.ndarray
            Array of thermal aspect ratios for which to compute the boil-off rates.
        lf_array : jnp.ndarray
            Array of liquid filling for which to compute the tank.
        t_final : float
            Final simulation time in seconds, used to set the time for the evaporation simulation.
        """
        LF_og = self.params['LF']
        opt_a_values = jnp.array([])
        opt_bor_values = jnp.array([])
        plt.figure()
        for LF in lf_array:
            self.tank.LF = LF
            self.params = self.make_params(self, self.tank)
            self.time = t_final
            BOR_values = jax.vmap(lambda a: self.thermal_objective_function(jnp.log(a)))(thermal_a_array)
            plt.plot(thermal_a_array, BOR_values, label=r'LF = ' + str(self.tank.LF) )
            thermal_aspect_ratio = self.thermal_optimize_grid_with_refinement(verbose=False, t_final=self.time, coarse_samples=100, fine_samples=1000,
                                   thermal_aspect_ratio_min=0.1, thermal_aspect_ratio_max=0.9, refinement_window=0.1)
            optimal_BOR = self.thermal_objective_function(jnp.log(thermal_aspect_ratio))
            opt_a_values = jnp.append(opt_a_values, thermal_aspect_ratio)
            opt_bor_values = jnp.append(opt_bor_values, optimal_BOR)
        print(f"optimal Thermal Aspect Ratio: {opt_a_values}")
        print(f"optimal BOR: {opt_bor_values}")
        plt.plot(opt_a_values, opt_bor_values, color='red', label='optimal values',linestyle='--')
        plt.xlabel('Thermal Aspect Ratio')
        plt.ylabel('Boil-Off Rate (BOR)')
        plt.title('Response Surface of Boil-Off Rate vs Thermal Aspect Ratio | t=' + str(t_final/3600) + ' h')
        plt.grid(True)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.axis('tight')
        self.tank.LF = LF_og
        self.params = self.make_params(self, self.tank)
        return plt.show()