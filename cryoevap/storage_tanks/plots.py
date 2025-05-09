import matplotlib.pyplot as plt
import numpy as np

def plot_tv(tank, t_unit='s', hour = 0):

    # Parámetros configurables
    fixed_hour = hour  # Hora fija (8 AM)
    n_days = 5         # Número de curvas a graficar

    # Factores de conversión de unidades de tiempo
    t_dict = {'s': 1, 'min': 60, 'h': 3600, 'd': 3600 * 24, 'w': 3600 * 24 * 7}

    # Tiempo total en días
    t_final_days = tank.sol.t[-1] / t_dict['d']

    # Seleccionar los tiempos de muestreo
    if t_final_days >= 5:
        selected_days = np.round(np.linspace(0, t_final_days, n_days)).astype(int)  # Días discretos
        selected_times = (selected_days * t_dict['d']) + fixed_hour  # Convertir a segundos y fijar hora
    else:
        selected_times = np.linspace(0, tank.sol.t[-1], n_days)  # Equiespaciado si < 5 días

    # Crear mapa de colores
    cmap = plt.get_cmap('cividis')
    norm = plt.Normalize(vmin=tank.sol.t[1], vmax=tank.sol.t[-1])

    # Crear figura
    fig, ax = plt.subplots(figsize=[7, 5])

    # Graficar condición inicial
    T_init = tank.sol.y[1:len(tank.z_grid) + 1, 0]
    ax.plot(T_init, tank.z_grid, '--k', label='Initial Condition')

    # Graficar curvas de temperatura para los tiempos seleccionados
    for i, t_sel in enumerate(selected_times):
        idx = np.searchsorted(tank.sol.t, t_sel)  # Buscar el índice más cercano

        if idx >= len(tank.sol.t):
            idx = len(tank.sol.t) - 1  # Evitar desbordamientos

        # Obtener perfil de temperatura
        T_v = tank.sol.y[1:len(tank.z_grid) + 1, idx]

        # Graficar con color correspondiente al tiempo
        ax.plot(T_v, tank.z_grid, color=cmap(norm(tank.sol.t[idx])))

        # Agregar etiqueta de tiempo al lado del gráfico
        ax.text(1.02, (i + 0.5) / n_days, f't={tank.sol.t[idx]/t_dict[t_unit]:.1f} {t_unit}', 
                transform=ax.transAxes, verticalalignment='center',
                bbox=dict(boxstyle='round,pad=0.5', edgecolor='none',
                        facecolor=cmap(norm(tank.sol.t[idx])), alpha=0.6))

    # Configuración del gráfico
    ax.grid(True)
    ax.set_ylabel(r'Dimensionless length $\zeta = z/l_V$')
    ax.set_xlabel('Temperature / K')

    plt.show()


def plot_V_L(tank, unit='m3', t_unit='s'):
    '''
    Plots liquid volume
    Inputs:
        Tank object with a sol object produced by the evaporate() function
        unit: Liquid volume units. Default: m3
        t_unit: Time unit. Default: s
        Options: m^3, L, mL 
    
    Returns:
        None:
    '''

    # Conversion factors for plotting
    unit_conv = {'m3': 1, 'L': 1e3, 'mL': 1e6}

        # Time unit conversion factor for plotting
    t_dict = {'s':1, 'min': 60, 'h':3600, 'd': 3600*24, 'w': 3600*24*7}

    # Create a colormap
    cmap = plt.get_cmap('cividis')

    plt.figure(figsize=[7, 5])
    # Access to the liquid volume
    plt.plot(tank.sol.t/t_dict[t_unit], tank.sol.y[0] * unit_conv[unit], color = cmap(1/6))
    plt.grid()
    plt.xlabel('Time / ' + t_unit)
    if unit == "m3":
        plt.ylabel(r'$V_L$ / $m^3$')
    else:
        plt.ylabel(r'$V_L$ / ' + unit)
    return

def plot_BOG(tank, unit='kg/h', t_unit = 's'):
    '''
    Plots boil-off gas and evaporation rate

    Inputs:
        tank: Tank object with a sol object produced by the evaporate() function
        unit: BOG units. Default: kg/h
        t_unit: Time units. Default: s
        Options: kg/h, kg/s, g/s 
    
    Returns:
        None:
    '''

    # Create a colormap
    cmap = plt.get_cmap('cividis')
    
    # Conversion factors for plotting
    unit_conv = {'kg/h': 3600, 'kg/s': 1, 'g/h' : 3600*1e3, 'g/s': 1000}

    # Time unit conversion factor for plotting
    t_dict = {'s':1, 'min': 60, 'h':3600, 'd': 3600*24, 'w': 3600*24*7}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h

    plt.figure(figsize=[7, 5])
    plt.plot(tank.sol.t/t_dict[t_unit], tank.data['B_L'] * unit_conv[unit], label=r'Evaporation rate, $\dot{B}_L$', color = cmap(1/6))
    plt.plot(tank.sol.t[1:]/t_dict[t_unit], tank.data['BOG'][1:] * unit_conv[unit], label=r'Boil-off gas rate, $\dot{B}$', color = cmap(5/6)) 
    plt.grid()
    plt.xlabel('Time / ' + t_unit)
    plt.ylabel(r'Mass flow $/$ ' + unit)
    plt.legend()
    return

def plot_Q(tank, unit='kW', t_unit = 's'):
    '''
    Plots vapour to liquid heat transfer rate

    Inputs:
        tank: Tank object with a sol object produced by the 
        evaporate() function
        unit: Q_VL units. Default: kW
        Options: kW, W
    
    Returns:
        None:
    '''

    # Create a colormap
    cmap = plt.get_cmap('cividis')

    # Conversion factors for plotting
    unit_conv = {'W': 1, 'kW': 1e-3}

    # Time unit conversion factor for plotting
    t_dict = {'s':1, 'min': 60, 'h':3600, 'd': 3600*24, 'w': 3600*24*7}


    fig, ax = plt.subplots(2, 2, figsize = [7,7])

    # Create space to breathe
    plt.subplots_adjust(wspace=0.5)

    # Q_L_in plot
    ax[0][0].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_L']* unit_conv[unit]), color = cmap(1/6))
    ax[0][0].set_ylabel(r"$\dot{Q}_L$ / " + unit)
    ax[0][0].set_xlabel("Time / " + t_unit)
    ax[0][0].grid()

    # Q_V_in plot
    ax[1][0].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_V'] * unit_conv[unit]), color = cmap(1/6))
    ax[1][0].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_Vw'] * unit_conv[unit]), label="Q_Vw",color = cmap(5/6))
    ax[1][0].set_ylabel(r"$\dot{Q}_V$ /  " + unit)
    ax[1][0].set_xlabel("Time / " + t_unit)
    ax[1][0].grid()

    # Q_VL plot
    ax[0][1].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_VL'] * unit_conv[unit]), label="Q_VL", color = cmap(1/6))
    ax[0][1].set_ylabel(r"$\dot{Q}_{VL}$ / " + unit)
    ax[0][1].set_xlabel("Time / " + t_unit)
    ax[0][1].grid()

    # Q_{V,w} plot
    #ax[1][1].plot(tank.sol.t, (tank.data['Q_Vw'] * unit_conv[unit]), label="Q_Vw",color = cmap(1/6))
    ax[1][1].plot(tank.sol.t/t_dict[t_unit], ( (tank.data['Q_Vw'] + tank.data['Q_VL'] + tank.data['Q_L'] + tank.Q_b(tank.sol.t))  *
                                unit_conv[unit]), label="Q_{tot}",color = cmap(1/6))
    # ax[1][1].set_ylabel("$\dot{Q}_{V,w}$ / " + unit)
    ax[1][1].set_ylabel(r"$\dot{Q}_{tot}$ / " + unit)
    ax[1][1].set_xlabel("Time / " + t_unit)
    ax[1][1].grid()

    # ax[0][0].set_title("Heat ingresses")

    # [axis.grid() for axis in ax]


def plot_tv_BOG(tank, t_unit = 's'):
    '''
    Plots average vapour temperature and boil-off gas
    temperature as a function of time

    Inputs:
        tank: Tank object with a sol object produced by the 
        evaporate() function
        t_unit: Time units. Default: s.
        Options: s, min, h, days.
    
    Returns:
        None:
    '''

    # Time unit conversion factor for plotting
    t_dict = {'s':1, 'min': 60, 'h':3600, 'd': 3600*24, 'w': 3600*24*7}

    # Ad-hoc plotting of average and boil-off gas temperature
    plt.figure(figsize=[7,5])
    plt.plot(tank.data['Time']/t_dict[t_unit], tank.data['Tv_avg'], label=r'$\overline{T}_V$')
    plt.plot(tank.data['Time']/t_dict[t_unit], tank.data['T_BOG'], label=r'$T_{BOG}$')
    plt.xlabel("Time / " + t_unit)
    plt.ylabel('Temperature / K')
    plt.legend()
    plt.grid()


############################################################################################################
##### NEW FUNCTIONS #####
############################################################################################################


def plot_tw(tank, t_unit='s', hour = 0):
    '''
    Plots the wall temperature profiles at specific times.

    Inputs:
        tank: Tank object with a sol object produced by the evaporate() function.
        t_unit: Time units. Options: 's' (seconds), 'min' (minutes), 'h' (hours), 
                'd' (days), or 'w' (weeks). Default: 's'.
        hour: Fixed hour (in seconds) to adjust the time labels on the plot. Default: 0.

    Returns:
        None:
    '''

    # Configurable parameters
    fixed_hour = hour  # Fixed hour (8 AM)
    n_days = 5         # Number of curves to plot

    # Time unit conversion factors
    t_dict = {'s': 1, 'min': 60, 'h': 3600, 'd': 3600 * 24, 'w': 3600 * 24 * 7}

    # Total time in days
    t_final_days = tank.sol.t[-1] / t_dict['d']

    # Select sampling times
    if t_final_days >= 5:
        selected_days = np.round(np.linspace(0, t_final_days, n_days)).astype(int)  # Discrete days
        selected_times = (selected_days * t_dict['d']) + fixed_hour  # Convert to seconds and fix hour
    else:
        selected_times = np.linspace(0, tank.sol.t[-1], n_days)  # Equally spaced if < 5 days

    # Create colormap
    cmap = plt.get_cmap('cividis')
    norm = plt.Normalize(vmin=tank.sol.t[1], vmax=tank.sol.t[-1])

    # Create figure
    fig, ax = plt.subplots(figsize=[7, 5])

    # Plot temperature curves for selected times
    for i, t_sel in enumerate(selected_times):
        idx = np.searchsorted(tank.sol.t, t_sel)  # Find the closest index

        # print(tank.sol.t[idx]/3600)
        if idx >= len(tank.sol.t):
            idx = len(tank.sol.t) - 1  # Avoid overflow

        # Get temperature profile
        T_w = tank.sol.y[len(tank.z_grid) + 1:, idx]

        # Plot with color corresponding to time
        ax.plot(tank.r_grid * (tank.d_o - tank.d_i) * 0.5, T_w, color=cmap(norm(tank.sol.t[idx])))

        # Add time label next to the plot
        ax.text(1.02, (i + 0.5) / n_days, f't={tank.sol.t[idx]/t_dict[t_unit]:.1f} {t_unit}', 
                transform=ax.transAxes, verticalalignment='center',
                bbox=dict(boxstyle='round,pad=0.5', edgecolor='none',
                        facecolor=cmap(norm(tank.sol.t[idx])), alpha=0.6))

    # Plot configuration
    ax.grid(True)
    ax.set_xlabel(r'Radius $r$ / m')
    ax.set_ylabel('Temperature / K')

    plt.show()


def plot_Q_w(tank, unit = 'W', t_unit = 's'):
    '''
    Plots boil-off gas and evaporation rate

    Inputs:
        tank: Tank object with a sol object produced by the evaporate() function
        unit: BOG units. Default: kg/h
        t_unit: Time units. Default: s
        Options: kg/h, kg/s, g/s 
    
    Returns:
        None:
    '''

    # Create a colormap
    cmap = plt.get_cmap('cividis')
    
    # Time unit conversion factor for plotting
    t_dict = {'s':1, 'min': 60, 'h':3600, 'd': 3600*24, 'w': 3600*24*7}

    unit_conv = {'W': 1, 'kW': 1e-3}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h
    plt.figure(figsize=[7, 5])
    plt.plot(tank.sol.t/t_dict[t_unit], tank.data['Q_env_w']*unit_conv[unit] , label=r'Heat flow from the enviroment to the wall, $\dot{Q}_{W,env}$', color = cmap(1/6))
    plt.plot(tank.sol.t/t_dict[t_unit], tank.data['Q_w_L'] *unit_conv[unit], label=r'Heat flow from the wall to the ammonia, $\dot{Q}_{WL}$', color = cmap(5/6)) 
    plt.grid()
    plt.xlabel('Time / ' + t_unit)
    plt.ylabel(r'Heat flow $/$ ' + unit)
    plt.legend()
    return

def plot_T_w_avg(tank, t_unit = 's'):
    '''
    Plots average wall temperature and the environmental temperature

    Inputs:
        tank: Tank object with a sol object produced by the evaporate() function
        t_unit: Time units. Default: s
    
    Returns:
        None:
    '''
        
    # Create a colormap
    cmap = plt.get_cmap('cividis')
    
    # Time unit conversion factor for plotting
    t_dict = {'s':1, 'min': 60, 'h':3600, 'd': 3600*24, 'w': 3600*24*7}

    t_span = np.linspace(tank.sol.t[0], tank.sol.t[-1], 1000)
    T_air  = [tank.T_env(t) for t in t_span]

    plt.figure(figsize=[7, 5])
    plt.plot(tank.sol.t/t_dict[t_unit], tank.data['Tw_avg'] , label='Average wall temperature', color = cmap(1/6))
    plt.plot(t_span/t_dict[t_unit], T_air, label='Environmental temperature', color = cmap(5/6)) 
    plt.grid()
    plt.xlabel('Time / ' + t_unit)
    plt.ylabel(r'Temperature $/ K $')
    plt.legend()
    return

