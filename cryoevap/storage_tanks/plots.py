import matplotlib.pyplot as plt

def plot_tv(tank, t_unit='s'):
    # Number of temperature profiles to visualise
    n_plots = int(tank.sol.t[-1] / tank.plot_interval)

    # Time unit conversion factor for plotting
    t_dict = {'s':1, 'min': 60, 'h':3600, 'd': 3600*24, 'w': 3600*24*7}

    # Step to move in the index
    plot_step = int(tank.plot_interval / tank.time_interval)

    # Create a colour map
    cmap = plt.get_cmap('cividis')

    # Normalize the colormap based on the time range
    norm = plt.Normalize(vmin=tank.sol.t[1], vmax=tank.sol.t[-1])

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Loop over each time step to  p
    for i in range(1, n_plots+1):

        # Get the temperature at this time step
        T_v = tank.sol.y[1:, i*plot_step]
        
        # Plot the temperature profile at this time step, with color indicating the time
        ax.plot(T_v, tank.z_grid, color=cmap(norm(tank.sol.t[i * plot_step])))

        # Add a text box with the time value at the right of the plot
        ax.text(1.02, ((i-1 + 0.15)*plot_step) / len(tank.sol.t), f't={tank.sol.t[i*plot_step]/t_dict[t_unit]:.0f} ' + t_unit, transform=ax.transAxes, 
                verticalalignment='center', bbox=dict(boxstyle='round,pad=0.5', edgecolor='none', facecolor=cmap(norm(tank.sol.t[i*plot_step])), alpha=0.6))

    # Add a grid
    ax.grid(True)

    # Add labels
    ax.set_ylabel(r'Dimensionless length $\zeta = z/l_V$')
    ax.set_xlabel('Temperature / K')


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

    # Access to the liquid volume
    plt.plot(tank.sol.t/t_dict[t_unit], tank.sol.y[0] * unit_conv[unit], color = cmap(1/6))
    plt.grid()
    plt.xlabel('Time / ' + t_unit)
    if unit == "m3":
        plt.ylabel('$V_L$ / $m^3$')
    else:
        plt.ylabel('$V_L$ / ' + unit)
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

    plt.figure(figsize=[4,4])
    plt.plot(tank.sol.t/t_dict[t_unit], tank.data['B_L'] * unit_conv[unit], label='Evaporation rate, $\dot{B}_L$', color = cmap(1/6))
    plt.plot(tank.sol.t[1:]/t_dict[t_unit], tank.data['BOG'][1:] * unit_conv[unit], label='Boil-off gas rate, $\dot{B}$', color = cmap(5/6)) 
    plt.grid()
    plt.xlabel('Time / ' + t_unit)
    plt.ylabel('Mass flow $/$ ' + unit)
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


    fig, ax = plt.subplots(2, 2, figsize = [6,6])

    # Create space to breathe
    plt.subplots_adjust(wspace=0.5)

    # Q_L_in plot
    ax[0][0].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_L']* unit_conv[unit]), color = cmap(1/6))
    ax[0][0].set_ylabel("$\dot{Q}_L$ / " + unit)
    ax[0][0].set_xlabel("Time / " + t_unit)
    ax[0][0].grid()

    # Q_V_in plot
    ax[1][0].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_V'] * unit_conv[unit]), color = cmap(1/6))
    ax[1][0].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_Vw'] * unit_conv[unit]), label="Q_Vw",color = cmap(5/6))
    ax[1][0].set_ylabel("$\dot{Q}_V$ /  " + unit)
    ax[1][0].set_xlabel("Time / " + t_unit)
    ax[1][0].grid()

    # Q_VL plot
    ax[0][1].plot(tank.sol.t/t_dict[t_unit], (tank.data['Q_VL'] * unit_conv[unit]), label="Q_VL", color = cmap(1/6))
    ax[0][1].set_ylabel("$\dot{Q}_{VL}$ / " + unit)
    ax[0][1].set_xlabel("Time / " + t_unit)
    ax[0][1].grid()

    # Q_{V,w} plot
    #ax[1][1].plot(tank.sol.t, (tank.data['Q_Vw'] * unit_conv[unit]), label="Q_Vw",color = cmap(1/6))
    ax[1][1].plot(tank.sol.t/t_dict[t_unit], ( (tank.data['Q_Vw'] + tank.data['Q_VL'] + tank.data['Q_L'] + tank.Q_b)  *
                                unit_conv[unit]), label="Q_{tot}",color = cmap(1/6))
    # ax[1][1].set_ylabel("$\dot{Q}_{V,w}$ / " + unit)
    ax[1][1].set_ylabel("$\dot{Q}_{tot}$ / " + unit)
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
    plt.figure(figsize=[4,4])
    plt.plot(tank.data['Time']/t_dict[t_unit], tank.data['Tv_avg'], label=r'$\overline{T}_V$')
    plt.plot(tank.data['Time']/t_dict[t_unit], tank.data['T_BOG'], label=r'$T_{BOG}$')
    plt.xlabel("Time / " + t_unit)
    plt.ylabel('Temperature / K')
    plt.legend()
    plt.grid()

