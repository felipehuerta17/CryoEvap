import matplotlib.pyplot as plt

def plot_tv(tank):
    # Number of temperature profiles to visualise
    n_plots = int(tank.sol.t[-1] / tank.plot_interval)

    # Step to move in the index
    plot_step = int(tank.plot_interval / tank.time_interval)

    # Create a colormap
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
        ax.text(1.02, ((i-1 + 0.15)*plot_step) / len(tank.sol.t), f't={tank.sol.t[i*plot_step]:.0f} s', transform=ax.transAxes, 
                verticalalignment='center', bbox=dict(boxstyle='round,pad=0.5', edgecolor='none', facecolor=cmap(norm(tank.sol.t[i*plot_step])), alpha=0.6))

    # Add a grid
    ax.grid(True)

    # Add labels
    ax.set_ylabel(r'Dimensionless length $\zeta = z/l_V$')
    ax.set_xlabel('Temperature / K')

    # Add a title
    ax.set_title('Vapour temperature profiles at different times')

    # Show the plot
    plt.show()

def plot_V_L(tank):
    '''
    Plots liquid volume
    Inputs:
        Tank object with a sol object produced by the evaporate() function
    
    Returns:
        None:
    '''
    # Access to the liquid volume
    plt.plot(tank.sol.t, tank.sol.y[0])
    plt.grid()
    plt.xlabel('Time / s')
    plt.ylabel('$ V_L / m^3$')
    return

def plot_BOG(tank, unit='kg/h'):
    '''
    Plots boil-off gas and evaporation rate

    Inputs:
        tank: Tank object with a sol object produced by the evaporate() function
        unit: BOG units. Default: kg/h
        Options: kg/h, kg/s, g/s 
    
    Returns:
        None:
    '''
    
    # Conversion factors for plotting
    unit_conv = {'kg/h': 3600, 'kg/s': 1, 'g/s': 1000}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h

    plt.plot(tank.sol.t, tank.data['B_L'] * unit_conv[unit], label='Evaporation rate, $\dot{B}_L$')
    plt.plot(tank.sol.t[1:], tank.data['BOG'][1:] * unit_conv[unit], label='Boil-off gas rate, $\dot{B}$') 
    plt.grid()
    plt.xlabel('Time / s')
    plt.ylabel('Mass flow $/$ ' + unit)
    plt.legend()
    plt.show()
    return

def plot_Q(tank, unit='kW'):
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

    # Conversion factors for plotting
    unit_conv = {'W': 1, 'kW': 1e-3}

    fig, ax = plt.subplots(1, 3, figsize = [9,3])

    # Create space to breathe
    plt.subplots_adjust(wspace=0.5)

    # Q_VL plot
    ax[0].plot(tank.sol.t, (tank.data['Q_VL'] * unit_conv[unit]), label="Q_VL")
    ax[0].set_ylabel("$\dot{Q}_{VL} / kW $")
    ax[0].set_xlabel("Time / s")

    # Q_L_in plot
    ax[1].plot(tank.sol.t, (tank.data['Q_L']* unit_conv[unit]))
    ax[1].set_ylabel("$\dot{Q}_{L,in} / kW $")
    ax[1].set_xlabel("Time / s")

    ax[2].plot(tank.sol.t, (tank.data['Q_V'] * unit_conv[unit]))
    ax[2].set_ylabel("$\dot{Q}_{V,in} / kW $")
    ax[2].set_xlabel("Time / s")

    ax[1].set_title("Heat ingresses")

    [axis.grid() for axis in ax]
    plt.show()

