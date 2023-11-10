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
        None: plo
    '''
    # Access to the liquid volume
    plt.plot(tank.sol.t, tank.sol.y[0])
    plt.grid()
    plt.xlabel('Time / s')
    plt.ylabel('$ V_L / m^3$')
    return

