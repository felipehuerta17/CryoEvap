import matplotlib.pyplot as plt

from matplotlib.ticker import FormatStrFormatter

import numpy as np
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
    ax.set_title('Vapour temperature profiles at different times - '+tank.Geo_v+' geometry')

    # Show the plot
    # plt.show()

def plot_V_L(tank, unit='m3'):
    '''
    Plots liquid volume
    Inputs:
        Tank object with a sol object produced by the evaporate() function
        unit: Liquid volume units. Default: m3
        Options: m^3, L, mL 
    
    Returns:
        None:
    '''

    # Conversion factors for plotting
    unit_conv = {'m3': 1, 'L': 1e3, 'mL': 1e6}

    # Create a colormap
    cmap = plt.get_cmap('cividis')


    # Access to the liquid volume
    plt.plot(tank.sol.t, tank.sol.y[0] * unit_conv[unit], color = cmap(1/6))
    plt.grid()
    plt.xlabel('Time / s')
    if unit == "m3":
        plt.ylabel('$V_L$ / $m^3$')
    else:
        plt.ylabel('$V_L$ / ' + unit)
    plt.title('Liquid volume over time - '+tank.Geo_l+' geometry')
    plt.show()
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

    # Create a colormap
    cmap = plt.get_cmap('cividis')
    
    # Conversion factors for plotting
    unit_conv = {'kg/h': 3600, 'kg/s': 1, 'g/h' : 3600*1e3, 'g/s': 1000}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h

    plt.plot(tank.sol.t, tank.data['B_L'] * unit_conv[unit], label='Evaporation rate, $\dot{B}_L$', color = cmap(1/6))
    plt.plot(tank.sol.t[1:], tank.data['BOG'][1:] * unit_conv[unit], label='Boil-off gas rate, $\dot{B}$', color = cmap(5/6)) 
    plt.grid()
    plt.xlabel('Time / s')
    plt.ylabel('Mass flow $/$ ' + unit)
    plt.title('BOG and evaporation rates over time - '+tank.Geo_l+' geometry')
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

    # Create a colormap
    cmap = plt.get_cmap('cividis')

    # Conversion factors for plotting
    unit_conv = {'W': 1, 'kW': 1e-3}

    fig, ax = plt.subplots(2, 2, figsize = [6,6])

    # Create space to breathe
    plt.subplots_adjust(wspace=0.5)

    # Q_L_in plot
    ax[0][0].plot(tank.sol.t, (tank.data['Q_L']* unit_conv[unit]), color = cmap(1/6))
    ax[0][0].set_ylabel("$\dot{Q}_L$ / " + unit)
    ax[0][0].set_xlabel("Time / s")
    ax[0][0].grid()

    # Q_V_in plot
    ax[0][1].plot(tank.sol.t, (tank.data['Q_V'] * unit_conv[unit]), color = cmap(1/6))
    ax[0][1].set_ylabel("$\dot{Q}_V$ /  " + unit)
    ax[0][1].set_xlabel("Time / s")
    ax[0][1].grid()

    # Q_VL plot
    ax[1][0].plot(tank.sol.t, (tank.data['Q_VL'] * unit_conv[unit]), label="Q_VL", color = cmap(1/6))
    ax[1][0].set_ylabel("$\dot{Q}_{VL}$ / " + unit)
    ax[1][0].set_xlabel("Time / s")
    ax[1][0].grid()

    # Q_{V,w} plot
    #ax[1][1].plot(tank.sol.t, (tank.data['Q_Vw'] * unit_conv[unit]), label="Q_Vw",color = cmap(1/6))
    ax[1][1].plot(tank.sol.t, ( (tank.data['Q_Vw'] + tank.data['Q_VL'] + tank.data['Q_L'])  *
                                unit_conv[unit]), label="Q_{tot}",color = cmap(1/6))
    # ax[1][1].set_ylabel("$\dot{Q}_{V,w}$ / " + unit)
    ax[1][1].set_ylabel("$\dot{Q}_{tot}$ / " + unit)
    ax[1][1].set_xlabel("Time / s")
    ax[1][1].grid()

    ax[0][0].set_title('Heat ingresses - '+tank.Geo_l+' geometry')

    # [axis.grid() for axis in ax]
    plt.show()

def plot_l_L(tank):
    cmap = plt.get_cmap('cividis')

    #unit_conv = {'m': 1, 'L': 1e3, 'mL': 1e6}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h

    plt.plot(tank.sol.t, tank.data["z"], color = cmap(1/6))
    plt.grid()
    plt.xlabel('Time / s')
    plt.ylabel('$l_L$ / $m$')
    plt.title('Liquid length over time - '+tank.Geo_l+' geometry')
    plt.show()
    return

def plot_A_T(tank):
    cmap = plt.get_cmap('cividis')

    #unit_conv = {'m': 1, 'L': 1e3, 'mL': 1e6}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h

    plt.plot(tank.sol.t, tank.data["A_T"], color = cmap(1/6))
    plt.grid()
    plt.xlabel('Time / s')
    plt.ylabel('$A_T$ / $m^2$')
    plt.show()
    return

def plot_LF(tank):
    cmap = plt.get_cmap('cividis')

    #unit_conv = {'m': 1, 'L': 1e3, 'mL': 1e6}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h

    plt.plot(tank.sol.t, tank.data["LF"], color = cmap(1/6))
    plt.grid()
    plt.xlabel('Time / s')
    plt.title('Percent liquid filling over time - '+tank.Geo_l+' geometry')
    plt.ylabel('$LF$')
    plt.show()
    return

def plot_rho_V_avg(tank):
    cmap = plt.get_cmap('cividis')

    #unit_conv = {'m': 1, 'L': 1e3, 'mL': 1e6}

    # Extract evaporation and BOG rates and convert to kg/h
    # Visualise evaporation and boil-off gas rate in kg/h

    plt.plot(tank.sol.t, tank.data["rho_V_avg"], color = cmap(1/6))
    plt.grid()
    plt.xlabel('Time / s')
    plt.ylabel(r'$\rho$ / kg/$m^3$')
    plt.title('Average vapour density over time - '+tank.Geo_v+' geometry')
    plt.show()
    return

def plot_vz(tank):
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
        height = np.roots([-np.pi/3,np.pi*tank.l/2,0,-tank.sol.y[0][i*plot_step]])[1]
        tank.z = height

        zed = tank.z_grid*(tank.l-height) + height

        v_z = tank.v_z*(height/zed)*(2*tank.d_i/2 - height)/(2*tank.d_i/2 - zed)
        #end = len(zed) - len(np.where(zed>tank.d_i*0.98)) - 1
        for j, val in np.ndenumerate(zed):
            if val>tank.d_i*0.98:
                v_z[j[0]] = v_z[j[0]-1]
        # Plot the temperature profile at this time step, with color indicating the time
        zed = np.insert(zed,0,0)
        v_z = np.insert(v_z,0,v_z[0])
        ax.plot(v_z, zed, color=cmap(norm(tank.sol.t[i * plot_step])))

        # Add a text box with the time value at the right of the plot
        ax.text(1.02, ((i-1 + 0.15)*plot_step) / len(tank.sol.t), f't={tank.sol.t[i*plot_step]:.0f} s', transform=ax.transAxes, 
                verticalalignment='center', bbox=dict(boxstyle='round,pad=0.5', edgecolor='none', facecolor=cmap(norm(tank.sol.t[i*plot_step])), alpha=0.6))

    # Add a grid
    ax.grid(True)

    # Add labels
    ax.set_ylabel('Tank Length / m')
    ax.set_xlabel('Velocity / m/s')
    ax.set_ylim(0,0.985*tank.d_i)

    ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.1e'))

    # Add a title
    ax.set_title('Velocity profiles at different times')

    # Show the plot
    # plt.show()


# Plots advective velocity profiles considering horizontal geometry

def plot_vz_h(tank): 
    import numpy as np
    from scipy.optimize import root_scalar
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    # Number of profiles
    n_plots = int(tank.sol.t[-1] / tank.plot_interval)
    plot_step = int(tank.plot_interval / tank.time_interval)

    cmap = plt.get_cmap('cividis')
    norm = plt.Normalize(vmin=tank.sol.t[1], vmax=tank.sol.t[-1])

    fig, ax = plt.subplots()

    eps = 1e-9  

    for i in range(1, n_plots + 1):
        V_L = float(tank.sol.y[0][i * plot_step])

        # 1) Altura (o "z" de líquido) según geometría del LÍQUIDO
        if tank.Geo_l == "horizontal":
            # usa tu misma función geométrica
            height = root_scalar(
                tank.calculate_z_height,
                args=(V_L,),
                method='brentq',
                bracket=[1e-9, tank.d_i - 1e-9]).root

        tank.z = height

        zed = tank.z_grid * (tank.l - height) + height
        # recortar un poco bajo d_i para evitar singularidad geométrica
        zcap = np.minimum(zed, 0.98 * tank.d_i)

        if tank.Geo_v == "horizontal":

            A_h = np.sqrt(np.maximum(tank.d_i * height - height**2, 0.0))
            A_z = np.sqrt(np.maximum(tank.d_i * zcap - zcap**2, 0.0))
            v_profile = tank.v_z * (np.maximum(A_h, eps) / np.maximum(A_z, eps))


            #Preparo para graficar
            #z_to_plot = np.insert(zed, 0, 0.0)
            #v_to_plot = np.insert(v_profile, 0, v_profile[0])

            #ax.plot(v_to_plot, z_to_plot, color=cmap(norm(tank.sol.t[i * plot_step])))
            ax.plot(v_profile, zed, color=cmap(norm(tank.sol.t[i * plot_step])))

        ax.text(
            1.02,
            ((i - 1 + 0.15) * plot_step) / len(tank.sol.t),
            f't={tank.sol.t[i * plot_step]:.0f} s',
            transform=ax.transAxes,
            verticalalignment='center',
            bbox=dict(
                boxstyle='round,pad=0.5',
                edgecolor='none',
                facecolor=cmap(norm(tank.sol.t[i * plot_step])),
                alpha=0.6),)

    ax.grid(True)
    ax.set_ylabel('Tank height / m')
    ax.set_xlabel('Velocity / m/s')
    ax.set_ylim(0, 0.985 * tank.d_i)
    ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.1e'))
    ax.set_title('Velocity profiles at different times')
    # plt.show()



def plot_tv_2x1(tank1, tank2, title1=None, title2=None):
    """
    Grafica los perfiles de temperatura de vapor de dos tanques en una figura 2x1,
    manteniendo el mismo formato de plot_tv (curvas coloreadas + cajitas de tiempo).

    Parameters
    ----------
    tank1, tank2 : Tank
        Objetos Tank con la solución ya calculada (tank.sol no None).
    title1, title2 : str, opcional
        Títulos de cada subplot. Si son None, se genera uno usando Geo_v.
    """

    import matplotlib.pyplot as plt
    import numpy as np

    if tank1.sol is None or tank2.sol is None:
        raise TypeError(
            "Ambos tanques deben tener 'sol' definido.\n"
            "Corre tank.evaporate(t_f) antes de usar plot_tv_2x1."
        )

    # Figura 2x1, con espacio a la derecha para las cajitas de tiempo
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 7))
    fig.subplots_adjust(right=0.80, hspace=0.35)

    # ===== Función interna: mismo estilo que plot_tv, pero en un axis dado =====
    def _plot_tv_on_axis(ax, tank, title):
        # Número de perfiles y paso
        n_plots   = int(tank.sol.t[-1] / tank.plot_interval)
        plot_step = int(tank.plot_interval / tank.time_interval)

        # Colormap y normalización como en plot_tv
        cmap = plt.get_cmap('cividis')
        norm = plt.Normalize(vmin=tank.sol.t[1], vmax=tank.sol.t[-1])

        for i in range(1, n_plots + 1):
            idx = i * plot_step
            t_i = tank.sol.t[idx]
            T_v = tank.sol.y[1:, idx]

            # Curva Tv(z)
            ax.plot(T_v, tank.z_grid, color=cmap(norm(t_i)))

            # Posición vertical normalizada para las cajitas de tiempo
            y_pos = (i - 0.5) / n_plots

            ax.text(
                1.02, y_pos,
                f"t={t_i:.0f} s",
                transform=ax.transAxes,
                verticalalignment="center",
                bbox=dict(
                    boxstyle="round,pad=0.5",
                    edgecolor="none",
                    facecolor=cmap(norm(t_i)),
                    alpha=0.6,
                ),
            )

        ax.grid(True)
        ax.set_ylabel(r"Dimensionless length $\zeta = z/l_V$")
        if title is None:
            title = (
                "Vapour temperature profiles at different times - "
                + tank.Geo_v
                + " geometry"
            )
        ax.set_title(title)

    # ---- Subplot superior: tank1 ----
    _plot_tv_on_axis(ax1, tank1, title1)

    # ---- Subplot inferior: tank2 ----
    _plot_tv_on_axis(ax2, tank2, title2)

    # Eje x común
    ax2.set_xlabel("Temperature / K")

    # Título general
    fig.suptitle("Vapour temperature profiles – comparison of two tanks", y=0.98)

    return fig, (ax1, ax2)


def plot_V_L_2x1(tank1, tank2, unit='m3', title1=None, title2=None):
    """
    Figura 2x1 con V_L(t) de dos tanques.
    Copia EXACTAMENTE el estilo de plot_V_L.
    """


    cmap = plt.get_cmap('cividis')
    unit_conv = {'m3': 1, 'L': 1e3, 'mL': 1e6}

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    fig.subplots_adjust(hspace=0.35)

    # --- Upper plot ---
    ax1.plot(tank1.sol.t, tank1.sol.y[0]*unit_conv[unit], color=cmap(1/6))
    ax1.grid()
    ax1.set_ylabel(f"$V_L$ / {unit}")
    if title1 is None:
        title1 = f"Liquid volume over time – {tank1.Geo_l} geometry"
    ax1.set_title(title1)

    # --- Lower plot ---
    ax2.plot(tank2.sol.t, tank2.sol.y[0]*unit_conv[unit], color=cmap(1/6))
    ax2.grid()
    ax2.set_ylabel(f"$V_L$ / {unit}")
    if title2 is None:
        title2 = f"Liquid volume over time – {tank2.Geo_l} geometry"
    ax2.set_title(title2)
    ax2.set_xlabel("Time / s")

    fig.suptitle("Liquid Volume Comparison – 2 Tanks", y=0.98)
    return fig, (ax1, ax2)


def plot_BOG_2x1(tank1, tank2, unit='kg/h', title1=None, title2=None):
    """
    Figura 2x1 con BOG(t) y B_L(t) para dos tanques.
    Copia EXACTAMENTE el estilo de plot_BOG.
    """
    cmap = plt.get_cmap('cividis')

    unit_conv = {'kg/h': 3600, 'kg/s':1, 'g/h':3600*1e3, 'g/s':1000}

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    fig.subplots_adjust(hspace=0.35)

    # -------- Tank 1 --------
    ax1.plot(tank1.sol.t, tank1.data['B_L']*unit_conv[unit],
             label='Evaporation rate, $\dot{B}_L$', color=cmap(1/6))
    ax1.plot(tank1.sol.t[1:], tank1.data['BOG'][1:]*unit_conv[unit],
             label='Boil-off gas rate, $\dot{B}$', color=cmap(5/6))
    ax1.grid()
    ax1.set_ylabel(f"Mass flow / {unit}")
    if title1 is None:
        title1 = f"BOG and evaporation rates – {tank1.Geo_l} geometry"
    ax1.set_title(title1)
    ax1.legend()

    # -------- Tank 2 --------
    ax2.plot(tank2.sol.t, tank2.data['B_L']*unit_conv[unit],
             label='Evaporation rate, $\dot{B}_L$', color=cmap(1/6))
    ax2.plot(tank2.sol.t[1:], tank2.data['BOG'][1:]*unit_conv[unit],
             label='Boil-off gas rate, $\dot{B}$', color=cmap(5/6))
    ax2.grid()
    ax2.set_ylabel(f"Mass flow / {unit}")
    if title2 is None:
        title2 = f"BOG and evaporation rates – {tank2.Geo_l} geometry"
    ax2.set_title(title2)
    ax2.set_xlabel("Time / s")
    ax2.legend()

    fig.suptitle("BOG Comparison – 2 Tanks", y=0.98)
    return fig, (ax1, ax2)

def plot_Q_2x1(tank1, tank2, unit='kW', title1=None, title2=None):
    """
    Figura 2x1, cada subplot contiene un layout 2x2 con:
    Q_L, Q_V, Q_VL, Q_tot
    Copia EXACTAMENTE el estilo de plot_Q.
    """

    cmap = plt.get_cmap('cividis')

    unit_conv = {'W': 1, 'kW': 1e-3}

    fig = plt.figure(figsize=(8, 10))

    # ======= Tank 1 (top 2x2) =======
    ax11 = fig.add_subplot(4, 2, 1)
    ax12 = fig.add_subplot(4, 2, 2)
    ax13 = fig.add_subplot(4, 2, 3)
    ax14 = fig.add_subplot(4, 2, 4)

    ax11.plot(tank1.sol.t, tank1.data['Q_L']*unit_conv[unit], color=cmap(1/6))
    ax12.plot(tank1.sol.t, tank1.data['Q_V']*unit_conv[unit], color=cmap(1/6))
    ax13.plot(tank1.sol.t, tank1.data['Q_VL']*unit_conv[unit], color=cmap(1/6))
    ax14.plot(tank1.sol.t, ((tank1.data['Q_Vw'] + tank1.data['Q_VL'] + tank1.data['Q_L'])*unit_conv[unit]), color=cmap(1/6))

    for ax in [ax11, ax12, ax13, ax14]:
        ax.grid()

    ax11.set_ylabel("$\dot{Q}_L$ / "+unit)
    ax12.set_ylabel("$\dot{Q}_V$ / "+unit)
    ax13.set_ylabel("$\dot{Q}_{VL}$ / "+unit)
    ax14.set_ylabel("$\dot{Q}_{tot}$ / "+unit)
    ax13.set_xlabel("Time / s")
    ax14.set_xlabel("Time / s")

    if title1 is None:
        title1 = f"Heat ingresses – {tank1.Geo_l} geometry"
    ax11.set_title(title1)

    # ======= Tank 2 (bottom 2x2) =======
    ax21 = fig.add_subplot(4, 2, 5)
    ax22 = fig.add_subplot(4, 2, 6)
    ax23 = fig.add_subplot(4, 2, 7)
    ax24 = fig.add_subplot(4, 2, 8)

    ax21.plot(tank2.sol.t, tank2.data['Q_L']*unit_conv[unit], color=cmap(1/6))
    ax22.plot(tank2.sol.t, tank2.data['Q_V']*unit_conv[unit], color=cmap(1/6))
    ax23.plot(tank2.sol.t, tank2.data['Q_VL']*unit_conv[unit], color=cmap(1/6))
    ax24.plot(tank2.sol.t, ((tank2.data['Q_Vw'] + tank2.data['Q_VL'] + tank2.data['Q_L'])*unit_conv[unit]), color=cmap(1/6))

    for ax in [ax21, ax22, ax23, ax24]:
        ax.grid()

    ax23.set_xlabel("Time / s")
    ax24.set_xlabel("Time / s")

    if title2 is None:
        title2 = f"Heat ingresses – {tank2.Geo_l} geometry"
    ax21.set_title(title2)

    fig.suptitle("Heat Ingress Comparison – 2 Tanks", y=0.99)
    fig.tight_layout()
    return fig