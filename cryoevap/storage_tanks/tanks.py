import numpy as np
from .tank import Tank

# Tank properties
def create_small(LF=0.97):
    Q_roof = 0  # [W]
    d_i = 1.604  # [m]
    d_o = 1.630  # [m]
    V_tank = 8  # [m^3]
    T_air = 288.15  # [K]
    U_V = 0.019  # W/[m^2*K]
    small_tank = Tank(d_i, d_o, V_tank, LF)
    small_tank.Q_VL = 0
    small_tank.set_HeatTransProps(U_V, U_V, Q_roof=Q_roof, T_air=T_air)
    return small_tank


def create_vsmall(LF=0.97):
    Q_roof = 0  # [W]
    d_i = 1.604 / 10  # [m]
    d_o = 1.630 / 10  # [m]
    V_tank = 0.08 * 0.06  # [m^3]
    T_air = 288.15  # [K]
    U_V = 0.019  # W/[m^2*K]
    tank = Tank(d_i, d_o, V_tank, LF)
    tank.Q_VL = 0
    tank.set_HeatTransProps(U_V, U_V, Q_roof=Q_roof, T_air=T_air)
    return tank


def create_large(LF=0.97):
    # Input tank properties
    Q_roof = 0  # W
    d_i = 76.4  # m
    d_o = 80  # m
    V_tank = 165000  # m^3
    T_air = 293.15  # K
    U_V = 0.037  # W/m^2K
    large_tank = Tank(d_i, d_o, V_tank, LF)
    Q_b = 60000  # W, heat ingress from the bottom
    return large_tank


def create_80m3(LF=0.97):
    # Input tank properties
    Q_roof = 0  # W
    d_i = 2.85  # m
    d_o = 3.15  # m
    V_tank = 80  # m^3
    T_air = 288.15  # K
    U_V = 0.011  # W/m^2K
    tank = Tank(d_i, d_o, V_tank, LF)
    tank.Q_VL = 0
    Q_b = 0  # W, heat ingress from the bottom
    tank.set_HeatTransProps(U_V, U_V, Q_roof=Q_roof, T_air=T_air)
    return tank


def create_test_tube(LF=0.97):
    Q_roof = 0  # [W]
    d_i = 1.604 / 25  # [m]
    d_o = 1.630 / 25  # [m]
    V_tank = np.pi * d_i ** 2 / 4 * 0.1583676  # [m^3]
    T_air = 288.15  # [K]
    U_V = 0.019  # W/[m^2*K]
    small_tank = Tank(d_i, d_o, V_tank, LF)
    small_tank.Q_VL = 0
    small_tank.set_HeatTransProps(U_V, U_V, Q_roof=Q_roof, T_air=T_air)
    return small_tank


def create_seo(LF=0.5):
    Q_roof = 0  # [W]
    V_tank = 6.75e-3  # m^3
    height = 0.213  # m
    d_i = np.sqrt(4 * 6.75e-3 / (np.pi * height))
    d_o = 0.201  # [m]

    T_air = 288.15  # [K]
    U_V = 0.042379  # W/[m^2*K]
    tank = Tank(d_i, d_o, V_tank, LF)
    tank.Q_VL = 0
    tank.set_HeatTransProps(U_V, U_V, Q_roof=Q_roof, T_air=T_air)
    return tank
