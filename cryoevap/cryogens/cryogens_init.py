# Initialize Methane and Nitrogen
from .cryogen import Cryogen


def methane_init():
    """Load methane thermopysical properites at 
    its saturation point at P = 100000PA Pa"""
    methane = Cryogen(
        "methane",
        100000,
        111.508,
        2.6205e4,
        128.1975,
        98.9155,
        8.2684e3,
        0.0117,
        35.7950,
        16.0428,
    )
    return methane


def nitrogen_init():
    """Load methane thermopysical properites at 
    its saturation point at P = 100000 Pa"""
    k_V = 7.1744e-3  # W/(m2K)
    rho_V = 162.65316  # mol/m^3
    rho_L = 2.8793e4  # mol/m^3
    cp_V = 31.4624  # J/molK
    T_L = 77.2435  # /K
    h_V = 2.2045e3  # J/kgK
    h_L = -3.3132e3  # J/kgK
    P = 100000  # Pa
    MW = 28.0134  # g/mol
    nitrogen = Cryogen("nitrogen", P, T_L, rho_L, rho_V, h_L, h_V, k_V, cp_V, MW)
    return nitrogen
