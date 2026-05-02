"""
Fig_Phase1.py
=============
Compares the LNH3 isobaric evaporation model with and without the
Phase-1 lumped metallic roof-shell ODE (T_roof state variable).

Physical model (Phase 1):
  - Lumped 0-D ODE for the metallic roof panel temperature T_roof.
  - Energy balance: external convection from T_env heats the roof;
    the roof drives the vapour Robin BC at z = z_R via U_roof.
  - Replaces the default Neumann (insulated) roof BC, which allowed
    unphysical unbounded vapour temperature growth.

Figures produced
----------------
Figures/Fig_Phase1.svg  – 3-panel comparison:
    (a) Vapour temperature profiles at selected snapshots
        (dashed = no roof model / Neumann, solid = Phase 1 roof shell)
    (b) BOG rate vs time
    (c) Roof shell temperature T_roof and top vapour node T_v[-1] vs time

Run from the Results/ directory:
    python Fig_Phase1.py
"""

import os, sys, pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.ndimage import uniform_filter1d

# Resolve the project root so the script can be executed from any cwd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from cryoevap.storage_tanks import Tank
from cryoevap.cryogens import Cryogen

# ---------------------------------------------------------------------------
# Global font sizes (match existing Results/*.py convention)
# ---------------------------------------------------------------------------
FS_LABEL  = 14
FS_TICKS  = 13
FS_LEGEND = 12

# ---------------------------------------------------------------------------
# Tank and cryogen parameters  (LNH3 large-scale, 5059 m³ at LF = 0.95)
# ---------------------------------------------------------------------------
e_w   = 18 * 0.0254           # Insulation thickness / m  (18 in perlite annulus)
d_i   = 10.24 * 2             # Internal diameter / m
d_o   = d_i + 2 * e_w         # External diameter / m
V_tank = 5059.889             # Volume / m³
LF    = 0.95
P     = 101325                # Operating pressure / Pa

T_air = 5.3286 + 273.15       # Mean environmental temperature / K
T_range = 15                  # Daily temperature swing / K
h_env = 396                   # External convection coefficient / W m⁻² K⁻¹

U_L = 0.087                   # Overall HTC (liquid side) / W m⁻² K⁻¹
U_V = 0.087                   # Overall HTC (vapour side) / W m⁻² K⁻¹
h_L = 124.96                  # Internal liquid convection coefficient / W m⁻² K⁻¹
eta_w = 0.70

# Insulation / perlite-filled annular wall (outer jacket)
k_ins  = 0.0411               # W m⁻¹ K⁻¹
rho_ins = 60                  # kg m⁻³
cp_ins  = 900                 # J kg⁻¹ K⁻¹

# Inner metallic wall (C-Mn steel)  [Phase 1]
e_wi  = 0.020                 # 20 mm inner steel shell / m
k_wi  = 50.0                  # W m⁻¹ K⁻¹
rho_wi = 7850.0               # kg m⁻³
cp_wi  = 490.0                # J kg⁻¹ K⁻¹

# Simulation duration and grid
EVAP_TIME    = 48 * 3600      # 48 h in seconds
TIME_INTERVAL = 0.5 * 3600   # recording interval
dz = 0.05                     # axial spacing / m
dr = 0.01                     # radial spacing / m

# Pickle cache paths (relative to Results/)
CACHE_DIR  = "Data"
PKL_NOWALL = os.path.join(CACHE_DIR, "phase1_nowall.pkl")
PKL_WALL   = os.path.join(CACHE_DIR, "phase1_wall.pkl")

# ---------------------------------------------------------------------------
# Helper – build and run one tank case
# ---------------------------------------------------------------------------

def build_tank(use_inner_wall: bool) -> Tank:
    tank = Tank(d_i, d_o, V_tank, LF)
    tank.set_HeatTransProps(
        U_L, U_V, T_air,
        Q_b_fixed=None, Q_roof=0,
        eta_w=eta_w,
        k_w=k_ins, rho_w=rho_ins, cp_w=cp_ins,
        h_L=h_L, T_init=True,
        k_wi=k_wi, rho_wi=rho_wi, cp_wi=cp_wi,
        e_wi=(e_wi if use_inner_wall else 0.0),
        h_env_roof=h_env,
    )
    nh3 = Cryogen(name="ammonia")
    nh3.set_coolprops(P)
    tank.cryogen = nh3

    # Grid
    l_V = tank.l * (1 - LF)
    e_ins = (d_o - d_i) / 2
    n_z = max(5, 1 + int(round(l_V / dz)))
    n_r = max(3, 1 + int(round(e_ins / dr)))
    tank.z_grid = np.linspace(0, 1, n_z)
    tank.r_grid = np.linspace(0, 1, n_r)

    # Environmental properties (daily cycle, no annual variation)
    tank.set_EnvironmentalProps(T_avg_day=T_air, T_range_day=T_range, h_env=h_env)

    tank.time_interval = TIME_INTERVAL
    return tank


def run_case(use_inner_wall: bool, cache_path: str) -> dict:
    """Run simulation (or load from cache) and return the data dict + extras."""
    if os.path.exists(cache_path):
        print(f"  Loading cached data from {cache_path}")
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    label = "Phase 1 (roof shell)" if use_inner_wall else "No roof model"
    print(f"  Running: {label} ...")
    tank = build_tank(use_inner_wall)
    tank.evaporate(EVAP_TIME)

    payload = {
        "data":   tank.data,
        "z_grid": tank.z_grid,
        "r_grid": tank.r_grid,
        "l":      tank.l,
        "LF":     LF,
        "e_wi":   tank.e_wi,
    }
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(cache_path, "wb") as f:
        pickle.dump(payload, f)
    print(f"    Saved → {cache_path}")
    return payload


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    os.makedirs("Figures", exist_ok=True)

    print("=== Phase 1 comparison ===")
    case_no  = run_case(use_inner_wall=False, cache_path=PKL_NOWALL)
    case_yes = run_case(use_inner_wall=True,  cache_path=PKL_WALL)

    d_no  = case_no["data"]
    d_yes = case_yes["data"]
    z_grid = case_yes["z_grid"]   # same grid for both cases

    t_h_no  = d_no["Time"]  / 3600
    t_h_yes = d_yes["Time"] / 3600

    # -----------------------------------------------------------------------
    # Colour map and snapshot times for the Tv profile panel
    # -----------------------------------------------------------------------
    N_SNAPS = 5
    cmap_tv = plt.get_cmap("cividis", N_SNAPS + 1)

    snap_times = np.linspace(0, EVAP_TIME, N_SNAPS)   # seconds

    def get_snap(data_dict, t_sec, z_grid_local):
        """Return (T_v profile, actual time index, time in h) for a snapshot."""
        idx = int(np.argmin(np.abs(data_dict["Time"] - t_sec)))
        n_z = len(z_grid_local)
        T_v = data_dict["T_V_raw"][:, idx]
        return T_v, idx, data_dict["Time"][idx] / 3600

    # -----------------------------------------------------------------------
    # BOG smoothing
    # -----------------------------------------------------------------------
    SMOOTH_W = 5   # uniform filter half-width in sample points

    def smooth(arr):
        return uniform_filter1d(arr, size=SMOOTH_W)

    bog_no  = smooth(d_no["BOG"]  * 3600)   # kg/h
    bog_yes = smooth(d_yes["BOG"] * 3600)

    # -----------------------------------------------------------------------
    # Figure layout: 3 panels
    # -----------------------------------------------------------------------
    fig = plt.figure(figsize=(16, 5), dpi=200)
    ax1 = fig.add_subplot(1, 3, 1)   # (a) Tv profiles
    ax2 = fig.add_subplot(1, 3, 2)   # (b) BOG rate
    ax3 = fig.add_subplot(1, 3, 3)   # (c) T_wi

    # ── Panel (a): Vapour temperature profiles ──────────────────────────────
    for k, t_sec in enumerate(snap_times):
        color = cmap_tv(k)

        T_v_no,  _, t_h = get_snap(d_no,  t_sec, z_grid)
        T_v_yes, _, _   = get_snap(d_yes, t_sec, z_grid)

        ax1.plot(T_v_no,  z_grid, color=color, linestyle="--", linewidth=1.5)
        ax1.plot(T_v_yes, z_grid, color=color, linestyle="-",  linewidth=1.5,
                 label=f"$t = {t_h:.0f}$ h")

    ax1.set_xlabel("Vapour temperature / K", fontsize=FS_LABEL)
    ax1.set_ylabel(r"Dimensionless height $\zeta = z/l_V$", fontsize=FS_LABEL)
    ax1.tick_params(labelsize=FS_TICKS)
    ax1.grid(True, alpha=0.4)

    # Legend for the snapshots
    snap_handles = [
        mlines.Line2D([], [], color=cmap_tv(k), linewidth=1.5,
                      label=f"$t = {snap_times[k]/3600:.0f}$ h")
        for k in range(N_SNAPS)
    ]
    style_handles = [
        mlines.Line2D([], [], color="k", linestyle="--", linewidth=1.5,
                      label="No roof model (Neumann)"),
        mlines.Line2D([], [], color="k", linestyle="-",  linewidth=1.5,
                      label="Phase 1: lumped roof shell"),
        mlines.Line2D([], [], color="grey", linestyle="--", linewidth=1.2,
                      label=r"$T_{env}$ (mean)"),
        mlines.Line2D([], [], color="grey", linestyle=":",  linewidth=1.2,
                      label=r"$T_{sat}$ (NH$_3$, 1 atm)"),
    ]
    ax1.legend(handles=snap_handles + style_handles,
               fontsize=FS_LEGEND - 1, loc="upper left",
               framealpha=0.85, ncol=1)
    # Reference: T_env (mean) — vapour should be bounded below this
    ax1.axvline(T_air, color="grey", linestyle="--", linewidth=1.2, alpha=0.7,
                label=r"$T_{env}$ (mean)")
    ax1.axvline(239.72, color="grey", linestyle=":", linewidth=1.2, alpha=0.7,
                label=r"$T_{sat}$ (NH$_3$, 1 atm)")

    ax1.set_title("(a) Vapour temperature profiles", fontsize=FS_LABEL)

    # ── Panel (b): BOG rate ─────────────────────────────────────────────────
    cmap_bog = plt.get_cmap("inferno", 4)

    ax2.plot(t_h_no,  bog_no,  color=cmap_bog(1), linestyle="--",
             linewidth=1.8, label="No roof model (Neumann)")
    ax2.plot(t_h_yes, bog_yes, color=cmap_bog(2), linestyle="-",
             linewidth=1.8, label="Phase 1: lumped roof shell")

    ax2.set_xlabel("Time / h", fontsize=FS_LABEL)
    ax2.set_ylabel(r"Boil-off gas rate / $\mathrm{kg\ h^{-1}}$", fontsize=FS_LABEL)
    ax2.tick_params(labelsize=FS_TICKS)
    ax2.legend(fontsize=FS_LEGEND, loc="upper right", framealpha=0.85)
    ax2.grid(True, alpha=0.4)
    ax2.set_title("(b) Boil-off gas rate", fontsize=FS_LABEL)

    # ── Panel (c): Roof shell T_roof and vapour top T_v[-1] ─────────────────
    cmap_wi = plt.get_cmap("inferno", 5)

    # Reference lines
    T_sat_nh3 = 239.72   # K  (1 atm)
    ax3.axhline(T_sat_nh3, color="grey", linestyle=":", linewidth=1.5,
                label=r"$T_{sat}$ (NH$_3$, 1 atm)")
    ax3.axhline(T_air, color="grey", linestyle="--", linewidth=1.5,
                label=r"$T_{env}$ (mean)")

    if case_yes["e_wi"] > 0 and "T_wi_raw" in d_yes:
        # Lumped roof shell temperature (new Phase 1 state variable)
        T_roof_arr = d_yes["T_wi_raw"][0, :]   # shape (n_t,)
        ax3.plot(t_h_yes, T_roof_arr,
                 color=cmap_wi(3), linewidth=2.0, label=r"$T_{roof}$ (roof shell, Phase 1)")

        # Top vapour node T_v[-1] — shows coupling between roof and vapour
        n_z_yes = len(case_yes["z_grid"])
        T_v_top_yes = d_yes["T_V_raw"][-1, :]   # last row = z_R node
        ax3.plot(t_h_yes, T_v_top_yes,
                 color=cmap_wi(2), linewidth=1.5, linestyle="-.",
                 label=r"$T_v(z_R)$ (vapour at roof, Phase 1)")

    # Top vapour node for the no-roof case (Neumann) — shows unphysical heating when missing
    n_z_no = len(case_no["z_grid"])
    T_v_top_no = d_no["T_V_raw"][-1, :]
    ax3.plot(t_h_no, T_v_top_no,
             color=cmap_wi(1), linewidth=1.5, linestyle="--",
             label=r"$T_v(z_R)$ (vapour at roof, no model)")

    ax3.set_xlabel("Time / h", fontsize=FS_LABEL)
    ax3.set_ylabel("Temperature / K", fontsize=FS_LABEL)
    ax3.tick_params(labelsize=FS_TICKS)
    ax3.legend(fontsize=FS_LEGEND - 1, loc="upper right", framealpha=0.85)
    ax3.grid(True, alpha=0.4)
    ax3.set_title(r"(c) Roof shell and vapour top temperature", fontsize=FS_LABEL)

    # ── Save ────────────────────────────────────────────────────────────────
    plt.tight_layout()
    out_path = "Figures/Fig_Phase1.svg"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"\nSaved → {out_path}")
    plt.show()
