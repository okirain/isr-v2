"""Comparison of the pumped NaK-78 loop and the passive sodium heat-pipe array.

The comparison is deliberately transparent about uncertainty.  The pumped-loop
mass is a component build-up with published or standard engineering estimates;
the heat-pipe mass comes from the geometry of the reference array and the
reflector/block mass penalty of the larger core.  Both are reported with a
contingency so that the conclusion is not an artefact of a single assumption.

Run:

    python -m src.isr_design.comparison
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from ..solver.limits import HeatPipeGeometry, limit_curve, operating_limits
from ..solver.transient_model import (
    HeatPipeTransientModel,
    TransientConfig,
    step_profile,
)
from ..thermophysical.properties import POTASSIUM, SODIUM
from .isr_core_coupling import (
    ArrayDesign,
    N_PINS,
    P_TOTAL,
    REACTOR_MASS,
    T_HOT_PIPE,
    UN_MELTING_LIMIT,
    design_array,
    evaluate_reference,
    fuel_temperature_profile,
    in_channel_case,
    make_geometry,
    scan_designs,
    select_reference,
)

RESULTS_DIR = Path(__file__).resolve().parents[2] / "results"
FIG_DIR = Path(__file__).resolve().parents[2] / "thesis_md" / "figures"


# ---------------------------------------------------------------------------
# Mass models
# ---------------------------------------------------------------------------
def pumped_loop_mass() -> pd.DataFrame:
    """Component build-up of the external NaK-78 coolant loop, kg.

    The in-pin rhenium channels (0.23 kg per pin) are already part of the
    Khandaq reactor mass and are not counted again.  The remaining items are
    engineering estimates for a 500 kWth space loop with a 6 kW(e) EM pump.
    """
    components = [
        ("External piping (SS316, 2 x 3 m)", 22.0, "3 m supply + 3 m return, 5 cm OD"),
        ("EM pump and controller", 15.0, "6 kW(e) input, 10 % efficiency"),
        ("Primary-to-secondary heat exchanger", 20.0, "compact plate-fin, 500 kWth"),
        ("Expansion / accumulator vessel", 8.0, "NaK-78, 1 bar cover gas"),
        ("NaK-78 inventory", 4.0, "external loop, excluding pin channels"),
        ("Insulation, supports and instrumentation", 10.0, "MLI, brackets, flowmeters"),
    ]
    df = pd.DataFrame(components, columns=["component", "mass_kg", "basis"])
    subtotal = df["mass_kg"].sum()
    df.loc[len(df)] = ["Contingency (15 %)", 0.15 * subtotal, "design maturity"]
    return df


def heat_pipe_system_mass(design: ArrayDesign) -> pd.DataFrame:
    """Component build-up of the passive heat-pipe system, kg."""
    components = [
        ("Sodium heat pipes (19 x 3.0 cm)", design.hp_mass_total,
         "Haynes 230 shell, Ni wick, Na inventory"),
        ("Headers and manifold", design.header_mass, "collector/condenser interface"),
        ("Core-block mass penalty", design.block_mass_added,
         "larger Mo-30Re block, 34 cm core"),
        ("Reflector mass penalty", design.reflector_mass_added,
         "BeO reflector scaled with core diameter"),
    ]
    df = pd.DataFrame(components, columns=["component", "mass_kg", "basis"])
    subtotal = df["mass_kg"].sum()
    df.loc[len(df)] = ["Contingency (10 %)", 0.10 * subtotal, "advanced wick"]
    return df


def mass_comparison(design: ArrayDesign) -> pd.DataFrame:
    pumped = pumped_loop_mass()["mass_kg"].sum()
    hp = heat_pipe_system_mass(design)["mass_kg"].sum()
    return pd.DataFrame({
        "architecture": ["Pumped NaK-78", "Passive Na heat pipe"],
        "loop_mass_kg": [pumped, hp],
        "reactor_mass_kg": [REACTOR_MASS, REACTOR_MASS],
        "total_mass_kg": [REACTOR_MASS + pumped, REACTOR_MASS + hp],
        "pump_power_kw_e": [6.0, 0.0],
    })


# ---------------------------------------------------------------------------
# Thermo-mechanical comparison
# ---------------------------------------------------------------------------
def fuel_profiles(design: ArrayDesign, T_hp: float = T_HOT_PIPE) -> dict:
    """Radial fuel temperature profiles for the two cooling architectures."""
    q_pin = P_TOTAL / N_PINS
    base = fuel_temperature_profile(q_pin, 850.0, boundary="inner")
    ref = evaluate_reference(design, T_hp=T_hp)
    v2 = ref.fuel
    return {
        "baseline": base,
        "v2": v2,
        "T_block": ref.T_block,
        "q_flux": ref.q_flux_block,
        "block_delta_T": ref.block_delta_T,
    }


def single_pipe_startup(
    design: ArrayDesign,
    t_end_h: float = 12.0,
    dt: float = 20.0,
) -> object:
    """Transient startup of one representative ISR heat pipe (26.3 kW)."""
    q_pipe = design.q_per_hp
    cfg = TransientConfig(
        geometry=design.geometry,
        n_pipes=1,
        n_axial=30,
        fuel_mass_total=155.65 / design.n_hp,
        fuel_cp=200.0,
        sink_mass_total=30.0,
        sink_cp=500.0,
        h_fuel_wall=500.0,
        h_vapor_evap=500.0,
        h_vapor_cond=500.0,
        h_sink=100.0,
        eps_area_mli=0.02,
        ua_reject=60.0,
        gravity=0.0,
        enable_limits=True,
        startup_temperature=700.0,
        startup_width=20.0,
        residual_vapor_conductance=1.0,
        pcs_power_gain=0.0,
    )
    # Controlled 4-hour power ramp, then full power
    power = step_profile([(0.0, 0.0), (4.0 * 3600.0, q_pipe),
                          (t_end_h * 3600.0, q_pipe)], ramp=4.0 * 3600.0)
    return HeatPipeTransientModel(cfg).simulate(
        t_end_h * 3600.0, power, dt=dt, sink_func=lambda t: 1.0
    )


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })
    return plt


def make_limit_figure(design: ArrayDesign) -> Path:
    plt = _style()
    T = np.linspace(750.0, 1300.0, 220)
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.1), constrained_layout=True)

    ax = axes[0]
    curves = limit_curve(T, design.geometry, SODIUM, g=0.0)
    for key, label, color in [
        ("capillary", "Capillary", "#1f4e79"),
        ("sonic", "Sonic", "#c00000"),
        ("entrainment", "Entrainment", "#2e7d32"),
        ("viscous", "Viscous", "#7b1fa2"),
        ("boiling", "Boiling", "#ef6c00"),
    ]:
        ax.semilogy(T, np.maximum(curves[key], 1.0), lw=1.0, color=color,
                    label=label)
    ax.semilogy(T, np.maximum(curves["q_max"], 1.0), lw=1.8, color="black",
                ls="--", label="Binding limit")
    ax.axhline(design.q_per_hp, color="0.35", lw=0.8, ls=":")
    ax.text(760, design.q_per_hp * 1.15, "26.3 kW/pipe", fontsize=7,
            color="0.35")
    ax.set_xlabel("Vapour temperature (K)")
    ax.set_ylabel("Heat transport (W)")
    ax.set_title("(a) Sodium, 3.0 cm arterial heat pipe", fontsize=8.5)
    ax.set_xlim(750, 1300)
    ax.set_ylim(1e2, 1e7)
    ax.grid(alpha=0.25, lw=0.5, which="both")
    ax.legend(loc="lower right", frameon=False, ncol=2, fontsize=6.3)

    ax = axes[1]
    for fluid, label, color in [(SODIUM, "Sodium", "#1f4e79"),
                                (POTASSIUM, "Potassium", "#c00000")]:
        c = limit_curve(T, design.geometry, fluid, g=0.0)
        ax.semilogy(T, np.maximum(c["q_max"], 1.0), lw=1.3, color=color,
                    label=label)
    ax.axhline(design.q_per_hp, color="0.35", lw=0.8, ls=":")
    ax.set_xlabel("Vapour temperature (K)")
    ax.set_ylabel("Maximum transport (W)")
    ax.set_title("(b) Working-fluid comparison", fontsize=8.5)
    ax.set_xlim(750, 1300)
    ax.set_ylim(1e2, 1e7)
    ax.grid(alpha=0.25, lw=0.5, which="both")
    ax.legend(loc="lower right", frameon=False)

    path = FIG_DIR / "fig_5_2_operating_limits.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def make_fuel_figure(design: ArrayDesign) -> Path:
    plt = _style()
    prof = fuel_profiles(design)
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.1), constrained_layout=True)

    ax = axes[0]
    ax.plot(prof["baseline"]["r"] * 100, prof["baseline"]["T"] - 273.15,
            lw=1.4, color="#1f4e79", label="Pumped NaK-78 (inner channel)")
    ax.plot(prof["v2"]["r"] * 100, prof["v2"]["T"] - 273.15,
            lw=1.4, color="#c00000", label="Passive Na heat pipe (outer block)")
    ax.axhline(UN_MELTING_LIMIT - 273.15, color="0.4", lw=0.8, ls="--")
    ax.text(0.5, 2300, "UN limit 2523 K", fontsize=7, color="0.4")
    ax.set_xlabel("Radial position (cm)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("(a) Radial fuel temperature", fontsize=8.5)
    ax.set_xlim(0.3, 1.7)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower left", frameon=False, fontsize=6.8)

    ax = axes[1]
    cases = [0, 1, 3]
    labels = ["Nominal", "1 pipe failed", "3 adjacent failed"]
    q = [design.q_per_hp / 1000.0 * (1 + n / 6) for n in cases]
    Tmax = [prof["v2"]["T_max"]]
    # recompute failure temperatures consistently with evaluate_reference
    ref = evaluate_reference(design)
    Tmax = ref.failure_cases["T_fuel_max_k"].tolist()
    bars = ax.bar(labels, Tmax, color=["#1f4e79", "#e6a23c", "#c00000"], width=0.55)
    ax.axhline(UN_MELTING_LIMIT, color="0.3", lw=0.9, ls="--")
    ax.set_ylabel("Peak fuel temperature (K)")
    ax.set_title("(b) Failure tolerance of the array", fontsize=8.5)
    ax.set_ylim(0, 2900)
    for b, v in zip(bars, Tmax):
        ax.text(b.get_x() + b.get_width() / 2, v + 40, f"{v:.0f} K",
                ha="center", fontsize=7)
    ax.grid(alpha=0.25, lw=0.5, axis="y")

    path = FIG_DIR / "fig_5_3_isr_fuel_temperature.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def make_startup_figure(design: ArrayDesign) -> Path:
    plt = _style()
    res = single_pipe_startup(design)
    fig, ax1 = plt.subplots(figsize=(4.6, 3.0), constrained_layout=True)
    t = res.t / 3600.0
    ax1.plot(t, res.T_fuel - 273.15, lw=1.3, color="#1f4e79", label="Fuel")
    ax1.plot(t, res.T_vapor - 273.15, lw=1.2, color="#c00000", ls="--",
             label="Vapour")
    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel("Temperature (°C)")
    ax1.set_ylim(0, 1200)
    ax1.grid(alpha=0.25, lw=0.5)
    ax2 = ax1.twinx()
    ax2.plot(t, res.Q_hp / 1000.0, lw=0.9, color="#2e7d32", alpha=0.75,
             label="Heat transported")
    ax2.set_ylabel("Heat transport (kW)", color="#2e7d32")
    ax2.tick_params(axis="y", colors="#2e7d32")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="lower right", frameon=False, fontsize=7)
    ax1.set_title("Startup of one 3.0 cm ISR-V2 heat pipe", fontsize=8.5)
    path = FIG_DIR / "fig_5_4_isr_startup.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def make_mass_figure(design: ArrayDesign) -> Path:
    plt = _style()
    pumped = pumped_loop_mass()
    hp = heat_pipe_system_mass(design)
    fig, ax = plt.subplots(figsize=(4.8, 3.0), constrained_layout=True)
    ax.barh(pumped["component"], pumped["mass_kg"], color="#1f4e79",
            alpha=0.85, label="Pumped NaK-78 loop")
    ax.barh(hp["component"], hp["mass_kg"], color="#c00000",
            alpha=0.85, label="Passive Na heat-pipe system")
    ax.set_xlabel("Mass (kg)")
    ax.grid(alpha=0.25, lw=0.5, axis="x")
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    ax.set_title("Balance-of-plant mass comparison", fontsize=8.5)
    path = FIG_DIR / "fig_5_5_mass_comparison.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path




def make_design_scan_figure(scan: pd.DataFrame, design: ArrayDesign) -> Path:
    """Mass and feasibility map of the heat-pipe array design space."""
    plt = _style()
    fig, ax = plt.subplots(figsize=(4.8, 3.2), constrained_layout=True)
    colors = {13: "#7b1fa2", 19: "#1f4e79", 37: "#2e7d32", 61: "#c00000"}
    for n, group in scan.groupby("n_hp"):
        group = group.sort_values("d_hp_cm")
        ax.plot(group["d_hp_cm"], group["total_added_mass_kg"], "-",
                color=colors.get(n, "0.4"), lw=1.0, alpha=0.8)
        feas = group[group["feasible"]]
        infeas = group[~group["feasible"]]
        ax.plot(feas["d_hp_cm"], feas["total_added_mass_kg"], "o", ms=4.5,
                color=colors.get(n, "0.4"), label=f"{n} pipes")
        ax.plot(infeas["d_hp_cm"], infeas["total_added_mass_kg"], "o", ms=4.5,
                mfc="white", mec=colors.get(n, "0.4"), mew=1.0)
    ax.axhline(90.85, color="0.35", lw=0.9, ls="--")
    ax.text(1.3, 93.5, "Pumped NaK-78 loop (90.9 kg)", fontsize=6.8, color="0.35")
    ax.plot([design.d_hp * 100], [design.total_added_mass], "*", ms=13,
            color="#e6a23c", mec="black", mew=0.4, zorder=5,
            label="Reference")
    ax.set_xlabel("Heat-pipe outer diameter (cm)")
    ax.set_ylabel("Added mass (kg)")
    ax.set_title("ISR-V2 array design space (arterial wick)", fontsize=8.5)
    ax.set_xlim(1.0, 5.2)
    ax.set_ylim(0, 480)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="upper left", frameon=False, fontsize=6.5, ncol=2)
    path = FIG_DIR / "fig_5_6_design_scan.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def make_sensitivity_figure(design: ArrayDesign) -> Path:
    """Sensitivity of peak fuel temperature to heat-pipe and block temperature."""
    plt = _style()
    T_hp = np.linspace(1000.0, 1120.0, 40)
    k_block = np.linspace(40.0, 120.0, 40)
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.0), constrained_layout=True)

    ax = axes[0]
    for kb, color in [(50.0, "#c00000"), (80.0, "#1f4e79"), (110.0, "#2e7d32")]:
        Tmax = [evaluate_reference(design, T_hp=t, k_block=kb).fuel["T_max"]
                for t in T_hp]
        ax.plot(T_hp, Tmax, lw=1.2, color=color, label=f"k_block = {kb:.0f} W/mK")
    ax.axhline(UN_MELTING_LIMIT, color="0.3", lw=0.8, ls="--")
    ax.set_xlabel("Heat-pipe temperature (K)")
    ax.set_ylabel("Peak fuel temperature (K)")
    ax.set_title("(a) Effect of heat-pipe temperature", fontsize=8.5)
    ax.set_ylim(1200, 2600)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="upper left", frameon=False, fontsize=6.5)

    ax = axes[1]
    for thp, color in [(1000.0, "#c00000"), (1050.0, "#1f4e79"),
                       (1100.0, "#2e7d32")]:
        Tmax = [evaluate_reference(design, T_hp=thp, k_block=kb).fuel["T_max"]
                for kb in k_block]
        ax.plot(k_block, Tmax, lw=1.2, color=color, label=f"T_hp = {thp:.0f} K")
    ax.axhline(UN_MELTING_LIMIT, color="0.3", lw=0.8, ls="--")
    ax.set_xlabel("Core-block conductivity (W/mK)")
    ax.set_ylabel("Peak fuel temperature (K)")
    ax.set_title("(b) Effect of block conductivity", fontsize=8.5)
    ax.set_ylim(1200, 2600)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower left", frameon=False, fontsize=6.5)

    path = FIG_DIR / "fig_5_7_sensitivity.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    scan = scan_designs(wick_type="arterial")
    scan.to_csv(RESULTS_DIR / "isr_v2_design_scan.csv", index=False)
    ref = select_reference(scan)
    ref_case = evaluate_reference(ref)
    ref_case.failure_cases.to_csv(RESULTS_DIR / "isr_v2_failure_cases.csv",
                                  index=False)

    in_channel = in_channel_case()
    screen = design_array(19, 0.030, wick_type="screen")
    mass = mass_comparison(ref)
    mass.to_csv(RESULTS_DIR / "isr_v2_mass_comparison.csv", index=False)
    pumped_loop_mass().to_csv(RESULTS_DIR / "pumped_loop_mass_breakdown.csv",
                              index=False)
    heat_pipe_system_mass(ref).to_csv(
        RESULTS_DIR / "heat_pipe_mass_breakdown.csv", index=False)

    summary = {
        "reference": ref_case.as_dict(),
        "in_channel": in_channel,
        "screen_3cm": screen.as_dict(),
        "mass": mass.to_dict(orient="records"),
        "pumped_loop_mass_kg": float(pumped_loop_mass()["mass_kg"].sum()),
        "heat_pipe_system_mass_kg": float(
            heat_pipe_system_mass(ref)["mass_kg"].sum()),
        "reactor_mass_kg": REACTOR_MASS,
    }
    (RESULTS_DIR / "isr_v2_summary.json").write_text(json.dumps(summary, indent=2))

    figs = {
        "limits": str(make_limit_figure(ref)),
        "fuel": str(make_fuel_figure(ref)),
        "startup": str(make_startup_figure(ref)),
        "mass": str(make_mass_figure(ref)),
        "design_scan": str(make_design_scan_figure(scan, ref)),
        "sensitivity": str(make_sensitivity_figure(ref)),
    }
    summary["figures"] = figs
    return summary


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))