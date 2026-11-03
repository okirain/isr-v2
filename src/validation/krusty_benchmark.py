"""Benchmark of the reduced-order thermal model against the KRUSTY 28-hour test.

The benchmark uses the measured fission power and the inferred PCS load state as
boundary conditions.  The model therefore predicts temperatures, not power: the
KRUSTY reactor is a load-following system whose fission power is set by
reactivity feedback, and that neutronic loop is outside the scope of a thermal
ROM.  This separation is stated explicitly in the thesis.

Run from the repository root:

    python -m src.validation.krusty_benchmark
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from ..solver.limits import HeatPipeGeometry
from ..solver.transient_model import (
    HeatPipeTransientModel,
    TransientConfig,
    step_profile,
)

DATA_DIR = Path(__file__).resolve().parent / "data"
RESULTS_DIR = Path(__file__).resolve().parents[2] / "results"
FIG_DIR = Path(__file__).resolve().parents[2] / "thesis_md" / "figures"

# KRUSTY heat-pipe geometry (Poston et al., 2020; Gibson et al., 2020)
KRUSTY_GEOMETRY = HeatPipeGeometry(
    name="KRUSTY Haynes-230/Na",
    L_total=1.00,
    L_evap=0.25,
    L_adiab=0.35,
    L_cond=0.40,
    r_outer=0.00635,      # 1.27 cm OD
    t_wall=0.00089,       # 0.089 cm wall
    t_wick=0.00050,       # thin Ni screen wick
    r_pore=15.0e-6,
    porosity=0.70,
    wire_diameter=50.0e-6,
    tilt_deg=90.0,
    k_wall=20.0,
)

# Reconstructed event sequence of the 28-hour nuclear system test.  Values are
# fission power in kW at the stated time in hours (Poston et al., 2020, Sec.
# IV-X).  The final segment is the scram at T = 28.03 h.
KRUSTY_POWER_SEGMENTS: list[tuple[float, float]] = [
    (0.00, 0.00), (0.15, 0.00), (0.17, 0.02), (0.25, 0.05), (0.33, 0.60),
    (0.38, 2.00), (0.41, 3.70), (0.45, 3.00), (0.63, 3.00), (1.00, 2.90),
    (1.55, 2.90), (1.70, 2.90), (8.00, 2.75), (8.10, 2.05), (9.08, 2.05),
    (9.10, 2.75), (10.02, 2.75), (10.05, 4.05), (11.00, 4.05),
    (11.05, 2.75), (12.00, 2.75), (12.10, 2.55), (12.50, 2.85),
    (13.00, 2.85), (13.10, 2.70), (13.50, 3.00), (14.00, 3.00),
    (14.05, 3.75), (15.00, 3.75), (15.05, 2.70), (16.00, 2.70),
    (16.05, 1.30), (16.20, 1.30), (17.00, 1.30), (17.05, 1.25),
    (18.00, 1.25), (18.05, 4.00), (18.10, 5.10), (18.16, 5.10),
    (18.20, 4.00), (19.00, 4.00), (19.05, 2.60), (20.00, 2.60),
    (20.05, 1.30), (22.70, 1.30), (22.75, 1.50), (24.00, 1.80),
    (25.00, 2.65), (27.00, 2.75), (27.05, 1.35), (28.03, 0.00),
]

# Fraction of nominal PCS heat-rejection conductance.  The converters and
# nitrogen-cooled simulators are the active load; when they are stalled the
# conductance collapses and the core must reject heat through the passive paths.
KRUSTY_SINK_SEGMENTS: list[tuple[float, float]] = [
    (0.00, 1.00), (20.00, 1.00), (20.05, 0.15), (22.70, 0.15),
    (22.75, 0.25), (24.00, 0.25), (24.05, 0.60), (25.00, 0.60),
    (25.05, 1.00), (27.00, 1.00), (27.05, 0.10), (28.03, 0.05),
    (40.00, 0.05),
]


def build_model(**overrides) -> HeatPipeTransientModel:
    """Create the calibrated KRUSTY thermal model."""
    cfg = TransientConfig(
        geometry=KRUSTY_GEOMETRY,
        n_pipes=8,
        n_axial=30,
        fuel_mass_total=42.9,
        fuel_cp=140.0,
        sink_mass_total=12.0,
        sink_cp=500.0,
        h_fuel_wall=60.0,
        h_vapor_evap=40.0,
        h_vapor_cond=40.0,
        h_sink=2.0,
        eps_area_mli=0.00536,
        ua_reject=3.9,
        gravity=0.0,
        enable_limits=True,
        startup_temperature=750.0,
        residual_vapor_conductance=0.2,
        pcs_power_gain=0.8,
        pcs_time_constant=300.0,
        passive_loss_boost=1.5,
    )
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return HeatPipeTransientModel(cfg)


def run_28h(t_end_h: float = 87.0, dt: float = 60.0) -> object:
    """Integrate the reconstructed 28-hour test and the 59-hour cooldown."""
    power = step_profile(
        [(t * 3600.0, p * 1000.0) for t, p in KRUSTY_POWER_SEGMENTS], ramp=30.0
    )
    sink = step_profile(
        [(t * 3600.0, s) for t, s in KRUSTY_SINK_SEGMENTS], ramp=30.0
    )
    model = build_model()
    return model.simulate(t_end_h * 3600.0, power, dt=dt, sink_func=sink)


def _load_state_points() -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / "krusty_statepoints.csv")


def compare(result, state: pd.DataFrame | None = None) -> pd.DataFrame:
    """Interpolate the model at the reported state points and compute errors."""
    if state is None:
        state = _load_state_points()
    t_h = result.t / 3600.0
    state = state.copy()
    state["model_temp_k"] = np.interp(state["time_h"], t_h, result.T_fuel)
    state["error_k"] = state["model_temp_k"] - state["fuel_temp_k"]
    state["relative_error_pct"] = 100.0 * state["error_k"] / state["fuel_temp_k"]
    return state


def _metrics(frame: pd.DataFrame, label: str) -> dict:
    err = frame["error_k"].to_numpy()
    rel = frame["relative_error_pct"].to_numpy()
    return {
        "set": label,
        "n_points": int(len(frame)),
        "mae_k": float(np.mean(np.abs(err))),
        "rmse_k": float(np.sqrt(np.mean(err ** 2))),
        "mape_pct": float(np.mean(np.abs(rel))),
        "max_ape_pct": float(np.max(np.abs(rel))),
    }


def _startup_time(result, threshold_k: float = 773.15) -> float:
    idx = np.where(result.T_fuel >= threshold_k)[0]
    if len(idx) == 0:
        return float("nan")
    return float(result.t[idx[0]] / 3600.0)


def make_figure(result, state: pd.DataFrame) -> Path:
    """Create the validation figure used in Chapter 5."""
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
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9), constrained_layout=True)
    t_h = result.t / 3600.0

    ax = axes[0]
    mask = t_h <= 20.0
    ax.plot(t_h[mask], result.T_fuel_c[mask], color="#1f4e79", lw=1.2,
            label="Model, fuel")
    ax.plot(t_h[mask], result.T_vapor_c[mask], color="#c00000", lw=1.0,
            ls="--", label="Model, vapour")
    pts = state[state["case_id"].isin(["startup_500C", "quasi_steady_1p7h", "nominal_8h"])]
    ax.plot(pts["time_h"], pts["fuel_temp_c"], "o", ms=4.5, mfc="white",
            mec="black", mew=0.8, label="KRUSTY")
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 950)
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("(a) Startup and load changes", fontsize=8.5)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(loc="lower right", frameon=False, fontsize=6.5)

    ax = axes[1]
    mask = (t_h >= 19.0) & (t_h <= 30.0)
    ax.plot(t_h[mask], result.T_fuel_c[mask], color="#1f4e79", lw=1.2)
    pts = state[state["case_id"].isin(["loss_cooling_20p5h", "stalled_24h", "pre_scram_27h"])]
    ax.plot(pts["time_h"], pts["fuel_temp_c"], "o", ms=4.5, mfc="white",
            mec="black", mew=0.8)
    ax.set_xlim(19, 30)
    ax.set_ylim(600, 900)
    ax.set_xlabel("Time (h)")
    ax.set_title("(b) Loss of active cooling", fontsize=8.5)
    ax.grid(alpha=0.25, lw=0.5)

    ax = axes[2]
    mask = t_h >= 28.0
    ax.plot(t_h[mask], result.T_fuel_c[mask], color="#1f4e79", lw=1.2)
    pts = state[state["case_id"].isin(["cooldown_31h", "cooldown_87h"])]
    ax.plot(pts["time_h"], pts["fuel_temp_c"], "o", ms=4.5, mfc="white",
            mec="black", mew=0.8)
    ax.set_xlim(28, 87)
    ax.set_ylim(0, 900)
    ax.set_xlabel("Time (h)")
    ax.set_title("(c) Post-scram cooldown", fontsize=8.5)
    ax.grid(alpha=0.25, lw=0.5)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    path = FIG_DIR / "fig_5_1_krusty_validation.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    result = run_28h()
    state = compare(result)

    out = pd.DataFrame({
        "time_h": result.t / 3600.0,
        "power_kw": [np.nan] * len(result.t),
        "T_fuel_k": result.T_fuel,
        "T_fuel_c": result.T_fuel_c,
        "T_vapor_k": result.T_vapor,
        "T_sink_k": result.T_sink,
        "Q_hp_w": result.Q_hp,
        "Q_limit_w": result.Q_limit,
        "limiting_mode": result.limiting_mode,
    })
    power = step_profile(
        [(t * 3600.0, p * 1000.0) for t, p in KRUSTY_POWER_SEGMENTS], ramp=30.0
    )
    out["power_kw"] = [power(t) / 1000.0 for t in result.t]
    out.to_csv(RESULTS_DIR / "krusty_validation.csv", index=False)
    state.to_csv(RESULTS_DIR / "krusty_statepoint_comparison.csv", index=False)

    envelope = state[state["include_in_envelope"] == 1]
    metrics = {
        "geometry": asdict(KRUSTY_GEOMETRY),
        "n_pipes": 8,
        "envelope": _metrics(envelope, "quasi-steady envelope"),
        "all_points": _metrics(state, "all reported points"),
        "startup_500c": {
            "measured_h": 0.63,
            "model_h": _startup_time(result),
            "relative_error_pct": 100.0 * (_startup_time(result) - 0.63) / 0.63,
        },
        "nominal_8h": {
            "model_k": float(state.loc[state.case_id == "nominal_8h", "model_temp_k"].iloc[0]),
            "measured_k": 1073.15,
        },
    }
    (RESULTS_DIR / "krusty_validation_metrics.json").write_text(
        json.dumps(metrics, indent=2)
    )
    fig = make_figure(result, state)
    metrics["figure"] = str(fig)
    return metrics


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
