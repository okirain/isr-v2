"""ISR core heat extraction and passive heat-pipe array design.

The baseline Indonesian Space Reactor (ISR) of Khandaq et al. (2020) is a
500 kWth, NaK-78 cooled fast reactor with 61 annular UN fuel pins.  This module
answers two questions:

1. Can the existing 6 mm central coolant channel be replaced one-for-one by a
   sodium heat pipe?
2. If not, what heat-pipe array can remove 500 kWth while keeping the UN
   centreline below 2523 K?

The fuel-temperature model is a one-dimensional radial finite-volume solution
of the steady heat-conduction equation with temperature-dependent UN
conductivity.  The heat-pipe limits come from the solver.limits module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import pi, sqrt

import numpy as np
import pandas as pd

from ..solver.limits import HeatPipeGeometry, LimitResult, operating_limits
from ..thermophysical.properties import SODIUM

# ---------------------------------------------------------------------------
# Reference data from Khandaq et al. (2020)
# ---------------------------------------------------------------------------
N_PINS = 61
P_TOTAL = 500.0e3              # W
FUEL_R_INNER = 0.003825        # m, UN annulus inner radius
FUEL_R_OUTER = 0.01630         # m, UN annulus outer radius
FUEL_LENGTH = 0.250            # m, active UN length
FUEL_MASS_PIN = 2.55           # kg
UN_DENSITY = 14.264e3          # kg m^-3, theoretical
UN_MELTING_LIMIT = 2523.0      # K, design limit used in this work
CORE_PITCH = 0.0340            # m, hexagonal fuel-pin pitch
CORE_BLOCK_MASS = 50.80        # kg
REFLECTOR_MASS = 194.37        # kg
REACTOR_MASS = 447.12          # kg, Khandaq et al. Table 2
T_INLET = 750.0                # K
T_OUTLET = 950.0               # K
T_HOT_PIPE = 1050.0            # K, reference sodium heat-pipe temperature


def un_thermal_conductivity(T):
    """Thermal conductivity of uranium mononitride in W m^-1 K^-1.

    A mild linear fit to the recommended UN data (15-23 W m^-1 K^-1 between
    800 and 2500 K).  The exact value is not critical because the Biot number of
    the pin is small; the thesis reports a sensitivity to this coefficient.
    """
    T = np.asarray(T, dtype=float)
    k = 21.0 - 0.0025 * (T - 800.0)
    return np.clip(k, 12.0, 25.0)


@dataclass(frozen=True)
class FuelPin:
    """Annular UN fuel pin of the ISR."""

    r_inner: float = FUEL_R_INNER
    r_outer: float = FUEL_R_OUTER
    length: float = FUEL_LENGTH
    mass: float = FUEL_MASS_PIN
    n_pins: int = N_PINS

    @property
    def volume(self) -> float:
        return pi * (self.r_outer ** 2 - self.r_inner ** 2) * self.length

    @property
    def density(self) -> float:
        return self.mass / self.volume

    @property
    def heat_generation(self) -> float:
        return P_TOTAL / (self.n_pins * self.volume)


def _face_conductance(r_left, r_right, k, length):
    """Cylindrical conduction conductance between two radii, W K^-1."""
    return 2.0 * pi * k * length / np.log(r_right / r_left)


def fuel_temperature_profile(
    q_pin: float,
    T_boundary: float,
    pin: FuelPin = FuelPin(),
    boundary: str = "outer",
    n: int = 160,
    n_iter: int = 25,
) -> dict:
    """Steady radial temperature profile of one annular fuel pin.

    Parameters
    ----------
    q_pin : float
        Total heat generated in the pin, W.
    T_boundary : float
        Temperature of the cooled boundary, K.
    boundary : {"outer", "inner"}
        "outer" for a pin cooled at its outer surface (heat-pipe array in the
        core block); "inner" for a pin cooled at its central channel (the
        original NaK-78 channel or a one-for-one in-channel heat pipe).

    Returns
    -------
    dict with r, T, T_max, T_min, T_mean and the radial temperature drop.
    """
    if boundary not in ("outer", "inner"):
        raise ValueError("boundary must be 'outer' or 'inner'")
    ri, ro = pin.r_inner, pin.r_outer
    # Cell faces, geometrically spaced for a better near-surface resolution
    faces = np.linspace(ri, ro, n + 1)
    r = 0.5 * (faces[:-1] + faces[1:])
    dr = np.diff(faces)
    volumes = pi * (faces[1:] ** 2 - faces[:-1] ** 2) * pin.length
    q_vol = q_pin / pin.volume
    q_nodes = q_vol * volumes

    T = np.full(n, float(T_boundary))
    for _ in range(n_iter):
        k = un_thermal_conductivity(T)
        # Face conductances
        g = np.array([
            _face_conductance(faces[i], faces[i + 1], 0.5 * (k[i] + k[i + 1]),
                              pin.length)
            for i in range(n - 1)
        ])
        A = np.zeros((n, n))
        b = q_nodes.copy()
        for i in range(n):
            if i > 0:
                A[i, i - 1] -= g[i - 1]
                A[i, i] += g[i - 1]
            if i < n - 1:
                A[i, i + 1] -= g[i]
                A[i, i] += g[i]
        if boundary == "inner":
            gb = _face_conductance(ri, r[0], k[0], pin.length)
            A[0, 0] += gb
            b[0] += gb * T_boundary
        else:
            gb = _face_conductance(r[-1], ro, k[-1], pin.length)
            A[-1, -1] += gb
            b[-1] += gb * T_boundary
        T_new = np.linalg.solve(A, b)
        if np.max(np.abs(T_new - T)) < 1e-6:
            T = T_new
            break
        T = T_new
    return {
        "r": r,
        "T": T,
        "T_max": float(T.max()),
        "T_min": float(T.min()),
        "T_mean": float(np.average(T, weights=volumes)),
        "delta_T": float(T.max() - T.min()),
    }


# ---------------------------------------------------------------------------
# Heat-pipe array design
# ---------------------------------------------------------------------------
@dataclass
class ArrayDesign:
    """A candidate ISR-V2 heat-pipe array."""

    n_hp: int
    d_hp: float
    q_per_hp: float
    q_max: float
    margin: float
    limiting_mode: str
    hp_mass_total: float
    core_diameter: float
    reflector_mass_added: float
    block_mass_added: float
    header_mass: float
    total_added_mass: float
    feasible: bool
    geometry: HeatPipeGeometry
    limits: LimitResult
    notes: str = ""

    def as_dict(self) -> dict:
        return {
            "n_hp": self.n_hp,
            "d_hp_cm": self.d_hp * 100.0,
            "q_per_hp_kw": self.q_per_hp / 1000.0,
            "q_max_kw": self.q_max / 1000.0,
            "margin": self.margin,
            "limiting_mode": self.limiting_mode,
            "hp_mass_total_kg": self.hp_mass_total,
            "core_diameter_cm": self.core_diameter * 100.0,
            "reflector_mass_added_kg": self.reflector_mass_added,
            "block_mass_added_kg": self.block_mass_added,
            "header_mass_kg": self.header_mass,
            "total_added_mass_kg": self.total_added_mass,
            "feasible": self.feasible,
            "notes": self.notes,
        }


WICK_PERMEABILITY = {
    "screen": None,          # Blake-Kozeny screen, ~6.4e-11 m^2
    "arterial": 1.0e-9,      # arterial/advanced wick, design assumption
    "groove": 5.0e-9,        # axial groove, upper-bound estimate
}


def make_geometry(
    d_hp: float,
    L_evap: float = 0.30,
    L_adiab: float = 0.10,
    L_cond: float = 0.40,
    t_wall: float = 0.0015,
    t_wick: float = 0.0015,
    r_pore: float = 10.0e-6,
    porosity: float = 0.70,
    wick_type: str = "screen",
) -> HeatPipeGeometry:
    """Build a cylindrical sodium heat-pipe geometry for the ISR core."""
    r_outer = 0.5 * d_hp
    if r_outer <= t_wall + t_wick:
        raise ValueError("Heat-pipe diameter is too small for the chosen wall/wick")
    if wick_type not in WICK_PERMEABILITY:
        raise ValueError(f"Unknown wick type {wick_type!r}")
    return HeatPipeGeometry(
        name=f"ISR-V2 Na HP D={d_hp*100:.1f} cm ({wick_type})",
        L_total=L_evap + L_adiab + L_cond,
        L_evap=L_evap,
        L_adiab=L_adiab,
        L_cond=L_cond,
        r_outer=r_outer,
        t_wall=t_wall,
        t_wick=t_wick,
        r_pore=r_pore,
        porosity=porosity,
        wire_diameter=50.0e-6,
        tilt_deg=0.0,
        k_wall=22.0,
        permeability_override=WICK_PERMEABILITY[wick_type],
    )


def heat_pipe_mass(geom: HeatPipeGeometry, rho_wall: float = 8900.0,
                   rho_wick: float = 8900.0, rho_fluid: float = 770.0) -> float:
    """Mass of one heat pipe (shell + dry wick + sodium inventory) in kg."""
    m_shell = rho_wall * geom.area_wall * geom.L_total
    m_wick = rho_wick * (1.0 - geom.porosity) * geom.area_wick * geom.L_total
    m_fluid = rho_fluid * geom.porosity * geom.area_wick * geom.L_total
    return m_shell + m_wick + m_fluid


def core_diameter(n_hp: int, d_hp: float, pitch: float = CORE_PITCH,
                  margin: float = 1.10) -> float:
    """Equivalent circular core diameter with interleaved fuel pins and pipes."""
    a_pin = sqrt(3.0) / 2.0 * pitch ** 2
    a_hp = sqrt(3.0) / 2.0 * (d_hp + 0.006) ** 2
    area = margin * (N_PINS * a_pin + n_hp * a_hp)
    return sqrt(4.0 * area / pi)


D_CORE_BASELINE = core_diameter(0, 0.0, margin=1.10)


def design_array(
    n_hp: int,
    d_hp: float,
    P_total: float = P_TOTAL,
    T_op: float = T_HOT_PIPE,
    safety_factor: float = 1.5,
    header_mass: float = 8.0,
    wick_type: str = "arterial",
) -> ArrayDesign:
    """Evaluate one heat-pipe array candidate."""
    geom = make_geometry(d_hp, wick_type=wick_type)
    q = P_total / n_hp
    lim = operating_limits(T_op, geom, SODIUM, g=0.0, q_operating=q)
    hp_mass = heat_pipe_mass(geom) * n_hp
    d_core = core_diameter(n_hp, d_hp)
    reflector_add = REFLECTOR_MASS * (d_core / D_CORE_BASELINE - 1.0)
    block_add = CORE_BLOCK_MASS * ((d_core / D_CORE_BASELINE) ** 2 - 1.0)
    total = hp_mass + header_mass + reflector_add + block_add
    feasible = (lim.q_max >= safety_factor * q) and (d_core <= 0.55)
    notes = []
    if lim.q_max < safety_factor * q:
        notes.append("transport margin below 1.5")
    if d_core > 0.55:
        notes.append("core diameter exceeds 55 cm")
    return ArrayDesign(
        n_hp=n_hp, d_hp=d_hp, q_per_hp=q, q_max=lim.q_max,
        margin=lim.q_max / q, limiting_mode=lim.limiting_mode,
        hp_mass_total=hp_mass, core_diameter=d_core,
        reflector_mass_added=reflector_add, block_mass_added=block_add,
        header_mass=header_mass, total_added_mass=total,
        feasible=feasible, geometry=geom, limits=lim,
        notes="; ".join(notes),
    )


def scan_designs(
    n_values=(13, 19, 37, 61),
    d_values=(0.0127, 0.020, 0.030, 0.040, 0.050),
    safety_factor: float = 1.5,
    wick_type: str = "arterial",
) -> pd.DataFrame:
    """Scan the discrete design space and return a comparison table."""
    rows = []
    for n in n_values:
        for d in d_values:
            try:
                design = design_array(n, d, safety_factor=safety_factor,
                                      wick_type=wick_type)
            except ValueError as exc:
                rows.append({"n_hp": n, "d_hp_cm": d * 100, "feasible": False,
                             "notes": str(exc)})
                continue
            row = design.as_dict()
            rows.append(row)
    return pd.DataFrame(rows)


def select_reference(scan: pd.DataFrame) -> ArrayDesign:
    """Pick the lightest feasible array with at least 19 independent pipes."""
    feasible = scan[(scan["feasible"]) & (scan["n_hp"] >= 19)]
    if feasible.empty:
        raise RuntimeError("No feasible array in the scanned design space")
    row = feasible.sort_values("total_added_mass_kg").iloc[0]
    return design_array(int(row["n_hp"]), float(row["d_hp_cm"]) / 100.0)


@dataclass
class ReferenceCase:
    """Thermal state of the selected ISR-V2 reference core."""

    design: ArrayDesign
    fuel: dict
    T_block: float
    T_hp: float
    q_flux_block: float
    block_delta_T: float
    fuel_margin: float
    failure_cases: pd.DataFrame = field(default_factory=pd.DataFrame)

    def as_dict(self) -> dict:
        return {
            "n_hp": self.design.n_hp,
            "d_hp_cm": self.design.d_hp * 100.0,
            "q_per_hp_kw": self.design.q_per_hp / 1000.0,
            "q_max_kw": self.design.q_max / 1000.0,
            "margin": self.design.margin,
            "core_diameter_cm": self.design.core_diameter * 100.0,
            "T_hp_k": self.T_hp,
            "T_block_k": self.T_block,
            "T_fuel_max_k": self.fuel["T_max"],
            "T_fuel_mean_k": self.fuel["T_mean"],
            "fuel_margin_k": self.fuel_margin,
            "fuel_margin_pct": 100.0 * self.fuel_margin / UN_MELTING_LIMIT,
            "block_delta_T_k": self.block_delta_T,
        }


def evaluate_reference(
    design: ArrayDesign,
    T_hp: float = T_HOT_PIPE,
    k_block: float = 80.0,
    spreading_fraction: float = 0.5,
    boundary: str = "outer",
) -> ReferenceCase:
    """Evaluate fuel temperature and failure tolerance of a reference array.

    The core-block temperature drop is estimated with a one-dimensional
    spreading resistance: q'' L / k_block, where L is a fraction of the
    heat-pipe pitch and q'' is the evaporator heat flux.
    """
    q = design.q_per_hp
    d = design.d_hp
    L_evap = design.geometry.L_evap
    q_flux = q / (pi * d * L_evap)
    L_spread = spreading_fraction * (d + 0.006)
    dT_block = q_flux * L_spread / k_block
    T_block = T_hp + dT_block
    # Each fuel pin generates P_TOTAL/N_PINS; one heat pipe serves several pins.
    q_pin = P_TOTAL / N_PINS
    fuel = fuel_temperature_profile(q_pin, T_block, boundary=boundary)
    margin = UN_MELTING_LIMIT - fuel["T_max"]

    # Failure: one pipe lost -> six neighbours share the load.  The local
    # spreading length doubles; three adjacent failures triple it.
    rows = []
    for n_failed, factor in [(0, 1.0), (1, 2.0), (3, 3.0)]:
        q_local = q * (1.0 + n_failed / 6.0)
        dT_block_f = dT_block * factor
        T_block_f = T_hp + dT_block_f
        fuel_f = fuel_temperature_profile(q_pin, T_block_f,
                                           boundary=boundary)
        rows.append({
            "failed_pipes": n_failed,
            "q_per_surviving_pipe_kw": q_local / 1000.0,
            "T_block_k": T_block_f,
            "T_fuel_max_k": fuel_f["T_max"],
            "margin_to_limit_k": UN_MELTING_LIMIT - fuel_f["T_max"],
            "within_limit": fuel_f["T_max"] < UN_MELTING_LIMIT,
        })
    return ReferenceCase(
        design=design, fuel=fuel, T_block=T_block, T_hp=T_hp,
        q_flux_block=q_flux, block_delta_T=dT_block,
        fuel_margin=margin, failure_cases=pd.DataFrame(rows),
    )


def in_channel_case(d_hp: float = 0.006, T_hp: float = T_HOT_PIPE) -> dict:
    """One-for-one replacement in the existing 6 mm central channel."""
    geom = make_geometry(d_hp, L_evap=0.25, L_adiab=0.05, L_cond=0.20,
                         t_wall=0.0005, t_wick=0.0005)
    q = P_TOTAL / N_PINS
    lim = operating_limits(T_hp, geom, SODIUM, g=0.0, q_operating=q)
    fuel = fuel_temperature_profile(q, T_hp, boundary="inner")
    return {
        "configuration": "V2-A in-channel",
        "d_hp_cm": d_hp * 100.0,
        "q_per_pin_kw": q / 1000.0,
        "q_max_kw": lim.q_max / 1000.0,
        "margin": lim.q_max / q,
        "limiting_mode": lim.limiting_mode,
        "T_fuel_max_k": fuel["T_max"],
        "feasible": lim.q_max >= 1.5 * q,
    }


def print_summary() -> None:
    scan = scan_designs()
    print(scan[["n_hp", "d_hp_cm", "q_per_hp_kw", "q_max_kw", "margin",
                "limiting_mode", "total_added_mass_kg", "feasible"]]
          .to_string(index=False))
    print("\nIn-channel case:")
    print(in_channel_case())


if __name__ == "__main__":
    print_summary()
