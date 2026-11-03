"""Operational limits of a high-temperature liquid-metal heat pipe.

The module implements the classical capillary, sonic, entrainment, viscous and
boiling limits used in heat-pipe design (Chi, 1976; Dunn and Reay, 2012;
Zohuri, 2016).  All limits are returned in watts for a single heat pipe.

For a microgravity application the gravitational pressure head is set to zero
and the capillary limit becomes the only hydrostatic constraint.  The same
model can be run with g = 9.81 m s^-2 and a thermosiphon-assisted return for
ground validation against KRUSTY.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, radians, sin

import numpy as np

from ..thermophysical.properties import (
    Fluid,
    SODIUM,
    dynamic_viscosity,
    get_fluid,
    latent_heat,
    liquid_density,
    saturation_pressure,
    saturation_temperature,
    surface_tension,
    thermal_conductivity_liquid,
    vapor_density,
    vapor_viscosity,
)


@dataclass(frozen=True)
class HeatPipeGeometry:
    """Axial and radial geometry of a cylindrical heat pipe.

    Parameters are in SI units.  The wick is a homogeneous annular screen with
    an effective pore radius r_pore and porosity porosity.
    """

    name: str
    L_total: float          # m
    L_evap: float           # m
    L_adiab: float          # m
    L_cond: float           # m
    r_outer: float          # m
    t_wall: float           # m
    t_wick: float           # m
    r_pore: float           # m
    porosity: float = 0.70
    wire_diameter: float = 5.0e-5
    contact_angle_deg: float = 0.0
    tilt_deg: float = 90.0
    k_wall: float = 20.0    # W m^-1 K^-1, Haynes 230 at 1000 K
    permeability_override: float | None = None  # m^2, for arterial/groove wicks

    def __post_init__(self):
        if abs(self.L_evap + self.L_adiab + self.L_cond - self.L_total) > 1e-9:
            raise ValueError("L_evap + L_adiab + L_cond must equal L_total")
        if self.r_outer <= self.t_wall + self.t_wick:
            raise ValueError("Vapour core radius must be positive")
        if not (0.0 < self.porosity < 1.0):
            raise ValueError("Porosity must lie in (0, 1)")

    @property
    def r_wall_inner(self) -> float:
        return self.r_outer - self.t_wall

    @property
    def r_vapor(self) -> float:
        return self.r_wall_inner - self.t_wick

    @property
    def area_vapor(self) -> float:
        return pi * self.r_vapor ** 2

    @property
    def area_wick(self) -> float:
        return pi * (self.r_wall_inner ** 2 - self.r_vapor ** 2)

    @property
    def area_wall(self) -> float:
        return pi * (self.r_outer ** 2 - self.r_wall_inner ** 2)

    @property
    def area_evap_outer(self) -> float:
        return 2.0 * pi * self.r_outer * self.L_evap

    @property
    def area_evap_inner(self) -> float:
        return 2.0 * pi * self.r_vapor * self.L_evap

    @property
    def area_cond_outer(self) -> float:
        return 2.0 * pi * self.r_outer * self.L_cond

    @property
    def effective_length(self) -> float:
        """Chi's effective length: half evaporator + adiabatic + half condenser."""
        return 0.5 * self.L_evap + self.L_adiab + 0.5 * self.L_cond

    @property
    def permeability(self) -> float:
        """Wick permeability in m^2.

        A Blake-Kozeny screen correlation is used by default; an arterial or
        grooved wick can supply a measured/design permeability directly.
        """
        if self.permeability_override is not None:
            return float(self.permeability_override)
        eps = self.porosity
        return self.wire_diameter ** 2 * eps ** 3 / (150.0 * (1.0 - eps) ** 2)

    @property
    def hydraulic_radius(self) -> float:
        """Hydraulic radius of the vapour core in m."""
        return 0.5 * self.r_vapor


@dataclass(frozen=True)
class LimitResult:
    """Container for one operating-limit evaluation."""

    temperature: float
    q_max: float
    limiting_mode: str
    capillary: float
    sonic: float
    entrainment: float
    viscous: float
    boiling: float
    dp_capillary: float
    dp_liquid: float
    dp_vapor: float
    dp_gravity: float
    q_operating: float = 0.0

    @property
    def margin(self) -> float:
        if self.q_operating <= 0:
            return float("inf")
        return self.q_max / self.q_operating


def _q_limit(value: float) -> float:
    return float(max(value, 0.0))


def capillary_limit(
    T: float,
    geom: HeatPipeGeometry,
    fluid: Fluid | str = SODIUM,
    g: float = 0.0,
    return_components: bool = False,
):
    """Capillary-limited heat transport in W.

    The available capillary pressure 2 sigma cos(theta) / r_pore must exceed the
    sum of the liquid wick, vapour core and gravitational pressure drops.
    """
    f = get_fluid(fluid)
    sigma = float(surface_tension(T, f))
    rho_l = float(liquid_density(T, f))
    mu_l = float(dynamic_viscosity(T, f))
    h_fg = float(latent_heat(T, f))
    rho_v = float(vapor_density(T, f))
    mu_v = float(vapor_viscosity(T, f))
    L_eff = geom.effective_length

    dp_cap = 2.0 * sigma * cos(radians(geom.contact_angle_deg)) / geom.r_pore
    dp_g = rho_l * g * L_eff * sin(radians(geom.tilt_deg))
    coef_l = mu_l * L_eff / (geom.permeability * geom.area_wick * rho_l * h_fg)
    coef_v = 8.0 * mu_v * L_eff / (pi * geom.r_vapor ** 4 * rho_v * h_fg)
    q = _q_limit((dp_cap - dp_g) / (coef_l + coef_v))
    dp_l = coef_l * q
    dp_v = coef_v * q
    if return_components:
        return q, {
            "dp_capillary": dp_cap,
            "dp_liquid": dp_l,
            "dp_vapor": dp_v,
            "dp_gravity": dp_g,
        }
    return q


def sonic_limit(T: float, geom: HeatPipeGeometry, fluid: Fluid | str = SODIUM) -> float:
    """Sonic (choked-vapour) limit in W.

    Q = A_v rho_v h_fg [ gamma R_v T / (2 (gamma + 1)) ]^(1/2)
    """
    f = get_fluid(fluid)
    rho_v = float(vapor_density(T, f))
    h_fg = float(latent_heat(T, f))
    Rv = f.R_specific
    gamma = f.gamma
    c = np.sqrt(gamma * Rv * T / (2.0 * (gamma + 1.0)))
    return _q_limit(geom.area_vapor * rho_v * h_fg * c)


def entrainment_limit(T: float, geom: HeatPipeGeometry, fluid: Fluid | str = SODIUM) -> float:
    """Entrainment (Weber) limit in W.

    Liquid droplets are entrained when the vapour Weber number exceeds unity.
    The characteristic length is the wick pore radius, following the
    conservative form used by Chi (1976).
    """
    f = get_fluid(fluid)
    sigma = float(surface_tension(T, f))
    rho_v = float(vapor_density(T, f))
    h_fg = float(latent_heat(T, f))
    q = geom.area_vapor * h_fg * np.sqrt(sigma * rho_v / (2.0 * geom.r_pore))
    return _q_limit(q)


def viscous_limit(T: float, geom: HeatPipeGeometry, fluid: Fluid | str = SODIUM) -> float:
    """Viscous (vapour-pressure) limit in W."""
    f = get_fluid(fluid)
    rho_v = float(vapor_density(T, f))
    h_fg = float(latent_heat(T, f))
    mu_v = float(vapor_viscosity(T, f))
    p_v = float(saturation_pressure(T, f))
    q = pi * geom.r_vapor ** 4 * rho_v * h_fg * p_v / (8.0 * mu_v * geom.effective_length)
    return _q_limit(q)


def boiling_limit(T: float, geom: HeatPipeGeometry, fluid: Fluid | str = SODIUM,
                  r_nucleation: float = 1.0e-6) -> float:
    """Wick boiling limit in W.

    The criterion is that the superheat needed to activate a vapour nucleus of
    radius r_nucleation must not be reached across the wick thickness.
    """
    f = get_fluid(fluid)
    sigma = float(surface_tension(T, f))
    rho_v = float(vapor_density(T, f))
    h_fg = float(latent_heat(T, f))
    k_l = float(thermal_conductivity_liquid(T, f))
    k_eff = geom.porosity * k_l + (1.0 - geom.porosity) * geom.k_wall
    dT_crit = 2.0 * sigma * T / (r_nucleation * rho_v * h_fg)
    q_flux = k_eff * dT_crit / geom.t_wick
    return _q_limit(q_flux * geom.area_evap_inner)


def operating_limits(
    T: float,
    geom: HeatPipeGeometry,
    fluid: Fluid | str = SODIUM,
    g: float = 0.0,
    q_operating: float = 0.0,
) -> LimitResult:
    """Evaluate all limits and return the binding one."""
    q_cap, dp = capillary_limit(T, geom, fluid, g=g, return_components=True)
    q_son = sonic_limit(T, geom, fluid)
    q_ent = entrainment_limit(T, geom, fluid)
    q_vis = viscous_limit(T, geom, fluid)
    q_boil = boiling_limit(T, geom, fluid)
    candidates = {
        "capillary": q_cap,
        "sonic": q_son,
        "entrainment": q_ent,
        "viscous": q_vis,
        "boiling": q_boil,
    }
    mode = min(candidates, key=candidates.get)
    return LimitResult(
        temperature=float(T),
        q_max=float(candidates[mode]),
        limiting_mode=mode,
        capillary=q_cap,
        sonic=q_son,
        entrainment=q_ent,
        viscous=q_vis,
        boiling=q_boil,
        dp_capillary=dp["dp_capillary"],
        dp_liquid=dp["dp_liquid"],
        dp_vapor=dp["dp_vapor"],
        dp_gravity=dp["dp_gravity"],
        q_operating=float(q_operating),
    )


def limit_curve(
    temperatures,
    geom: HeatPipeGeometry,
    fluid: Fluid | str = SODIUM,
    g: float = 0.0,
) -> dict:
    """Evaluate all limits over an array of temperatures."""
    T = np.asarray(temperatures, dtype=float)
    rows = [operating_limits(t, geom, fluid, g=g) for t in T]
    return {
        "T": T,
        "capillary": np.array([r.capillary for r in rows]),
        "sonic": np.array([r.sonic for r in rows]),
        "entrainment": np.array([r.entrainment for r in rows]),
        "viscous": np.array([r.viscous for r in rows]),
        "boiling": np.array([r.boiling for r in rows]),
        "q_max": np.array([r.q_max for r in rows]),
        "limiting_mode": [r.limiting_mode for r in rows],
    }


def required_temperature(
    q_target: float,
    geom: HeatPipeGeometry,
    fluid: Fluid | str = SODIUM,
    g: float = 0.0,
    T_lo: float | None = None,
    T_hi: float = 1400.0,
) -> float:
    """Lowest temperature at which the pipe can transport q_target.

    A monotone bracket is constructed by scanning; the first crossing is refined
    with Brent's method.  Returns nan when the target is unreachable.
    """
    from scipy.optimize import brentq

    f = get_fluid(fluid)
    lo = (f.T_melt + 5.0) if T_lo is None else T_lo
    grid = np.linspace(lo, T_hi, 240)
    q = np.array([operating_limits(t, geom, f, g=g).q_max for t in grid])
    above = np.where(q >= q_target)[0]
    if len(above) == 0:
        return float("nan")
    i = int(above[0])
    if i == 0:
        return float(grid[0])
    a, b = grid[i - 1], grid[i]
    return float(brentq(
        lambda t: operating_limits(t, geom, f, g=g).q_max - q_target,
        a, b, xtol=1e-6, rtol=1e-10, maxiter=100,
    ))
