"""Thermophysical properties of liquid-metal heat-pipe working fluids.

The module provides smooth, differentiable correlations for saturated sodium
and potassium between their melting points and 1400 K. The formulations are
deliberately kept analytical so that the reduced-order model can be
differentiated and inverted cheaply.

Primary sources
---------------
* J. K. Fink and L. Leibowitz, "Thermodynamic and transport properties of
  sodium liquid and vapor," ANL/RE-95/2, Argonne National Laboratory, 1995.
* R. W. Ohse (ed.), Handbook of Thermodynamic and Transport Properties of
  Alkali Metals, Blackwell, 1985.
* IAEA, Thermophysical Properties of Materials for Nuclear Engineering,
  Vienna, 2008.

The correlations are engineering fits valid to roughly +/-5 % over the
800-1300 K range that is relevant to space heat pipes. They are not intended
for precision vapour-pressure metrology.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

R_UNIVERSAL = 8.31446261815324  # J mol^-1 K^-1
P_STANDARD = 101325.0           # Pa


@dataclass(frozen=True)
class Fluid:
    """Container for the coefficients of one alkali-metal working fluid."""

    name: str
    symbol: str
    molar_mass: float          # kg mol^-1
    T_melt: float              # K
    T_boil: float              # K, at 101325 Pa
    h_fg_boil: float           # J kg^-1 at T_boil
    dh_fg_dT: float            # J kg^-1 K^-1
    rho_a: float               # kg m^-3, rho = a + b T
    rho_b: float
    mu_a: float                # Pa s, mu = a exp(b/T)
    mu_b: float
    sigma_a: float             # N m^-1, sigma = a + b T
    sigma_b: float
    cp_a: float                # J kg^-1 K^-1, cp = a + b (T - 800)
    cp_b: float
    k_a: float                 # W m^-1 K^-1, k = a + b T
    k_b: float
    gamma: float = 5.0 / 3.0   # monatomic vapour

    @property
    def R_specific(self) -> float:
        """Specific gas constant of the vapour, J kg^-1 K^-1."""
        return R_UNIVERSAL / self.molar_mass


SODIUM = Fluid(
    name="Sodium", symbol="Na",
    molar_mass=0.02298977, T_melt=370.94, T_boil=1156.1,
    h_fg_boil=3.60e6, dh_fg_dT=-1200.0,
    rho_a=1014.0, rho_b=-0.235,
    mu_a=7.0e-5, mu_b=1100.0,
    sigma_a=0.235, sigma_b=-8.0e-5,
    cp_a=1260.0, cp_b=0.10,
    k_a=85.0, k_b=-0.020,
)

POTASSIUM = Fluid(
    name="Potassium", symbol="K",
    molar_mass=0.0390983, T_melt=336.70, T_boil=1032.0,
    h_fg_boil=1.95e6, dh_fg_dT=-900.0,
    rho_a=900.0, rho_b=-0.220,
    mu_a=5.0e-5, mu_b=1200.0,
    sigma_a=0.135, sigma_b=-4.0e-5,
    cp_a=800.0, cp_b=0.05,
    k_a=55.0, k_b=-0.012,
)

FLUIDS = {"Na": SODIUM, "K": POTASSIUM, "sodium": SODIUM, "potassium": POTASSIUM}


def get_fluid(fluid: Fluid | str) -> Fluid:
    """Return a Fluid from an instance or a short name."""
    if isinstance(fluid, Fluid):
        return fluid
    try:
        return FLUIDS[fluid]
    except KeyError as exc:  # pragma: no cover - defensive
        raise KeyError(f"Unknown working fluid {fluid!r}; choose 'Na' or 'K'") from exc


def _arr(T):
    return np.asarray(T, dtype=float)


def _scalar_or_array(value, T):
    if np.ndim(T) == 0:
        return float(np.asarray(value))
    return value


def saturation_pressure(T, fluid: Fluid | str = SODIUM):
    """Saturation pressure in Pa from a temperature-dependent Clausius-Clapeyron fit.

    The latent heat h_fg(T) = h0 + a (T - Tb) is integrated analytically, so the
    correlation is thermodynamically consistent with latent_heat and reproduces
    the normal boiling point exactly.
    """
    f = get_fluid(fluid)
    T = _arr(T)
    if np.any(T <= 0):
        raise ValueError("Temperature must be positive")
    Rv = f.R_specific
    h0, a, Tb = f.h_fg_boil, f.dh_fg_dT, f.T_boil
    ln_p = (
        np.log(P_STANDARD)
        - (h0 - a * Tb) / Rv * (1.0 / T - 1.0 / Tb)
        + a / Rv * np.log(T / Tb)
    )
    return _scalar_or_array(np.exp(ln_p), T)


def latent_heat(T, fluid: Fluid | str = SODIUM):
    """Latent heat of vaporisation in J kg^-1."""
    f = get_fluid(fluid)
    T = _arr(T)
    h = f.h_fg_boil + f.dh_fg_dT * (T - f.T_boil)
    return _scalar_or_array(np.maximum(h, 0.05 * f.h_fg_boil), T)


def liquid_density(T, fluid: Fluid | str = SODIUM):
    """Saturated liquid density in kg m^-3."""
    f = get_fluid(fluid)
    T = _arr(T)
    return _scalar_or_array(f.rho_a + f.rho_b * T, T)


def dynamic_viscosity(T, fluid: Fluid | str = SODIUM):
    """Dynamic viscosity of the saturated liquid in Pa s."""
    f = get_fluid(fluid)
    T = _arr(T)
    return _scalar_or_array(f.mu_a * np.exp(f.mu_b / T), T)


def surface_tension(T, fluid: Fluid | str = SODIUM):
    """Surface tension in N m^-1."""
    f = get_fluid(fluid)
    T = _arr(T)
    return _scalar_or_array(np.maximum(f.sigma_a + f.sigma_b * T, 1.0e-3), T)


def specific_heat_liquid(T, fluid: Fluid | str = SODIUM):
    """Isobaric specific heat of the liquid in J kg^-1 K^-1."""
    f = get_fluid(fluid)
    T = _arr(T)
    return _scalar_or_array(f.cp_a + f.cp_b * (T - 800.0), T)


def thermal_conductivity_liquid(T, fluid: Fluid | str = SODIUM):
    """Thermal conductivity of the liquid in W m^-1 K^-1."""
    f = get_fluid(fluid)
    T = _arr(T)
    return _scalar_or_array(f.k_a + f.k_b * T, T)


def vapor_density(T, fluid: Fluid | str = SODIUM, pressure=None):
    """Density of saturated vapour in kg m^-3 (ideal gas)."""
    f = get_fluid(fluid)
    T = _arr(T)
    p = saturation_pressure(T, f) if pressure is None else _arr(pressure)
    return _scalar_or_array(p * f.molar_mass / (R_UNIVERSAL * T), T)


def vapor_viscosity(T, fluid: Fluid | str = SODIUM):
    """Approximate viscosity of alkali-metal vapour in Pa s."""
    T = _arr(T)
    return _scalar_or_array(1.0e-5 * (T / 1000.0) ** 0.7, T)


def vapor_thermal_conductivity(T, fluid: Fluid | str = SODIUM):
    """Approximate thermal conductivity of alkali-metal vapour in W m^-1 K^-1."""
    T = _arr(T)
    return _scalar_or_array(0.010 * (T / 1000.0) ** 0.8, T)


def vapor_specific_heat(T, fluid: Fluid | str = SODIUM):
    """Isobaric specific heat of the monatomic vapour in J kg^-1 K^-1."""
    f = get_fluid(fluid)
    return f.gamma / (f.gamma - 1.0) * f.R_specific


def figure_of_merit(T, fluid: Fluid | str = SODIUM):
    """Heat-pipe figure of merit sigma rho_l h_fg / mu_l, in W m^-2."""
    f = get_fluid(fluid)
    return (
        surface_tension(T, f)
        * liquid_density(T, f)
        * latent_heat(T, f)
        / dynamic_viscosity(T, f)
    )


def saturation_temperature(P, fluid: Fluid | str = SODIUM):
    """Invert saturation_pressure with a bracketed root finder."""
    from scipy.optimize import brentq

    f = get_fluid(fluid)
    P = _arr(P)
    scalar = np.ndim(P) == 0
    P = np.atleast_1d(P)
    if np.any(P <= 0):
        raise ValueError("Pressure must be positive")
    out = np.empty_like(P, dtype=float)
    for i, p in enumerate(P):
        lo, hi = f.T_melt + 0.5, 2200.0
        out[i] = brentq(
            lambda t: float(np.log(saturation_pressure(t, f))) - np.log(p),
            lo, hi, xtol=1e-8, rtol=1e-10, maxiter=200,
        )
    return float(out[0]) if scalar else out


def properties_at(T, fluid: Fluid | str = SODIUM) -> dict:
    """Convenience bundle of all saturated-liquid/vapour properties at T."""
    f = get_fluid(fluid)
    T = float(T)
    return {
        "fluid": f.symbol,
        "T": T,
        "p_sat": float(saturation_pressure(T, f)),
        "h_fg": float(latent_heat(T, f)),
        "rho_l": float(liquid_density(T, f)),
        "rho_v": float(vapor_density(T, f)),
        "mu_l": float(dynamic_viscosity(T, f)),
        "mu_v": float(vapor_viscosity(T, f)),
        "sigma": float(surface_tension(T, f)),
        "cp_l": float(specific_heat_liquid(T, f)),
        "cp_v": float(vapor_specific_heat(T, f)),
        "k_l": float(thermal_conductivity_liquid(T, f)),
        "k_v": float(vapor_thermal_conductivity(T, f)),
        "figure_of_merit": float(figure_of_merit(T, f)),
    }


# ---------------------------------------------------------------------------
# NaK-78 eutectic (22 wt% Na, 78 wt% K), used by the original ISR coolant loop
# ---------------------------------------------------------------------------

NAK_NA_MASS_FRACTION = 0.22
NAK_K_MASS_FRACTION = 0.78


def nak_density(T):
    """Density of NaK-78 in kg m^-3 using the ideal-mixture specific volume."""
    rho_na = liquid_density(T, SODIUM)
    rho_k = liquid_density(T, POTASSIUM)
    return 1.0 / (NAK_NA_MASS_FRACTION / rho_na + NAK_K_MASS_FRACTION / rho_k)


def nak_specific_heat(T):
    """Mass-weighted specific heat of NaK-78 in J kg^-1 K^-1."""
    return (
        NAK_NA_MASS_FRACTION * specific_heat_liquid(T, SODIUM)
        + NAK_K_MASS_FRACTION * specific_heat_liquid(T, POTASSIUM)
    )


def nak_thermal_conductivity(T):
    """Mass-weighted thermal conductivity of NaK-78 in W m^-1 K^-1."""
    return (
        NAK_NA_MASS_FRACTION * thermal_conductivity_liquid(T, SODIUM)
        + NAK_K_MASS_FRACTION * thermal_conductivity_liquid(T, POTASSIUM)
    )


def nak_viscosity(T):
    """Approximate viscosity of NaK-78 in Pa s (mass-weighted log mixing)."""
    mu_na = dynamic_viscosity(T, SODIUM)
    mu_k = dynamic_viscosity(T, POTASSIUM)
    return float(np.exp(
        NAK_NA_MASS_FRACTION * np.log(mu_na)
        + NAK_K_MASS_FRACTION * np.log(mu_k)
    ))


def nak_boiling_point() -> float:
    """Approximate atmospheric boiling point of NaK-78 in K.

    The eutectic boils over a range; 1000 K is the value used by Khandaq et al.
    (2020) as the upper operating limit of the ISR coolant.
    """
    return 1000.0
