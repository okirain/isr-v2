"""Unit tests for the thermophysical property correlations."""

import numpy as np
import pytest

from src.thermophysical.properties import (
    POTASSIUM,
    SODIUM,
    dynamic_viscosity,
    latent_heat,
    liquid_density,
    nak_density,
    saturation_pressure,
    saturation_temperature,
    surface_tension,
    vapor_density,
)


def test_sodium_boiling_point():
    assert saturation_pressure(1156.1, SODIUM) == pytest.approx(101325.0, rel=1e-6)


def test_potassium_boiling_point():
    assert saturation_pressure(1032.0, POTASSIUM) == pytest.approx(101325.0, rel=1e-6)


def test_saturation_temperature_inverse():
    T = saturation_temperature(101325.0, SODIUM)
    assert T == pytest.approx(1156.1, abs=0.02)


def test_pressure_monotonic():
    T = np.linspace(800.0, 1300.0, 50)
    p = saturation_pressure(T, SODIUM)
    assert np.all(np.diff(p) > 0)


def test_liquid_density_range():
    assert liquid_density(800.0, SODIUM) == pytest.approx(826.0, rel=0.01)
    assert liquid_density(1200.0, SODIUM) == pytest.approx(732.0, rel=0.01)


def test_latent_heat_positive_and_decreasing():
    h1 = latent_heat(800.0, SODIUM)
    h2 = latent_heat(1200.0, SODIUM)
    assert h1 > 3.5e6
    assert h1 > h2 > 3.0e6


def test_surface_tension_positive():
    assert surface_tension(800.0, SODIUM) == pytest.approx(0.171, rel=0.01)
    assert surface_tension(1300.0, SODIUM) > 0.1


def test_viscosity_decreases_with_temperature():
    assert dynamic_viscosity(800.0, SODIUM) > dynamic_viscosity(1200.0, SODIUM)


def test_vapor_density_ideal_gas():
    T = 1073.0
    p = saturation_pressure(T, SODIUM)
    rho = vapor_density(T, SODIUM)
    assert rho == pytest.approx(p * SODIUM.molar_mass / (8.31446261815324 * T), rel=1e-12)


def test_nak_density():
    # Khandaq et al. (2020) report 0.732 g/cm3 at 850 K.
    assert nak_density(850.0) == pytest.approx(733.0, rel=0.01)
