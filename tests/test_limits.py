"""Unit tests for the heat-pipe operating-limit solver."""

import dataclasses

import numpy as np
import pytest

from src.isr_design.isr_core_coupling import make_geometry
from src.solver.limits import (
    HeatPipeGeometry,
    capillary_limit,
    entrainment_limit,
    limit_curve,
    operating_limits,
    required_temperature,
    sonic_limit,
)
from src.thermophysical.properties import SODIUM


@pytest.fixture
def geom():
    return make_geometry(0.030, wick_type="arterial")


def test_limit_components_positive(geom):
    r = operating_limits(1050.0, geom, SODIUM, g=0.0)
    assert r.q_max > 0
    for value in (r.capillary, r.sonic, r.entrainment, r.viscous, r.boiling):
        assert value > 0
    assert r.q_max == pytest.approx(min(
        r.capillary, r.sonic, r.entrainment, r.viscous, r.boiling))


def test_gravity_reduces_capillary_limit(geom):
    vertical = dataclasses.replace(geom, tilt_deg=90.0)
    q0 = capillary_limit(900.0, vertical, SODIUM, g=0.0)
    qg = capillary_limit(900.0, vertical, SODIUM, g=9.81)
    assert qg < q0


def test_sonic_increases_with_temperature(geom):
    assert sonic_limit(1200.0, geom, SODIUM) > sonic_limit(800.0, geom, SODIUM)


def test_entrainment_increases_with_temperature(geom):
    assert entrainment_limit(1200.0, geom, SODIUM) > entrainment_limit(800.0, geom, SODIUM)


def test_required_temperature(geom):
    q = 0.5 * operating_limits(1050.0, geom, SODIUM).q_max
    T = required_temperature(q, geom, SODIUM)
    assert 750.0 < T < 1050.0


def test_limit_curve_shapes(geom):
    T = np.linspace(800.0, 1300.0, 20)
    c = limit_curve(T, geom, SODIUM)
    assert c["q_max"].shape == T.shape
    assert np.all(c["q_max"] > 0)


def test_invalid_geometry():
    with pytest.raises(ValueError):
        HeatPipeGeometry("bad", 1.0, 0.3, 0.2, 0.6, 0.001, 0.0008, 0.0003, 1e-5)
