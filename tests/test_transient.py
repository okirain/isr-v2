"""Unit tests for the transient thermal network."""

import numpy as np
import pytest

from src.isr_design.isr_core_coupling import make_geometry
from src.solver.transient_model import (
    HeatPipeTransientModel,
    TransientConfig,
    step_profile,
)


def test_step_profile_hold_and_ramp():
    p = step_profile([(0.0, 0.0), (100.0, 10.0)], ramp=10.0)
    assert p(50.0) == pytest.approx(0.0)
    assert p(95.0) == pytest.approx(5.0)
    assert p(100.0) == pytest.approx(10.0)


def test_step_profile_full_interval_ramp():
    p = step_profile([(0.0, 0.0), (100.0, 10.0), (200.0, 10.0)], ramp=100.0)
    assert p(50.0) == pytest.approx(5.0)
    assert p(100.0) == pytest.approx(10.0)


def test_model_heats_and_transports():
    geom = make_geometry(0.030, wick_type="arterial")
    cfg = TransientConfig(geometry=geom, n_pipes=1, n_axial=20,
                          fuel_mass_total=8.0, sink_mass_total=20.0,
                          h_sink=100.0, ua_reject=60.0, eps_area_mli=0.02,
                          residual_vapor_conductance=1.0,
                          pcs_power_gain=0.0)
    model = HeatPipeTransientModel(cfg)
    res = model.simulate(4 * 3600.0, lambda t: 20000.0, dt=30.0)
    assert np.all(np.isfinite(res.T_fuel))
    assert res.T_fuel[-1] > 600.0
    assert res.Q_hp[-1] > 1000.0


def test_limits_cap_transport():
    geom = make_geometry(0.0127, wick_type="screen")
    cfg = TransientConfig(geometry=geom, n_pipes=1, n_axial=20,
                          fuel_mass_total=8.0, sink_mass_total=20.0,
                          h_sink=100.0, ua_reject=60.0, eps_area_mli=0.02,
                          residual_vapor_conductance=1.0,
                          pcs_power_gain=0.0, enable_limits=True)
    res = HeatPipeTransientModel(cfg).simulate(
        2 * 3600.0, lambda t: 20000.0, dt=30.0)
    finite = np.isfinite(res.Q_limit)
    assert finite.any()
    assert np.nanmax(res.Q_hp) <= np.nanmax(res.Q_limit) + 1.0
