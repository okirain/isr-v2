"""Unit tests for the ISR core-coupling and comparison models."""

import numpy as np
import pytest

from src.isr_design.comparison import mass_comparison
from src.isr_design.isr_core_coupling import (
    UN_MELTING_LIMIT,
    design_array,
    evaluate_reference,
    fuel_temperature_profile,
    in_channel_case,
    select_reference,
    scan_designs,
)


def test_inner_cooled_profile():
    p = fuel_temperature_profile(8200.0, 850.0, boundary="inner")
    assert p["T_min"] == pytest.approx(850.0, abs=3.0)
    assert p["T_max"] > 900.0
    assert p["T_max"] > p["T_mean"] > 850.0


def test_outer_cooled_profile():
    p = fuel_temperature_profile(8200.0, 1200.0, boundary="outer")
    assert p["T_min"] == pytest.approx(1200.0, abs=1.0)
    assert p["T_max"] > 1250.0


def test_in_channel_is_infeasible():
    case = in_channel_case()
    assert case["feasible"] is False
    assert case["margin"] < 0.2


def test_arterial_reference_feasible():
    design = design_array(19, 0.030, wick_type="arterial")
    assert design.feasible
    assert design.margin >= 1.5
    ref = evaluate_reference(design)
    assert ref.fuel["T_max"] < UN_MELTING_LIMIT
    assert ref.fuel_margin > 0


def test_failure_tolerance():
    design = design_array(19, 0.030, wick_type="arterial")
    ref = evaluate_reference(design)
    assert ref.failure_cases["within_limit"].all()
    assert ref.failure_cases["T_fuel_max_k"].is_monotonic_increasing


def test_screen_fails_at_reference_power():
    design = design_array(19, 0.030, wick_type="screen")
    assert not design.feasible


def test_scan_and_selection():
    scan = scan_designs(n_values=(19, 37), d_values=(0.02, 0.03),
                        wick_type="arterial")
    assert len(scan) == 4
    ref = select_reference(scan)
    assert ref.feasible


def test_mass_comparison_is_close():
    design = design_array(19, 0.030, wick_type="arterial")
    m = mass_comparison(design)
    assert len(m) == 2
    assert abs(m["total_mass_kg"].iloc[0] - m["total_mass_kg"].iloc[1]) < 20.0
