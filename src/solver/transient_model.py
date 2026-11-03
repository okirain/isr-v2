"""One-dimensional transient thermal network for a liquid-metal heat pipe.

The model resolves the axial wall/vapour domain with a finite-volume mesh and
couples it to a lumped fuel node and a lumped heat-rejection node. The vapour
is treated as a single isothermal node, which is the standard and physically
appropriate reduced-order assumption for an alkali-metal heat pipe: the vapour
core is three to four orders of magnitude lighter than the wall and reaches
quasi-steady state in seconds.

Governing equations
-------------------
Fuel:
    C_f dT_f/dt = P - Q_MLI - sum_i h_fw,i (T_f - T_w,i)

Wall node i:
    C_w dT_w,i/dt = h_fw,i (T_f - T_w,i)
                    + G (T_w,i-1 - 2 T_w,i + T_w,i+1)
                    - Q_v,i - Q_sink,i - Q_loss,i

Vapour:
    C_v dT_v/dt = sum_ev Q_v,i - sum_cd Q_v,i

Heat rejection:
    C_s dT_s/dt = sum_cd Q_sink,i - UA_rej (T_s - T_rej)

The evaporator and condenser vapour conductances are multiplied by a smooth
start-up activation so that a frozen pipe conducts only through its wall. When
the operating limits are enabled the transported power is capped by the minimum
of the capillary, sonic, entrainment, viscous and boiling limits.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from ..thermophysical.properties import (
    Fluid,
    SODIUM,
    get_fluid,
    liquid_density,
    specific_heat_liquid,
)
from .limits import HeatPipeGeometry, operating_limits


@dataclass
class TransientConfig:
    """Configuration of the reduced-order transient model (all per pipe)."""

    geometry: HeatPipeGeometry
    fluid: Fluid = SODIUM
    n_pipes: int = 8
    n_axial: int = 30

    fuel_mass_total: float = 42.9       # kg, fuel + clamps (6.0 kJ/K)
    fuel_cp: float = 140.0              # J kg^-1 K^-1
    wall_mass_per_pipe: float = 0.25    # kg, Haynes 230 shell
    wall_cp: float = 500.0              # J kg^-1 K^-1
    sink_mass_total: float = 12.0       # kg, PCS/heat-rejection mass
    sink_cp: float = 500.0              # J kg^-1 K^-1

    h_fuel_wall: float = 60.0           # W K^-1 per pipe, fuel-to-wall
    h_vapor_evap: float = 40.0          # W K^-1 per pipe
    h_vapor_cond: float = 40.0          # W K^-1 per pipe
    h_sink: float = 2.0                 # W K^-1 per pipe, condenser-to-sink
    eps_area_mli: float = 0.00536       # m^2, effective radiative area of MLI
    ua_reject: float = 3.9              # W K^-1 total, radiator/rejector
    passive_loss_boost: float = 2.0     # extra PCS radiative loss when cooling is off
    pcs_power_reference: float = 2750.0 # W, nominal PCS load
    pcs_power_gain: float = 0.0         # conductance gain with PCS load
    pcs_time_constant: float = 300.0    # s, PCS thermal/control lag

    T_ambient: float = 300.0
    T_reject: float = 300.0

    gravity: float = 0.0
    enable_limits: bool = True
    startup_temperature: float = 750.0
    startup_width: float = 25.0
    vapor_capacitance: float = 5.0
    residual_vapor_conductance: float = 0.05
    initial_temperature: float = 300.0
    k_wall_axial: float = 25.0


@dataclass
class SimulationResult:
    """Time histories produced by the simulate method."""

    t: np.ndarray
    T_fuel: np.ndarray
    T_vapor: np.ndarray
    T_sink: np.ndarray
    T_wall: np.ndarray
    Q_hp: np.ndarray
    Q_limit: np.ndarray
    Q_sink: np.ndarray
    limiting_mode: list
    config: TransientConfig
    metadata: dict = field(default_factory=dict)

    @property
    def T_fuel_c(self) -> np.ndarray:
        return self.T_fuel - 273.15

    @property
    def T_vapor_c(self) -> np.ndarray:
        return self.T_vapor - 273.15

    def to_dict(self) -> dict:
        return {
            "t": self.t.tolist(),
            "T_fuel": self.T_fuel.tolist(),
            "T_vapor": self.T_vapor.tolist(),
            "T_sink": self.T_sink.tolist(),
            "Q_hp": self.Q_hp.tolist(),
            "Q_limit": self.Q_limit.tolist(),
        }


class HeatPipeTransientModel:
    """Finite-volume transient model of one equivalent heat pipe."""

    def __init__(self, config: TransientConfig):
        self.cfg = config
        self.geom = config.geometry
        self.fluid = get_fluid(config.fluid)
        self._build_mesh()

    # ------------------------------------------------------------------
    # Mesh and capacities
    # ------------------------------------------------------------------
    def _build_mesh(self) -> None:
        c = self.cfg
        g = self.geom
        n = c.n_axial
        self.n = n
        self.dz = g.L_total / n
        self.z = (np.arange(n) + 0.5) * self.dz
        self.mask_evap = self.z < g.L_evap
        self.mask_adiab = (self.z >= g.L_evap) & (self.z < g.L_evap + g.L_adiab)
        self.mask_cond = self.z >= g.L_evap + g.L_adiab
        self.n_evap = int(self.mask_evap.sum())
        self.n_cond = int(self.mask_cond.sum())
        if self.n_evap == 0 or self.n_cond == 0:
            raise ValueError("Mesh must contain evaporator and condenser nodes")

        # Wick and working-fluid mass, used for the wall heat capacity
        rho_ni = 8900.0
        m_wick = rho_ni * (1.0 - g.porosity) * g.area_wick * g.L_total
        T_ref = 0.5 * (self.fluid.T_melt + self.fluid.T_boil)
        m_na = g.porosity * g.area_wick * g.L_total * float(
            liquid_density(T_ref, self.fluid)
        )
        cp_l = float(specific_heat_liquid(T_ref, self.fluid))
        c_wall_total = c.wall_mass_per_pipe * c.wall_cp + m_wick * 460.0 + m_na * cp_l
        self.C_w = np.full(n, c_wall_total / n)

        self.C_f = c.fuel_mass_total / c.n_pipes * c.fuel_cp
        self.C_s = c.sink_mass_total / c.n_pipes * c.sink_cp
        self.C_v = c.vapor_capacitance

        # Axial wall conduction
        a_cond = g.area_wall + 0.5 * g.area_wick
        self.G_ax = c.k_wall_axial * a_cond / self.dz

        # Distributed couplings
        self.h_fw = np.zeros(n)
        self.h_ve = np.zeros(n)
        self.h_vc = np.zeros(n)
        self.h_snk = np.zeros(n)

        self.h_fw[self.mask_evap] = c.h_fuel_wall / self.n_evap
        self.h_ve[self.mask_evap] = c.h_vapor_evap / self.n_evap
        self.h_vc[self.mask_cond] = c.h_vapor_cond / self.n_cond
        self.h_snk[self.mask_cond] = c.h_sink / self.n_cond

        # 80 % of the passive loss is assigned to the fuel, 20 % to the wall.
        # The loss is radiative, Q = eps A sigma (T^4 - T_amb^4), so the
        # effective conductance falls naturally during cooldown.
        self.eps_area_fuel = 0.8 * c.eps_area_mli / c.n_pipes
        self.eps_area_wall = np.full(n, 0.2 * c.eps_area_mli / c.n_pipes / n)
        self.ua_rej = c.ua_reject / c.n_pipes

        self._last_limits = None

    # ------------------------------------------------------------------
    # Right-hand side
    # ------------------------------------------------------------------
    def _activation(self, T_v: float) -> float:
        c = self.cfg
        return float(1.0 / (1.0 + np.exp(-(T_v - c.startup_temperature) / c.startup_width)))

    def _rhs(
        self,
        t: float,
        y: np.ndarray,
        power_func: Callable[[float], float],
        sink_func: Callable[[float], float] | None,
    ) -> np.ndarray:
        c = self.cfg
        n = self.n
        T_f = y[0]
        T_w = y[1:1 + n]
        T_v = y[1 + n]
        T_s = y[2 + n]
        P_filt = y[3 + n]

        P_total = max(float(power_func(t)), 0.0)
        P = P_total / c.n_pipes
        sink_state = 1.0 if sink_func is None else max(float(sink_func(t)), 0.0)
        # The PCS conductance grows with the Stirling stroke / simulator flow,
        # but the mechanical/thermal response has a first-order lag.
        power_scale = max(0.05, 1.0 + c.pcs_power_gain
                          * (P_filt / c.pcs_power_reference - 1.0))
        sink_scale = sink_state * power_scale

        act = self._activation(T_v)
        h_ve_res = c.residual_vapor_conductance / max(self.n_evap, 1)
        h_vc_res = c.residual_vapor_conductance / max(self.n_cond, 1)
        h_ve_2p = self.h_ve * act
        h_vc_2p = self.h_vc * act

        loss_scale = 1.0 + c.passive_loss_boost * (1.0 - min(sink_state, 1.0))
        sigma = 5.670374419e-8
        q_fw = self.h_fw * (T_f - T_w)
        q_loss_fuel = (loss_scale * self.eps_area_fuel * sigma
                       * (T_f ** 4 - c.T_ambient ** 4))
        q_loss_wall = (loss_scale * self.eps_area_wall * sigma
                       * (T_w ** 4 - c.T_ambient ** 4))

        # Conduction through the frozen wall/vapour space is not subject to the
        # two-phase transport limit; only the evaporating/condensing component is.
        q_ve_res = h_ve_res * (T_w - T_v)
        q_vc_res = h_vc_res * (T_v - T_w)
        q_ve_2p = h_ve_2p * (T_w - T_v)
        q_vc_2p = h_vc_2p * (T_v - T_w)
        q_ve = q_ve_res + q_ve_2p
        q_vc = q_vc_res + q_vc_2p
        q_sink = self.h_snk * sink_scale * (T_w - T_s)

        q_hp = float(q_ve_2p.sum())
        q_lim = np.inf
        if c.enable_limits:
            lim = operating_limits(T_v, self.geom, self.fluid, g=c.gravity,
                                   q_operating=q_hp)
            self._last_limits = lim
            q_lim = lim.q_max
            if q_hp > q_lim > 0.0:
                scale = q_lim / q_hp
                q_ve = q_ve_res + q_ve_2p * scale
                q_vc = q_vc_res + q_vc_2p * scale
                q_hp = q_lim

        dT_f = (P - q_loss_fuel - q_fw.sum()) / self.C_f

        dT_w = np.zeros(n)
        dT_w += q_fw
        dT_w -= q_ve
        dT_w += q_vc
        dT_w -= q_sink
        dT_w -= q_loss_wall

        # Axial conduction (adiabatic ends)
        q_ax = np.zeros(n)
        q_ax[1:] += self.G_ax * (T_w[:-1] - T_w[1:])
        q_ax[:-1] += self.G_ax * (T_w[1:] - T_w[:-1])
        dT_w += q_ax
        dT_w /= self.C_w

        dT_v = (q_ve.sum() - q_vc.sum()) / self.C_v
        dT_s = (q_sink.sum() - self.ua_rej * sink_scale * (T_s - c.T_reject)) / self.C_s
        dP_filt = (P_total - P_filt) / max(c.pcs_time_constant, 1.0)

        return np.concatenate(([dT_f], dT_w, [dT_v], [dT_s], [dP_filt]))

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    def _diagnostics(self, t: float, y: np.ndarray,
                     power_func: Callable[[float], float],
                     sink_func: Callable[[float], float] | None) -> dict:
        c = self.cfg
        n = self.n
        T_f, T_w, T_v, T_s = y[0], y[1:1 + n], y[1 + n], y[2 + n]
        P_filt = y[3 + n]
        act = self._activation(T_v)
        h_ve_res = c.residual_vapor_conductance / max(self.n_evap, 1)
        h_vc_res = c.residual_vapor_conductance / max(self.n_cond, 1)
        q_ve_res = h_ve_res * (T_w - T_v)
        q_vc_res = h_vc_res * (T_v - T_w)
        q_ve_2p = self.h_ve * act * (T_w - T_v)
        q_vc_2p = self.h_vc * act * (T_v - T_w)
        sink_state = 1.0 if sink_func is None else max(float(sink_func(t)), 0.0)
        power_scale = max(0.05, 1.0 + c.pcs_power_gain
                          * (P_filt / c.pcs_power_reference - 1.0))
        sink_scale = sink_state * power_scale
        q_sink = self.h_snk * sink_scale * (T_w - T_s)
        q_hp = float(q_ve_2p.sum())
        q_lim = np.nan
        mode = "unlimited"
        if c.enable_limits:
            lim = operating_limits(T_v, self.geom, self.fluid, g=c.gravity,
                                   q_operating=q_hp)
            q_lim = lim.q_max
            mode = lim.limiting_mode
            if q_hp > q_lim > 0.0:
                q_hp = q_lim
        return {
            "Q_hp": q_hp,
            "Q_limit": q_lim,
            "Q_sink": float(q_sink.sum()),
            "limiting_mode": mode,
            "T_vapor": float(T_v),
            "T_fuel": float(T_f),
            "T_sink": float(T_s),
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def simulate(
        self,
        t_end: float,
        power_func: Callable[[float], float],
        dt: float = 10.0,
        sink_func: Callable[[float], float] | None = None,
        method: str = "BDF",
        rtol: float = 1e-5,
        atol: float = 1e-4,
        t_eval: np.ndarray | None = None,
    ) -> SimulationResult:
        """Integrate the thermal network from a uniform initial state."""
        c = self.cfg
        n = self.n
        y0 = np.full(3 + n + 1, c.initial_temperature, dtype=float)
        y0[3 + n] = 0.0  # filtered PCS power starts at zero
        if t_eval is None:
            t_eval = np.arange(0.0, t_end, dt)
            if len(t_eval) == 0 or t_eval[-1] < t_end:
                t_eval = np.append(t_eval, t_end)
        t_eval = np.asarray(t_eval, dtype=float)
        t_eval = t_eval[(t_eval >= 0.0) & (t_eval <= t_end)]

        sol = solve_ivp(
            lambda t, y: self._rhs(t, y, power_func, sink_func),
            (0.0, float(t_end)),
            y0,
            t_eval=t_eval,
            method=method,
            rtol=rtol,
            atol=atol,
            max_step=dt,
        )
        if not sol.success:
            sol = solve_ivp(
                lambda t, y: self._rhs(t, y, power_func, sink_func),
                (0.0, float(t_end)),
                y0,
                t_eval=t_eval,
                method="Radau",
                rtol=rtol,
                atol=atol,
                max_step=dt,
            )
        if not sol.success:
            raise RuntimeError(f"Transient integration failed: {sol.message}")

        diag = [self._diagnostics(t, sol.y[:, i], power_func, sink_func)
                for i, t in enumerate(sol.t)]
        return SimulationResult(
            t=sol.t,
            T_fuel=sol.y[0],
            T_wall=sol.y[1:1 + n].T,
            T_vapor=sol.y[1 + n],
            T_sink=sol.y[2 + n],
            Q_hp=np.array([d["Q_hp"] for d in diag]),
            Q_limit=np.array([d["Q_limit"] for d in diag], dtype=float),
            Q_sink=np.array([d["Q_sink"] for d in diag]),
            limiting_mode=[d["limiting_mode"] for d in diag],
            config=c,
            metadata={"method": method, "n_axial": n, "n_pipes": c.n_pipes},
        )

    def steady_state(
        self,
        power: float,
        sink_scale: float = 1.0,
        t_end: float = 7200.0,
        dt: float = 10.0,
    ) -> SimulationResult:
        """Approximate steady state by marching a constant-power case."""
        return self.simulate(
            t_end=t_end,
            power_func=lambda t: power,
            dt=dt,
            sink_func=lambda t: sink_scale,
        )


def step_profile(segments: list[tuple[float, float]], ramp: float = 10.0):
    """Return a smooth stepwise function from (time, value) pairs.

    The transitions are linear over ramp seconds, which keeps the ODE solver
    well behaved while preserving the intended step magnitudes.
    """
    segs = sorted(segments)
    times = np.array([s[0] for s in segs], dtype=float)
    values = np.array([s[1] for s in segs], dtype=float)

    def func(t: float) -> float:
        if t <= times[0]:
            return float(values[0])
        if t >= times[-1]:
            return float(values[-1])
        i = int(np.searchsorted(times, t, side="right") - 1)
        t0, t1 = times[i], times[i + 1]
        v0, v1 = values[i], values[i + 1]
        if ramp <= 0:
            return float(v0 if t < t1 else v1)
        if (t1 - t0) <= ramp:
            # The ramp spans the whole interval.
            return float(v0 + (v1 - v0) * (t - t0) / (t1 - t0))
        # Hold v0 until the last 'ramp' seconds of the interval, then blend.
        if t >= t1 - ramp:
            return float(v1 - (v1 - v0) * (t1 - t) / ramp)
        return float(v0)

    return func
