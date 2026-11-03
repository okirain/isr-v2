# LAMPIRAN A KODE PROGRAM

Lampiran ini memuat bagian inti dari kode program yang dikembangkan. Kode lengkap, termasuk pengujian unit, data benchmark, dan skrip reproduksi, tersedia di repositori proyek pada direktori src/. Kode ditulis dalam Python 3.10+ dan menggunakan NumPy, SciPy, Matplotlib, dan pandas.

## A.1 Struktur Repositori

~~~text
ISR-v2/
├── src/
│   ├── thermophysical/
│   │   ├── __init__.py
│   │   └── properties.py
│   ├── solver/
│   │   ├── __init__.py
│   │   ├── limits.py
│   │   └── transient_model.py
│   ├── validation/
│   │   ├── __init__.py
│   │   ├── krusty_benchmark.py
│   │   └── data/
│   │       ├── krusty_statepoints.csv
│   │       └── README.md
│   └── isr_design/
│       ├── __init__.py
│       ├── isr_core_coupling.py
│       └── comparison.py
├── tests/
│   ├── conftest.py
│   ├── test_properties.py
│   ├── test_limits.py
│   ├── test_transient.py
│   └── test_isr.py
├── scripts/
│   └── run_all.py
├── results/
├── thesis_md/
├── requirements.txt
└── README.md
~~~

## A.2 Properti Termofisika dan Tekanan Saturasi

~~~python
R_UNIVERSAL = 8.31446261815324
P_STANDARD = 101325.0

@dataclass(frozen=True)
class Fluid:
    name: str
    symbol: str
    molar_mass: float
    T_melt: float
    T_boil: float
    h_fg_boil: float
    dh_fg_dT: float
    rho_a: float
    rho_b: float
    mu_a: float
    mu_b: float
    sigma_a: float
    sigma_b: float
    cp_a: float
    cp_b: float
    k_a: float
    k_b: float
    gamma: float = 5.0 / 3.0

    @property
    def R_specific(self) -> float:
        return R_UNIVERSAL / self.molar_mass


def saturation_pressure(T, fluid=SODIUM):
    f = get_fluid(fluid)
    T = np.asarray(T, dtype=float)
    Rv = f.R_specific
    h0, a, Tb = f.h_fg_boil, f.dh_fg_dT, f.T_boil
    ln_p = (
        np.log(P_STANDARD)
        - (h0 - a * Tb) / Rv * (1.0 / T - 1.0 / Tb)
        + a / Rv * np.log(T / Tb)
    )
    return np.exp(ln_p)
~~~

## A.3 Batas Kapiler

~~~python
def capillary_limit(T, geom, fluid=SODIUM, g=0.0, return_components=False):
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
    q = max(0.0, (dp_cap - dp_g) / (coef_l + coef_v))
    if return_components:
        return q, {
            "dp_capillary": dp_cap,
            "dp_liquid": coef_l * q,
            "dp_vapor": coef_v * q,
            "dp_gravity": dp_g,
        }
    return q
~~~

## A.4 Right-Hand Side Model Transien

~~~python
def _rhs(self, t, y, power_func, sink_func):
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
    power_scale = max(0.05, 1.0 + c.pcs_power_gain
                      * (P_filt / c.pcs_power_reference - 1.0))
    sink_scale = sink_state * power_scale

    act = self._activation(T_v)
    h_ve_res = c.residual_vapor_conductance / max(self.n_evap, 1)
    h_vc_res = c.residual_vapor_conductance / max(self.n_cond, 1)
    h_ve_2p = self.h_ve * act
    h_vc_2p = self.h_vc * act

    q_fw = self.h_fw * (T_f - T_w)
    loss_scale = 1.0 + c.passive_loss_boost * (1.0 - min(sink_state, 1.0))
    sigma = 5.670374419e-8
    q_loss_fuel = loss_scale * self.eps_area_fuel * sigma * (T_f ** 4 - c.T_ambient ** 4)
    q_loss_wall = loss_scale * self.eps_area_wall * sigma * (T_w ** 4 - c.T_ambient ** 4)

    q_ve = h_ve_res * (T_w - T_v) + h_ve_2p * (T_w - T_v)
    q_vc = h_vc_res * (T_v - T_w) + h_vc_2p * (T_v - T_w)
    q_sink = self.h_snk * sink_scale * (T_w - T_s)

    q_hp = float((h_ve_2p * (T_w - T_v)).sum())
    if c.enable_limits:
        lim = operating_limits(T_v, self.geom, self.fluid, g=c.gravity,
                               q_operating=q_hp)
        if q_hp > lim.q_max > 0.0:
            scale = lim.q_max / q_hp
            q_ve = h_ve_res * (T_w - T_v) + h_ve_2p * (T_w - T_v) * scale
            q_vc = h_vc_res * (T_v - T_w) + h_vc_2p * (T_v - T_w) * scale

    dT_f = (P - q_loss_fuel - q_fw.sum()) / self.C_f

    dT_w = np.zeros(n)
    dT_w += q_fw - q_ve + q_vc - q_sink - q_loss_wall
    q_ax = np.zeros(n)
    q_ax[1:] += self.G_ax * (T_w[:-1] - T_w[1:])
    q_ax[:-1] += self.G_ax * (T_w[1:] - T_w[:-1])
    dT_w += q_ax
    dT_w /= self.C_w

    dT_v = (q_ve.sum() - q_vc.sum()) / self.C_v
    dT_s = (q_sink.sum() - self.ua_rej * sink_scale * (T_s - c.T_reject)) / self.C_s
    dP_filt = (P_total - P_filt) / max(c.pcs_time_constant, 1.0)

    return np.concatenate(([dT_f], dT_w, [dT_v], [dT_s], [dP_filt]))
~~~

## A.5 Solver Konduksi Radial Bahan Bakar

~~~python
def fuel_temperature_profile(q_pin, T_boundary, pin=FuelPin(),
                             boundary="outer", n=160, n_iter=25):
    ri, ro = pin.r_inner, pin.r_outer
    faces = np.linspace(ri, ro, n + 1)
    r = 0.5 * (faces[:-1] + faces[1:])
    volumes = pi * (faces[1:] ** 2 - faces[:-1] ** 2) * pin.length
    q_vol = q_pin / pin.volume
    q_nodes = q_vol * volumes

    T = np.full(n, float(T_boundary))
    for _ in range(n_iter):
        k = un_thermal_conductivity(T)
        g = np.array([
            _face_conductance(faces[i], faces[i + 1],
                              0.5 * (k[i] + k[i + 1]), pin.length)
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
    return {"r": r, "T": T, "T_max": float(T.max()),
            "T_min": float(T.min()),
            "T_mean": float(np.average(T, weights=volumes))}
~~~

## A.6 Menjalankan Seluruh Analisis

~~~bash
python -m venv --system-site-packages .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m pytest -q
python scripts/run_all.py
~~~

Perintah terakhir akan menghasilkan seluruh berkas CSV, JSON, dan gambar di direktori results/ dan thesis_md/figures/.
