# ISR-V2: Reduced-Order Space Heat-Pipe Simulator

**Author:** Muhammad Naufal Tamimi (20/456315/TK/50445)  
**Department:** Teknik Nuklir dan Teknik Fisika, Fakultas Teknik, Universitas Gadjah Mada  
**Thesis:** *Pengembangan Model Reduced-Order untuk Simulasi Transien Heat Pipe Logam Cair pada Kondisi Mikrogravitasi dan Aplikasinya pada Desain Ulang Sistem Pendingin Indonesian Space Reactor (ISR-V2)*

This repository contains a fast, verified 1D/lumped-parameter reduced-order
model (ROM) of high-temperature sodium and potassium heat pipes operating in
microgravity, its validation against the NASA Kilopower/KRUSTY 28-hour test,
and its application to a passive heat-pipe redesign of the Indonesian Space
Reactor (ISR).

## Results at a glance

| Quantity | Result |
|----------|--------|
| KRUSTY validation, six quasi-steady points | MAPE 3.13%, max error 7.20% (all points 9.59%) |
| KRUSTY nominal 8 h fuel temperature | 1081.8 K (measured 1073.15 K, +0.81%) |
| KRUSTY time to 500 °C | 0.70 h (measured 0.63 h, +4.2 min) |
| ISR-V2 reference array | 19 sodium heat pipes, 3.0 cm OD, arterial wick |
| Power per heat pipe | 26.3 kW, transport margin 1.85 |
| Peak UN fuel temperature | 1369.8 K (limit 2523 K, margin 45.7%) |
| Three adjacent heat pipes failed | 1794.9 K, still below the limit |
| Cooling-system mass | 89.29 kg (pumped NaK-78 loop: 90.85 kg) |
| Pump power eliminated | 6 kW(e) |

## What is implemented

### Thermophysical properties (src/thermophysical/properties.py)

Smooth correlations for saturated Na and K from 800 to 1300 K:

* saturation pressure from a temperature-dependent Clausius-Clapeyron fit
* latent heat, liquid density, dynamic viscosity, surface tension
* liquid specific heat and thermal conductivity
* vapour density, viscosity, conductivity and specific heat
* NaK-78 mixture properties for the baseline coolant loop

### Operating limits (src/solver/limits.py)

* capillary limit with Darcy wick flow, Poiseuille vapour flow and gravity
* sonic (choked-vapour) limit
* entrainment (Weber) limit
* viscous limit
* wick boiling limit
* combined limiter and inverse temperature solver

### Transient model (src/solver/transient_model.py)

* finite-volume axial wall mesh
* lumped fuel, vapour and heat-rejection nodes
* radiative MLI loss and PCS load-following conductance
* smooth frozen-start / thaw activation
* BDF integration with Radau fallback
* per-pipe or whole-core scaling

### Validation (src/validation/krusty_benchmark.py)

* KRUSTY Haynes-230/Na heat-pipe geometry
* reconstructed 28-hour power and sink schedules
* six-point quasi-steady validation envelope plus cooldown
* MAPE, RMSE and maximum-error metrics
* publication-quality figure

### ISR design (src/isr_design/)

* 1D radial UN fuel conduction with temperature-dependent conductivity
* heat-pipe array design scan over count, diameter and wick type
* in-channel one-for-one replacement study
* core-block and reflector mass penalties
* failure analysis for one and three adjacent failed pipes
* pumped-loop versus heat-pipe mass comparison

## Repository layout

~~~text
ISR-v2/
├── src/
│   ├── thermophysical/properties.py
│   ├── solver/limits.py
│   ├── solver/transient_model.py
│   ├── validation/krusty_benchmark.py
│   ├── validation/data/krusty_statepoints.csv
│   └── isr_design/
│       ├── isr_core_coupling.py
│       └── comparison.py
├── tests/
├── scripts/run_all.py
├── results/
├── thesis_md/
├── requirements.txt
└── README.md
~~~

## Installation

Python 3.10 or newer is required. The code is tested with CPython 3.12.

~~~bash
cd ISR-v2
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
~~~

If the system Python already provides NumPy, SciPy and Matplotlib, a virtual
environment with the option --system-site-packages also works.

## Usage

Run the unit tests:

~~~bash
python -m pytest -q
~~~

Regenerate every result and figure used in the thesis:

~~~bash
python scripts/run_all.py
~~~

Run only the KRUSTY validation:

~~~bash
python -m src.validation.krusty_benchmark
~~~

Run only the ISR-V2 design and comparison:

~~~bash
python -m src.isr_design.comparison
~~~

The scripts write CSV and JSON files to results/ and figures to
thesis_md/figures/.

## Validation notes

The KRUSTY state points are transcribed from the open-access papers of the
KRUSTY special issue (Poston et al., *Nuclear Technology*, vol. 206, sup1,
2020). The original time-series figures are not distributed as machine-readable
data, so the benchmark uses numeric values stated explicitly in the text. The
28-hour power and sink schedules are reconstructions of the documented event
sequence, not a digitisation of the figures.

KRUSTY's heat pipes operated as gravity-assisted thermosiphons. The validation
therefore targets the thermal response and energy balance, not the microgravity
capillary limit. The microgravity limit model is validated analytically and
through parametric studies.

The model is thermal-only. Fission power is a prescribed boundary condition;
reactivity feedback and point kinetics are outside the scope of this work. This
is why the load-following power oscillations of KRUSTY are not reproduced.

## Thesis documents

The complete Indonesian thesis is in thesis_md/:

* front matter: title, plagiarism statement, approval page, preface
* table of contents, lists of tables, figures and abbreviations
* abstract (English) and intisari (Indonesian)
* Chapters I–VI
* bibliography in IEEE style
* Appendix A (program code) and Appendix B (material tables)

All figures referenced by the thesis are in thesis_md/figures/.

## Key caveats

* The reference ISR-V2 design assumes an arterial wick with permeability
  1.0e-9 m2. A conventional screen wick cannot transport 26.3 kW in a 3.0 cm
  pipe; this is the main technology-development item.
* The mass comparison is an engineering estimate with 10–15% contingency, not
  a detailed structural design.
* The fuel temperature is a 1D radial estimate with a spreading resistance in
  the Mo-30Re core block. It is not a 3D conjugate heat-transfer solution.
* The design is not coupled to neutronics. The 500 kWth power and uniform pin
  distribution are inputs, not a criticality calculation.
