# KRUSTY benchmark data

The CSV files in this directory are **reduced benchmark state points**, not raw
instrumentation.  They were transcribed from the open-access papers of the
KRUSTY special issue:

* D. I. Poston, M. A. Gibson, T. Godfroy, and P. R. McClure, "KRUSTY Reactor
  Design," *Nuclear Technology*, vol. 206, no. sup1, pp. S13-S30, 2020.
* D. I. Poston, M. A. Gibson, R. G. Sanchez, and P. R. McClure, "Results of the
  KRUSTY Nuclear System Test," *Nuclear Technology*, vol. 206, no. sup1,
  pp. S89-S117, 2020.
* M. A. Gibson et al., "Heat Transport and Power Conversion of the Kilopower
  Reactor Test," *Nuclear Technology*, vol. 206, no. sup1, pp. S31-S42, 2020.

The original time-series figures are not distributed as machine-readable data.
The state points below were selected where the paper gives an explicit numeric
value (text or table), so that the comparison is traceable.  The 28-hour power
and sink schedules in `krusty_benchmark.py` are reconstructions of the event
sequence described in Sec. IV-X of the system-test paper; they are not a
digitisation of the figures.

Columns
-------
case_id             short identifier
time_h              time from test start, hours
fuel_temp_c         reported fuel thermocouple temperature, degC
fuel_temp_k         same value in K
description         physical meaning of the state point
source              paper and section
include_in_envelope 1 if the point is used for the <=10 % validation metric
