# ABSTRACT

The Indonesian Space Reactor (ISR) is a 500 kWth NaK-78-cooled fast reactor concept with 61 uranium mononitride fuel pins. The active electromagnetic-pump loop introduces mass, parasitic power, vibration, and a single-point failure. Sodium heat pipes offer a passive alternative, but their microgravity transient behaviour has not been openly modelled in sufficient detail.

This work develops a one-dimensional/lumped-parameter reduced-order model in Python for the thermal dynamics of sodium and potassium heat pipes between 800 and 1300 K. The model couples nodal energy balances for the fuel, wall, vapour, and sink with capillary, sonic, entrainment, viscous, and boiling limits. It is validated against the 28-hour KRUSTY test and then applied to a redesign of the ISR with 19 arterial-wick sodium heat pipes of 3.0 cm diameter.

Validation gives a MAPE of 3.13% and a maximum error of 7.20% over six quasi-steady points; the 8-hour nominal temperature is 1081.8 K versus 1073.2 K. The ISR-V2 configuration transports 26.3 kW per heat pipe with a margin of 1.85, reduces the peak fuel temperature to 1370 K (45.7% margin to 2523 K), and remains below the limit with three failed heat pipes (1795 K). The system mass is nearly neutral (536.4 kg versus 538.0 kg), while the 6 kW(e) pump and its single-point failure are eliminated.

**Keywords:** heat pipe, reduced-order model, microgravity, space reactor, KRUSTY
