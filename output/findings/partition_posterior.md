# The phase-line comparisons as posterior probabilities

Produced by `analyses/86_partition_posterior.py`: 2000 posterior draws of every assemblage's class proportions (Dirichlet(counts + 1/2)), 300 alternative divisions of each kind, 20 of each compared per draw. Replaces the ensemble percentiles of analyses 74 and 83 as the reported quantity (rule 18).

## The primary analysis's phases, basin set

28 sites in 3 groups. Plug-in F_ST 0.0123.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0121 | 0.0101 to 0.0144 |
| median F_ST, same-size divisions around random centers | 0.0112 | 0.0091 to 0.0144 |
| median F_ST, same-size divisions made compact | 0.0103 | 0.0079 to 0.0127 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.57**
- P(phase lines separate the pottery better than a compact division | data) = **0.72**

Boundary excess at the lines (similarity lost across a line beyond what distance predicts, river distance), 500 posterior draws, 10 alternatives of each kind per draw:

| quantity | posterior median | 95% credible interval |
|---|---|---|
| boundary excess at the phase lines (plug-in +7.7) | +7.3 | +4.6 to +9.8 |
| median, same-size divisions around random centers | +6.3 | +2.6 to +10.2 |
| median, same-size divisions made compact | +6.3 | +3.3 to +9.1 |

- P(larger boundary excess at the phase lines than at a random-center division | data) = **0.58**
- P(larger boundary excess at the phase lines than at a compact division | data) = **0.65**

## The Parkin phase against the rest of the basin

28 sites in 2 groups. Plug-in F_ST 0.0067.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0066 | 0.0052 to 0.0083 |
| median F_ST, same-size divisions around random centers | 0.0062 | 0.0049 to 0.0081 |
| median F_ST, same-size divisions made compact | 0.0057 | 0.0042 to 0.0071 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.52**
- P(phase lines separate the pottery better than a compact division | data) = **0.94**

Boundary excess at the lines (similarity lost across a line beyond what distance predicts, river distance), 500 posterior draws, 10 alternatives of each kind per draw:

| quantity | posterior median | 95% credible interval |
|---|---|---|
| boundary excess at the phase lines (plug-in +5.7) | +5.9 | +3.5 to +8.3 |
| median, same-size divisions around random centers | +5.9 | +2.6 to +9.9 |
| median, same-size divisions made compact | +5.9 | +2.4 to +9.5 |

- P(larger boundary excess at the phase lines than at a random-center division | data) = **0.43**
- P(larger boundary excess at the phase lines than at a compact division | data) = **0.49**

Agreement (adjusted Rand index) between the Parkin-versus-rest division and the two spatial clusters the site layout supports: 0.857.

## Mainfort's (2003) phases, his table

29 sites in 5 groups. Plug-in F_ST 0.0470.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0459 | 0.0413 to 0.0505 |
| median F_ST, same-size divisions around random centers | 0.0297 | 0.0252 to 0.0362 |
| median F_ST, same-size divisions made compact | 0.0315 | 0.0266 to 0.0382 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.96**
- P(phase lines separate the pottery better than a compact division | data) = **0.94**
