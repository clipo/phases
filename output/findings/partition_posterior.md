# The phase-line comparisons as posterior probabilities

Produced by `analyses/86_partition_posterior.py`: 2000 posterior draws of every assemblage's class proportions (Dirichlet(counts + 1/2)), 300 alternative divisions of each kind, 20 of each compared per draw. Replaces the ensemble percentiles of analyses 74 and 83 as the reported quantity (rule 18).

## The phases of Figure 1, basin set

28 sites in 3 groups. Plug-in F_ST 0.0110.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0109 | 0.0091 to 0.0128 |
| median F_ST, same-size divisions around random centers | 0.0146 | 0.0116 to 0.0180 |
| median F_ST, same-size divisions made compact | 0.0162 | 0.0113 to 0.0223 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.20**
- P(phase lines separate the pottery better than a compact division | data) = **0.05**

Boundary excess at the lines (similarity lost across a line beyond what distance predicts, river distance), 500 posterior draws, 10 alternatives of each kind per draw:

| quantity | posterior median | 95% credible interval |
|---|---|---|
| boundary excess at the phase lines (plug-in +10.7) | +9.9 | +5.1 to +14.8 |
| median, same-size divisions around random centers | +7.4 | -0.3 to +14.2 |
| median, same-size divisions made compact | +5.3 | -8.3 to +13.1 |

- P(larger boundary excess at the phase lines than at a random-center division | data) = **0.59**
- P(larger boundary excess at the phase lines than at a compact division | data) = **0.59**

## Mainfort's (2003) phases, his table

29 sites in 5 groups. Plug-in F_ST 0.0470.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0459 | 0.0413 to 0.0505 |
| median F_ST, same-size divisions around random centers | 0.0266 | 0.0225 to 0.0333 |
| median F_ST, same-size divisions made compact | 0.0282 | 0.0233 to 0.0377 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.99**
- P(phase lines separate the pottery better than a compact division | data) = **1.00**
