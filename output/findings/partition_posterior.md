# The phase-line comparisons as posterior probabilities

Produced by `analyses/86_partition_posterior.py`: 2000 posterior draws of every assemblage's class proportions (Dirichlet(counts + 1/2)), 300 alternative divisions of each kind, 20 of each compared per draw. Replaces the ensemble percentiles of analyses 74 and 83 as the reported quantity (rule 18).

## The phases of Figure 1, basin set

28 sites in 3 groups. Plug-in F_ST 0.0110.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0109 | 0.0091 to 0.0128 |
| median F_ST, same-size divisions around random centers | 0.0145 | 0.0116 to 0.0180 |
| median F_ST, same-size divisions made compact | 0.0165 | 0.0113 to 0.0226 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.20**
- P(phase lines separate the pottery better than a compact division | data) = **0.05**

Boundary excess at the lines (similarity lost across a line beyond what distance predicts, river distance), 500 posterior draws, 10 alternatives of each kind per draw:

| quantity | posterior median | 95% credible interval |
|---|---|---|
| boundary excess at the phase lines (plug-in +10.0) | +9.4 | +4.6 to +13.9 |
| median, same-size divisions around random centers | +6.9 | -0.2 to +13.6 |
| median, same-size divisions made compact | +4.6 | -8.4 to +12.2 |

- P(larger boundary excess at the phase lines than at a random-center division | data) = **0.59**
- P(larger boundary excess at the phase lines than at a compact division | data) = **0.61**

## The Parkin phase against the rest of the basin

28 sites in 2 groups. Plug-in F_ST 0.0062.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0061 | 0.0051 to 0.0074 |
| median F_ST, same-size divisions around random centers | 0.0064 | 0.0052 to 0.0091 |
| median F_ST, same-size divisions made compact | 0.0076 | 0.0053 to 0.0111 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.27**
- P(phase lines separate the pottery better than a compact division | data) = **0.00**

Boundary excess at the lines (similarity lost across a line beyond what distance predicts, river distance), 500 posterior draws, 10 alternatives of each kind per draw:

| quantity | posterior median | 95% credible interval |
|---|---|---|
| boundary excess at the phase lines (plug-in +13.7) | +13.0 | +7.8 to +17.6 |
| median, same-size divisions around random centers | +10.7 | +5.7 to +15.3 |
| median, same-size divisions made compact | +15.1 | +8.7 to +19.1 |

- P(larger boundary excess at the phase lines than at a random-center division | data) = **0.53**
- P(larger boundary excess at the phase lines than at a compact division | data) = **0.05**

Agreement (adjusted Rand index) between the Parkin-versus-rest division and the two spatial clusters the site layout supports: 1.000.

## Mainfort's (2003) phases, his table

29 sites in 5 groups. Plug-in F_ST 0.0470.

| quantity | posterior median | 95% credible interval |
|---|---|---|
| F_ST under the phase scheme | 0.0459 | 0.0413 to 0.0505 |
| median F_ST, same-size divisions around random centers | 0.0297 | 0.0252 to 0.0362 |
| median F_ST, same-size divisions made compact | 0.0315 | 0.0266 to 0.0382 |

- P(phase lines separate the pottery better than a division around random centers | data) = **0.96**
- P(phase lines separate the pottery better than a compact division | data) = **0.94**
