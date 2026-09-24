# Do the published phases behave as groups? (placebo 11)

PLACEBO 11: the groups are NOT the phases but a division with the phases' sizes around random centers (adjusted Rand with the phases 0.411); 'phase' below means this division.

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the placebo division (Kent 8, Parkin 11, Walls 9).

Observed: between-phase cultural F_ST **0.0138**; boundary excess at phase lines **+11.1**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0094 [0.0017, 0.0323] | 27% | +5.1 [-2.4, +15.4] | 13% | yes |
| 0.03 | 0.1 | 0.0098 [0.0019, 0.1266] | 41% | +5.8 [-2.5, +49.8] | 28% | yes |
| 0.03 | 0.2 | 0.0214 [0.0023, 0.1639] | 62% | +13.2 [-1.3, +62.9] | 53% | yes |
| 0.03 | 0.3 | 0.0499 [0.0030, 0.1718] | 79% | +26.1 [+1.0, +72.7] | 74% | no |
| 0.1 | 0 | 0.0058 [0.0013, 0.0245] | 13% | +2.7 [-2.0, +11.8] | 3% | yes |
| 0.1 | 0.1 | 0.0080 [0.0016, 0.0765] | 32% | +3.9 [-2.1, +28.6] | 19% | yes |
| 0.1 | 0.2 | 0.0133 [0.0016, 0.1141] | 48% | +6.8 [-1.9, +37.3] | 36% | yes |
| 0.1 | 0.3 | 0.0282 [0.0028, 0.1110] | 70% | +16.1 [-0.2, +43.5] | 58% | yes |
| 0.5 | 0 | 0.0034 [0.0007, 0.0149] | 3% | +1.6 [-3.2, +7.0] | 0% | yes |
| 0.5 | 0.1 | 0.0044 [0.0008, 0.0344] | 10% | +1.7 [-2.9, +11.1] | 3% | yes |
| 0.5 | 0.2 | 0.0075 [0.0011, 0.0481] | 30% | +3.5 [-1.8, +17.3] | 12% | yes |
| 0.5 | 0.3 | 0.0129 [0.0013, 0.0525] | 47% | +6.1 [-1.6, +21.5] | 24% | no |
| 1 | 0 | 0.0031 [0.0008, 0.0137] | 3% | +1.1 [-3.4, +6.5] | 0% | yes |
| 1 | 0.1 | 0.0040 [0.0006, 0.0245] | 8% | +1.7 [-2.6, +10.7] | 2% | yes |
| 1 | 0.2 | 0.0059 [0.0008, 0.0302] | 21% | +2.6 [-2.5, +13.4] | 5% | yes |
| 1 | 0.3 | 0.0083 [0.0009, 0.0396] | 38% | +4.3 [-2.2, +18.1] | 12% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0062**; boundary excess at the Parkin line **+13.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0045 [0.0004, 0.0234] | 38% | +3.6 [-5.1, +14.7] | 4% |
| 0.03 | 0.1 | 0.0047 [0.0005, 0.0882] | 41% | +3.9 [-5.0, +31.2] | 14% |
| 0.03 | 0.2 | 0.0094 [0.0007, 0.1064] | 63% | +6.7 [-5.8, +42.4] | 32% |
| 0.03 | 0.3 | 0.0255 [0.0010, 0.1283] | 76% | +13.7 [-4.2, +49.8] | 50% |
| 0.1 | 0 | 0.0035 [0.0003, 0.0191] | 32% | +1.8 [-6.2, +12.9] | 2% |
| 0.1 | 0.1 | 0.0049 [0.0006, 0.0673] | 43% | +2.4 [-5.9, +24.6] | 9% |
| 0.1 | 0.2 | 0.0085 [0.0003, 0.0976] | 56% | +4.7 [-5.8, +33.0] | 25% |
| 0.1 | 0.3 | 0.0217 [0.0013, 0.1014] | 75% | +10.5 [-4.5, +38.4] | 41% |
| 0.5 | 0 | 0.0021 [0.0003, 0.0156] | 14% | +0.8 [-6.4, +9.3] | 0% |
| 0.5 | 0.1 | 0.0032 [0.0003, 0.0349] | 23% | +1.0 [-6.9, +15.6] | 3% |
| 0.5 | 0.2 | 0.0063 [0.0004, 0.0494] | 50% | +3.1 [-5.6, +20.7] | 11% |
| 0.5 | 0.3 | 0.0121 [0.0006, 0.0544] | 67% | +6.0 [-4.9, +27.5] | 24% |
| 1 | 0 | 0.0019 [0.0002, 0.0110] | 12% | +0.5 [-6.5, +7.5] | 0% |
| 1 | 0.1 | 0.0026 [0.0002, 0.0264] | 20% | +1.2 [-5.3, +16.5] | 5% |
| 1 | 0.2 | 0.0045 [0.0003, 0.0309] | 42% | +3.0 [-6.2, +19.9] | 8% |
| 1 | 0.3 | 0.0077 [0.0004, 0.0395] | 53% | +4.5 [-5.0, +23.2] | 14% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0138 against a median of 0.0031 (95 percent range 0.0008 to 0.0137); 3% of runs reach it. Boundary excess at phase lines: observed +11.1 against +1.1 (-3.4 to +6.5); 0% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.87** (prior 0.75).
- P(local innovation, share > 0) = **0.75** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.92 | 0.81 |
| 1% | phase F_ST and boundary excess only | 0.85 | 0.77 |
| 2% | all six summaries | 0.88 | 0.78 |
| 2% | phase F_ST and boundary excess only | 0.86 | 0.80 |
| 5% | all six summaries | 0.87 | 0.75 |
| 5% | phase F_ST and boundary excess only | 0.86 | 0.75 |
| 10% | all six summaries | 0.81 | 0.73 |
| 10% | phase F_ST and boundary excess only | 0.85 | 0.77 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.375 | 0.985 |
| 0.5 | 0.269 | 0.971 |
| 0.1 | 0.107 | 0.869 |
| 0.03 | 0.050 | 0.667 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.37 | | 0 | 0.25 |
| 0.1 | 0.31 | | 0.1 | 0.29 |
| 0.5 | 0.19 | | 0.2 | 0.25 |
| 1 | 0.13 | | 0.3 | 0.21 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.79 | 0.65 |
| 0.1 | 0.2 | 0.79 | 0.64 |
| 0.03 | 0 | 0.84 | 0.70 |
| 0.03 | 0.2 | 0.86 | 0.68 |
| 1 (no boundary) | 0 | 0.69 | 0.57 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
