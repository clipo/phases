# Do the published phases behave as groups? (placebo 11)

PLACEBO 11: the groups are NOT the phases but a division with the phases' sizes around random centers (adjusted Rand with the phases 0.264); 'phase' below means this division.

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the placebo division (Kent 11, Parkin 12, Walls 5).

Observed: between-phase cultural F_ST **0.0080**; boundary excess at phase lines **+5.1**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0075 [0.0015, 0.0222] | 44% | +4.8 [-2.3, +14.2] | 47% | yes |
| 0.03 | 0.1 | 0.0084 [0.0015, 0.0839] | 51% | +6.1 [-1.4, +41.3] | 56% | yes |
| 0.03 | 0.2 | 0.0187 [0.0021, 0.1290] | 70% | +13.7 [-0.8, +53.5] | 74% | yes |
| 0.03 | 0.3 | 0.0408 [0.0030, 0.1343] | 84% | +24.7 [-0.1, +62.5] | 89% | yes |
| 0.1 | 0 | 0.0043 [0.0009, 0.0150] | 19% | +2.5 [-3.3, +10.1] | 24% | yes |
| 0.1 | 0.1 | 0.0048 [0.0008, 0.0594] | 30% | +3.4 [-2.7, +26.8] | 37% | yes |
| 0.1 | 0.2 | 0.0092 [0.0011, 0.0650] | 57% | +6.0 [-2.5, +30.2] | 55% | yes |
| 0.1 | 0.3 | 0.0214 [0.0017, 0.0800] | 72% | +12.6 [-2.1, +37.2] | 75% | no |
| 0.5 | 0 | 0.0024 [0.0005, 0.0092] | 5% | +0.8 [-4.4, +6.3] | 6% | yes |
| 0.5 | 0.1 | 0.0030 [0.0006, 0.0193] | 13% | +0.9 [-4.8, +9.2] | 11% | yes |
| 0.5 | 0.2 | 0.0049 [0.0009, 0.0340] | 35% | +2.9 [-2.9, +13.5] | 32% | yes |
| 0.5 | 0.3 | 0.0089 [0.0010, 0.0355] | 54% | +4.4 [-3.0, +17.2] | 44% | no |
| 1 | 0 | 0.0021 [0.0005, 0.0095] | 6% | +0.6 [-3.5, +5.3] | 4% | yes |
| 1 | 0.1 | 0.0026 [0.0005, 0.0188] | 11% | +1.2 [-3.7, +8.3] | 11% | yes |
| 1 | 0.2 | 0.0050 [0.0007, 0.0252] | 34% | +2.3 [-4.2, +11.3] | 24% | yes |
| 1 | 0.3 | 0.0065 [0.0008, 0.0315] | 44% | +2.9 [-2.7, +13.8] | 34% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0067**; boundary excess at the Parkin line **+5.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0028 [0.0003, 0.0146] | 19% | +1.4 [-3.2, +8.1] | 8% |
| 0.03 | 0.1 | 0.0032 [0.0004, 0.0463] | 30% | +2.0 [-3.0, +23.3] | 21% |
| 0.03 | 0.2 | 0.0084 [0.0006, 0.0853] | 56% | +5.3 [-3.1, +32.7] | 48% |
| 0.03 | 0.3 | 0.0171 [0.0006, 0.0805] | 71% | +10.0 [-1.7, +36.3] | 62% |
| 0.1 | 0 | 0.0026 [0.0003, 0.0127] | 14% | +1.0 [-3.9, +7.4] | 8% |
| 0.1 | 0.1 | 0.0031 [0.0003, 0.0429] | 25% | +1.2 [-3.1, +16.7] | 16% |
| 0.1 | 0.2 | 0.0065 [0.0002, 0.0479] | 48% | +3.2 [-3.2, +22.3] | 34% |
| 0.1 | 0.3 | 0.0137 [0.0007, 0.0615] | 66% | +6.2 [-2.0, +25.7] | 53% |
| 0.5 | 0 | 0.0021 [0.0002, 0.0119] | 13% | +0.3 [-3.8, +6.1] | 3% |
| 0.5 | 0.1 | 0.0026 [0.0003, 0.0187] | 19% | +0.3 [-4.0, +8.2] | 8% |
| 0.5 | 0.2 | 0.0045 [0.0005, 0.0336] | 36% | +1.7 [-3.8, +11.9] | 18% |
| 0.5 | 0.3 | 0.0080 [0.0006, 0.0386] | 56% | +3.0 [-2.8, +16.1] | 28% |
| 1 | 0 | 0.0017 [0.0002, 0.0112] | 7% | +0.0 [-3.9, +5.9] | 3% |
| 1 | 0.1 | 0.0024 [0.0002, 0.0193] | 16% | +0.1 [-3.3, +8.6] | 7% |
| 1 | 0.2 | 0.0044 [0.0003, 0.0267] | 38% | +1.3 [-3.8, +10.2] | 15% |
| 1 | 0.3 | 0.0062 [0.0004, 0.0338] | 48% | +2.4 [-3.6, +14.0] | 25% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0080 against a median of 0.0021 (95 percent range 0.0005 to 0.0095); 6% of runs reach it. Boundary excess at phase lines: observed +5.1 against +0.6 (-3.5 to +5.3); 4% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.75** (prior 0.75).
- P(local innovation, share > 0) = **0.70** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.73 | 0.62 |
| 1% | phase F_ST and boundary excess only | 0.79 | 0.65 |
| 2% | all six summaries | 0.76 | 0.67 |
| 2% | phase F_ST and boundary excess only | 0.78 | 0.67 |
| 5% | all six summaries | 0.75 | 0.70 |
| 5% | phase F_ST and boundary excess only | 0.80 | 0.71 |
| 10% | all six summaries | 0.75 | 0.69 |
| 10% | phase F_ST and boundary excess only | 0.78 | 0.72 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.378 | 0.993 |
| 0.5 | 0.264 | 0.987 |
| 0.1 | 0.101 | 0.937 |
| 0.03 | 0.052 | 0.816 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.24 | | 0 | 0.30 |
| 0.1 | 0.27 | | 0.1 | 0.32 |
| 0.5 | 0.25 | | 0.2 | 0.23 |
| 1 | 0.25 | | 0.3 | 0.14 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.76 | 0.68 |
| 0.1 | 0.2 | 0.75 | 0.62 |
| 0.03 | 0 | 0.83 | 0.72 |
| 0.03 | 0.2 | 0.95 | 0.77 |
| 1 (no boundary) | 0 | 0.69 | 0.57 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
