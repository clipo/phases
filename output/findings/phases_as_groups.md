# Do the published phases behave as groups?

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the published phases (Kent 11, Parkin 12, Walls 5).

Observed: between-phase cultural F_ST **0.0123**; boundary excess at phase lines **+7.7**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0062 [0.0013, 0.0237] | 22% | +6.5 [+0.6, +17.8] | 38% | yes |
| 0.03 | 0.1 | 0.0085 [0.0012, 0.1530] | 36% | +9.1 [+0.3, +61.9] | 56% | yes |
| 0.03 | 0.2 | 0.0215 [0.0022, 0.1827] | 58% | +20.0 [+1.3, +78.8] | 75% | yes |
| 0.03 | 0.3 | 0.0491 [0.0026, 0.1977] | 75% | +35.5 [+2.7, +93.7] | 89% | yes |
| 0.1 | 0 | 0.0052 [0.0010, 0.0213] | 11% | +4.3 [-1.3, +13.1] | 18% | yes |
| 0.1 | 0.1 | 0.0063 [0.0012, 0.1123] | 28% | +5.2 [-0.8, +38.8] | 35% | yes |
| 0.1 | 0.2 | 0.0117 [0.0014, 0.1113] | 49% | +9.4 [+0.6, +44.8] | 56% | yes |
| 0.1 | 0.3 | 0.0301 [0.0022, 0.1406] | 70% | +20.8 [+1.3, +54.3] | 79% | yes |
| 0.5 | 0 | 0.0028 [0.0004, 0.0121] | 3% | +1.8 [-2.7, +8.0] | 3% | yes |
| 0.5 | 0.1 | 0.0034 [0.0007, 0.0266] | 11% | +2.1 [-2.8, +10.0] | 6% | yes |
| 0.5 | 0.2 | 0.0075 [0.0009, 0.0584] | 34% | +3.6 [-1.7, +17.6] | 25% | yes |
| 0.5 | 0.3 | 0.0170 [0.0013, 0.0614] | 59% | +7.0 [-0.8, +20.3] | 45% | no |
| 1 | 0 | 0.0024 [0.0006, 0.0116] | 2% | +1.2 [-3.3, +6.7] | 1% | yes |
| 1 | 0.1 | 0.0030 [0.0006, 0.0331] | 11% | +1.7 [-2.2, +8.8] | 5% | yes |
| 1 | 0.2 | 0.0056 [0.0006, 0.0407] | 30% | +2.9 [-1.9, +11.4] | 13% | yes |
| 1 | 0.3 | 0.0084 [0.0009, 0.0430] | 39% | +3.8 [-2.0, +14.3] | 17% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0067**; boundary excess at the Parkin line **+5.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0043 [0.0006, 0.0209] | 33% | +2.7 [-5.4, +14.5] | 28% |
| 0.03 | 0.1 | 0.0057 [0.0005, 0.1417] | 46% | +3.4 [-8.6, +67.0] | 36% |
| 0.03 | 0.2 | 0.0127 [0.0011, 0.1752] | 62% | +7.7 [-11.2, +91.0] | 55% |
| 0.03 | 0.3 | 0.0412 [0.0010, 0.1848] | 79% | +22.4 [-16.9, +96.3] | 69% |
| 0.1 | 0 | 0.0036 [0.0004, 0.0206] | 25% | +1.6 [-4.7, +11.1] | 17% |
| 0.1 | 0.1 | 0.0054 [0.0006, 0.1082] | 41% | +2.5 [-5.0, +49.6] | 30% |
| 0.1 | 0.2 | 0.0090 [0.0006, 0.1085] | 57% | +4.9 [-6.6, +55.7] | 47% |
| 0.1 | 0.3 | 0.0260 [0.0011, 0.1359] | 77% | +14.0 [-8.1, +69.5] | 64% |
| 0.5 | 0 | 0.0019 [0.0002, 0.0118] | 9% | +0.4 [-3.5, +7.4] | 4% |
| 0.5 | 0.1 | 0.0027 [0.0003, 0.0261] | 20% | +0.6 [-4.3, +13.2] | 12% |
| 0.5 | 0.2 | 0.0061 [0.0003, 0.0567] | 47% | +1.6 [-4.8, +25.9] | 32% |
| 0.5 | 0.3 | 0.0152 [0.0007, 0.0599] | 66% | +6.0 [-4.4, +29.6] | 52% |
| 1 | 0 | 0.0017 [0.0002, 0.0112] | 7% | +0.0 [-3.9, +5.9] | 3% |
| 1 | 0.1 | 0.0024 [0.0003, 0.0317] | 19% | +0.5 [-3.8, +13.9] | 9% |
| 1 | 0.2 | 0.0050 [0.0003, 0.0393] | 40% | +1.6 [-4.5, +17.4] | 25% |
| 1 | 0.3 | 0.0068 [0.0004, 0.0413] | 50% | +2.3 [-4.1, +19.9] | 32% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0123 against a median of 0.0024 (95 percent range 0.0006 to 0.0116); 2% of runs reach it. Boundary excess at phase lines: observed +7.7 against +1.2 (-3.3 to +6.7); 1% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.75** (prior 0.75).
- P(local innovation, share > 0) = **0.68** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.73 | 0.60 |
| 1% | phase F_ST and boundary excess only | 0.79 | 0.75 |
| 2% | all six summaries | 0.78 | 0.62 |
| 2% | phase F_ST and boundary excess only | 0.79 | 0.74 |
| 5% | all six summaries | 0.75 | 0.68 |
| 5% | phase F_ST and boundary excess only | 0.79 | 0.79 |
| 10% | all six summaries | 0.73 | 0.66 |
| 10% | phase F_ST and boundary excess only | 0.81 | 0.79 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.250 | 0.862 |
| 0.5 | 0.164 | 0.757 |
| 0.1 | 0.049 | 0.384 |
| 0.03 | 0.017 | 0.158 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.25 | | 0 | 0.32 |
| 0.1 | 0.29 | | 0.1 | 0.33 |
| 0.5 | 0.21 | | 0.2 | 0.19 |
| 1 | 0.25 | | 0.3 | 0.17 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.79 | 0.64 |
| 0.1 | 0.2 | 0.82 | 0.62 |
| 0.03 | 0 | 0.82 | 0.70 |
| 0.03 | 0.2 | 0.94 | 0.71 |
| 1 (no boundary) | 0 | 0.66 | 0.56 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
