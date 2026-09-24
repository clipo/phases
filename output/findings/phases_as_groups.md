# Do the published phases behave as groups?

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the published phases (Kent 8, Parkin 11, Walls 9).

Observed: between-phase cultural F_ST **0.0110**; boundary excess at phase lines **+10.0**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0093 [0.0017, 0.0337] | 40% | +7.6 [-0.3, +20.5] | 34% | yes |
| 0.03 | 0.1 | 0.0109 [0.0023, 0.1326] | 50% | +10.3 [-1.0, +57.4] | 52% | yes |
| 0.03 | 0.2 | 0.0265 [0.0021, 0.1916] | 67% | +19.0 [-0.1, +77.1] | 65% | yes |
| 0.03 | 0.3 | 0.0607 [0.0026, 0.1962] | 80% | +34.8 [+2.0, +84.3] | 79% | no |
| 0.1 | 0 | 0.0068 [0.0014, 0.0233] | 28% | +4.4 [-5.4, +15.1] | 15% | yes |
| 0.1 | 0.1 | 0.0074 [0.0014, 0.1174] | 37% | +6.3 [-2.8, +38.9] | 27% | yes |
| 0.1 | 0.2 | 0.0158 [0.0021, 0.1229] | 63% | +10.8 [-1.8, +44.4] | 53% | yes |
| 0.1 | 0.3 | 0.0357 [0.0032, 0.1492] | 78% | +20.6 [+0.4, +54.6] | 72% | no |
| 0.5 | 0 | 0.0042 [0.0009, 0.0212] | 10% | +2.8 [-4.5, +10.8] | 4% | yes |
| 0.5 | 0.1 | 0.0050 [0.0010, 0.0350] | 22% | +2.8 [-6.0, +18.8] | 10% | yes |
| 0.5 | 0.2 | 0.0096 [0.0012, 0.0642] | 45% | +5.7 [-4.4, +26.5] | 32% | yes |
| 0.5 | 0.3 | 0.0172 [0.0026, 0.0707] | 65% | +8.8 [-2.8, +28.7] | 47% | no |
| 1 | 0 | 0.0035 [0.0007, 0.0145] | 4% | +1.5 [-6.0, +7.8] | 0.3% | yes |
| 1 | 0.1 | 0.0044 [0.0009, 0.0282] | 16% | +1.9 [-6.8, +12.7] | 5% | yes |
| 1 | 0.2 | 0.0077 [0.0008, 0.0420] | 39% | +4.0 [-5.4, +21.2] | 23% | yes |
| 1 | 0.3 | 0.0116 [0.0015, 0.0500] | 52% | +6.3 [-4.9, +23.7] | 33% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0062**; boundary excess at the Parkin line **+13.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0055 [0.0006, 0.0300] | 45% | +5.2 [-3.7, +21.1] | 11% |
| 0.03 | 0.1 | 0.0074 [0.0009, 0.1152] | 55% | +6.9 [-2.8, +59.6] | 19% |
| 0.03 | 0.2 | 0.0168 [0.0009, 0.1791] | 70% | +10.4 [-3.3, +85.4] | 41% |
| 0.03 | 0.3 | 0.0379 [0.0011, 0.1834] | 82% | +19.0 [-2.0, +86.3] | 57% |
| 0.1 | 0 | 0.0046 [0.0004, 0.0192] | 38% | +2.8 [-6.5, +15.9] | 5% |
| 0.1 | 0.1 | 0.0051 [0.0005, 0.1154] | 42% | +3.6 [-5.5, +51.7] | 13% |
| 0.1 | 0.2 | 0.0108 [0.0010, 0.1155] | 65% | +5.8 [-4.2, +56.2] | 24% |
| 0.1 | 0.3 | 0.0265 [0.0013, 0.1448] | 79% | +11.9 [-3.5, +66.0] | 47% |
| 0.5 | 0 | 0.0026 [0.0002, 0.0194] | 22% | +1.7 [-5.8, +12.3] | 2% |
| 0.5 | 0.1 | 0.0033 [0.0002, 0.0286] | 32% | +1.1 [-6.4, +15.1] | 4% |
| 0.5 | 0.2 | 0.0069 [0.0004, 0.0563] | 53% | +3.4 [-6.9, +25.8] | 14% |
| 0.5 | 0.3 | 0.0134 [0.0008, 0.0669] | 69% | +5.6 [-6.3, +29.8] | 24% |
| 1 | 0 | 0.0019 [0.0002, 0.0110] | 12% | +0.5 [-6.5, +7.5] | 0% |
| 1 | 0.1 | 0.0028 [0.0002, 0.0247] | 24% | +1.0 [-5.9, +11.1] | 2% |
| 1 | 0.2 | 0.0054 [0.0003, 0.0392] | 47% | +2.9 [-6.8, +21.4] | 9% |
| 1 | 0.3 | 0.0089 [0.0006, 0.0472] | 59% | +4.3 [-7.7, +25.2] | 16% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0110 against a median of 0.0035 (95 percent range 0.0007 to 0.0145); 4% of runs reach it. Boundary excess at phase lines: observed +10.0 against +1.5 (-6.0 to +7.8); 0.3% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.79** (prior 0.75).
- P(local innovation, share > 0) = **0.65** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.83 | 0.73 |
| 1% | phase F_ST and boundary excess only | 0.81 | 0.69 |
| 2% | all six summaries | 0.82 | 0.70 |
| 2% | phase F_ST and boundary excess only | 0.80 | 0.70 |
| 5% | all six summaries | 0.79 | 0.65 |
| 5% | phase F_ST and boundary excess only | 0.82 | 0.69 |
| 10% | all six summaries | 0.75 | 0.67 |
| 10% | phase F_ST and boundary excess only | 0.83 | 0.73 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.239 | 0.991 |
| 0.5 | 0.178 | 0.982 |
| 0.1 | 0.092 | 0.917 |
| 0.03 | 0.053 | 0.769 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.25 | | 0 | 0.35 |
| 0.1 | 0.28 | | 0.1 | 0.26 |
| 0.5 | 0.25 | | 0.2 | 0.21 |
| 1 | 0.21 | | 0.3 | 0.18 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.81 | 0.62 |
| 0.1 | 0.2 | 0.84 | 0.64 |
| 0.03 | 0 | 0.81 | 0.66 |
| 0.03 | 0.2 | 0.95 | 0.75 |
| 1 (no boundary) | 0 | 0.66 | 0.53 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
