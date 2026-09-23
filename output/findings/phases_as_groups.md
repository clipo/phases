# Do the published phases behave as groups?

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.002, mixing 0.005), interaction length 24 km along rivers. Groups are the published phases (Kent 8, Parkin 11, Walls 9).

Observed: between-phase cultural F_ST **0.0110**; boundary excess at phase lines **+10.7**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0062 [0.0012, 0.0223] | 20% | +2.8 [-7.6, +11.6] | 5% | yes |
| 0.03 | 0.1 | 0.0068 [0.0013, 0.1675] | 29% | +4.0 [-6.8, +58.3] | 18% | yes |
| 0.03 | 0.2 | 0.0187 [0.0016, 0.2238] | 59% | +12.0 [-5.7, +85.2] | 53% | yes |
| 0.03 | 0.3 | 0.0841 [0.0028, 0.2252] | 79% | +42.6 [-1.6, +98.8] | 76% | no |
| 0.1 | 0 | 0.0057 [0.0012, 0.0222] | 20% | +2.3 [-6.4, +11.6] | 3% | yes |
| 0.1 | 0.1 | 0.0081 [0.0010, 0.1613] | 37% | +3.7 [-8.3, +48.9] | 20% | yes |
| 0.1 | 0.2 | 0.0177 [0.0015, 0.1858] | 60% | +10.8 [-6.2, +74.2] | 50% | yes |
| 0.1 | 0.3 | 0.0730 [0.0034, 0.2011] | 82% | +33.8 [-4.6, +80.8] | 78% | no |
| 0.5 | 0 | 0.0057 [0.0010, 0.0212] | 16% | +1.9 [-9.2, +10.0] | 2% | yes |
| 0.5 | 0.1 | 0.0065 [0.0010, 0.1188] | 28% | +2.2 [-10.1, +39.1] | 15% | yes |
| 0.5 | 0.2 | 0.0110 [0.0013, 0.1571] | 50% | +5.8 [-6.3, +47.7] | 36% | yes |
| 0.5 | 0.3 | 0.0541 [0.0015, 0.1631] | 75% | +17.1 [-7.4, +58.1] | 60% | no |
| 1 | 0 | 0.0049 [0.0012, 0.0164] | 11% | +0.8 [-9.2, +8.0] | 0% | yes |
| 1 | 0.1 | 0.0067 [0.0012, 0.1065] | 29% | +1.9 [-8.1, +30.5] | 13% | yes |
| 1 | 0.2 | 0.0155 [0.0016, 0.1422] | 56% | +5.8 [-7.9, +44.7] | 36% | yes |
| 1 | 0.3 | 0.0435 [0.0019, 0.1564] | 79% | +17.2 [-5.1, +54.1] | 61% | no |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0110 against a median of 0.0049 (95 percent range 0.0012 to 0.0164); 11% of runs reach it. Boundary excess at phase lines: observed +10.7 against +0.8 (-9.2 to +8.0); 0% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.68** (prior 0.75).
- P(local innovation, share > 0) = **0.70** (prior 0.75).

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.17 | | 0 | 0.30 |
| 0.1 | 0.20 | | 0.1 | 0.27 |
| 0.5 | 0.30 | | 0.2 | 0.28 |
| 1 | 0.32 | | 0.3 | 0.14 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.75 | 0.69 |
| 0.1 | 0.2 | 0.79 | 0.72 |
| 0.03 | 0 | 0.78 | 0.72 |
| 0.03 | 0.2 | 0.79 | 0.70 |
| 1 (no boundary) | 0 | 0.75 | 0.69 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
