# Do the published phases behave as groups? (placebo 12)

PLACEBO 12: the groups are NOT the phases but a division with the phases' sizes around random centers (adjusted Rand with the phases 0.171); 'phase' below means this division.

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the placebo division (Kent 11, Parkin 12, Walls 5).

Observed: between-phase cultural F_ST **0.0112**; boundary excess at phase lines **+1.3**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0075 [0.0016, 0.0262] | 26% | +4.1 [-1.1, +14.9] | 81% | yes |
| 0.03 | 0.1 | 0.0092 [0.0016, 0.0857] | 39% | +5.3 [-1.4, +48.0] | 84% | yes |
| 0.03 | 0.2 | 0.0219 [0.0021, 0.1301] | 64% | +13.3 [-0.8, +61.2] | 92% | yes |
| 0.03 | 0.3 | 0.0419 [0.0034, 0.1349] | 76% | +26.9 [+0.6, +71.1] | 96% | no |
| 0.1 | 0 | 0.0042 [0.0009, 0.0145] | 7% | +1.9 [-1.7, +7.7] | 59% | yes |
| 0.1 | 0.1 | 0.0045 [0.0010, 0.0458] | 24% | +2.7 [-1.9, +20.4] | 68% | yes |
| 0.1 | 0.2 | 0.0081 [0.0012, 0.0674] | 45% | +5.3 [-2.2, +27.1] | 78% | yes |
| 0.1 | 0.3 | 0.0180 [0.0015, 0.0645] | 62% | +10.9 [-0.7, +34.5] | 89% | yes |
| 0.5 | 0 | 0.0025 [0.0006, 0.0105] | 2% | +0.6 [-2.8, +4.2] | 35% | yes |
| 0.5 | 0.1 | 0.0028 [0.0005, 0.0192] | 8% | +0.6 [-2.5, +6.9] | 39% | yes |
| 0.5 | 0.2 | 0.0043 [0.0010, 0.0350] | 22% | +1.4 [-3.1, +9.6] | 52% | yes |
| 0.5 | 0.3 | 0.0080 [0.0009, 0.0371] | 37% | +2.6 [-2.0, +11.8] | 64% | no |
| 1 | 0 | 0.0022 [0.0005, 0.0100] | 1% | +0.1 [-3.2, +4.1] | 28% | yes |
| 1 | 0.1 | 0.0028 [0.0005, 0.0168] | 7% | +0.8 [-2.3, +5.5] | 39% | yes |
| 1 | 0.2 | 0.0044 [0.0005, 0.0230] | 21% | +1.1 [-2.9, +7.8] | 45% | yes |
| 1 | 0.3 | 0.0069 [0.0009, 0.0293] | 31% | +1.6 [-2.2, +8.6] | 55% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0067**; boundary excess at the Parkin line **+5.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0027 [0.0003, 0.0140] | 19% | -0.9 [-7.2, +4.9] | 1% |
| 0.03 | 0.1 | 0.0036 [0.0004, 0.0337] | 26% | -1.4 [-8.4, +4.3] | 0% |
| 0.03 | 0.2 | 0.0088 [0.0006, 0.0654] | 56% | -2.4 [-10.7, +5.6] | 3% |
| 0.03 | 0.3 | 0.0146 [0.0005, 0.0657] | 68% | -2.6 [-11.1, +6.8] | 3% |
| 0.1 | 0 | 0.0020 [0.0002, 0.0134] | 10% | -0.7 [-5.7, +4.6] | 0.7% |
| 0.1 | 0.1 | 0.0028 [0.0003, 0.0348] | 26% | -1.0 [-8.6, +5.8] | 3% |
| 0.1 | 0.2 | 0.0052 [0.0003, 0.0450] | 43% | -1.2 [-10.4, +6.3] | 3% |
| 0.1 | 0.3 | 0.0086 [0.0006, 0.0494] | 56% | -2.0 [-12.1, +5.3] | 2% |
| 0.5 | 0 | 0.0019 [0.0003, 0.0122] | 9% | -0.0 [-3.6, +5.5] | 3% |
| 0.5 | 0.1 | 0.0025 [0.0002, 0.0164] | 17% | -0.0 [-4.5, +6.0] | 3% |
| 0.5 | 0.2 | 0.0041 [0.0003, 0.0333] | 31% | +0.2 [-6.2, +8.1] | 8% |
| 0.5 | 0.3 | 0.0070 [0.0005, 0.0338] | 52% | -0.0 [-8.1, +8.3] | 11% |
| 1 | 0 | 0.0017 [0.0002, 0.0112] | 7% | +0.0 [-3.9, +5.9] | 3% |
| 1 | 0.1 | 0.0023 [0.0002, 0.0169] | 19% | +0.1 [-4.2, +8.3] | 7% |
| 1 | 0.2 | 0.0040 [0.0002, 0.0267] | 32% | +0.5 [-5.0, +9.2] | 10% |
| 1 | 0.3 | 0.0067 [0.0006, 0.0313] | 50% | +1.1 [-5.0, +9.9] | 16% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0112 against a median of 0.0022 (95 percent range 0.0005 to 0.0100); 1% of runs reach it. Boundary excess at phase lines: observed +1.3 against +0.1 (-3.2 to +4.1); 28% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.70** (prior 0.75).
- P(local innovation, share > 0) = **0.66** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.62 | 0.65 |
| 1% | phase F_ST and boundary excess only | 0.67 | 0.69 |
| 2% | all six summaries | 0.65 | 0.66 |
| 2% | phase F_ST and boundary excess only | 0.66 | 0.79 |
| 5% | all six summaries | 0.70 | 0.66 |
| 5% | phase F_ST and boundary excess only | 0.65 | 0.78 |
| 10% | all six summaries | 0.71 | 0.64 |
| 10% | phase F_ST and boundary excess only | 0.65 | 0.75 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.343 | 0.889 |
| 0.5 | 0.236 | 0.801 |
| 0.1 | 0.077 | 0.445 |
| 0.03 | 0.027 | 0.194 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.14 | | 0 | 0.34 |
| 0.1 | 0.30 | | 0.1 | 0.30 |
| 0.5 | 0.27 | | 0.2 | 0.18 |
| 1 | 0.30 | | 0.3 | 0.18 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.75 | 0.62 |
| 0.1 | 0.2 | 0.79 | 0.63 |
| 0.03 | 0 | 0.82 | 0.69 |
| 0.03 | 0.2 | 0.93 | 0.78 |
| 1 (no boundary) | 0 | 0.69 | 0.58 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
