# Do the published phases behave as groups? (placebo 13)

PLACEBO 13: the groups are NOT the phases but a division with the phases' sizes around random centers (adjusted Rand with the phases 0.399); 'phase' below means this division.

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the placebo division (Kent 8, Parkin 11, Walls 9).

Observed: between-phase cultural F_ST **0.0179**; boundary excess at phase lines **+14.2**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0089 [0.0017, 0.0385] | 19% | +5.1 [-0.3, +14.5] | 3% | no |
| 0.03 | 0.1 | 0.0112 [0.0014, 0.1259] | 27% | +5.9 [-0.3, +34.2] | 18% | yes |
| 0.03 | 0.2 | 0.0240 [0.0025, 0.1894] | 59% | +12.9 [+0.0, +53.1] | 47% | yes |
| 0.03 | 0.3 | 0.0565 [0.0026, 0.1834] | 74% | +23.3 [+1.6, +60.3] | 67% | no |
| 0.1 | 0 | 0.0061 [0.0010, 0.0341] | 11% | +2.7 [-1.9, +11.3] | 0.7% | yes |
| 0.1 | 0.1 | 0.0076 [0.0016, 0.1032] | 24% | +3.6 [-1.6, +26.4] | 12% | yes |
| 0.1 | 0.2 | 0.0136 [0.0017, 0.1324] | 42% | +6.1 [-1.5, +33.1] | 28% | yes |
| 0.1 | 0.3 | 0.0367 [0.0024, 0.1548] | 69% | +15.0 [-0.0, +37.9] | 52% | no |
| 0.5 | 0 | 0.0039 [0.0008, 0.0158] | 2% | +1.6 [-2.6, +8.3] | 0% | yes |
| 0.5 | 0.1 | 0.0044 [0.0008, 0.0287] | 9% | +1.8 [-2.3, +12.0] | 0.7% | yes |
| 0.5 | 0.2 | 0.0078 [0.0013, 0.0550] | 28% | +3.2 [-1.7, +18.4] | 7% | yes |
| 0.5 | 0.3 | 0.0161 [0.0016, 0.0611] | 45% | +6.4 [-0.9, +19.0] | 15% | no |
| 1 | 0 | 0.0029 [0.0007, 0.0127] | 0.7% | +1.3 [-2.9, +6.4] | 0% | yes |
| 1 | 0.1 | 0.0038 [0.0007, 0.0233] | 6% | +1.4 [-1.7, +10.7] | 0% | yes |
| 1 | 0.2 | 0.0068 [0.0010, 0.0420] | 17% | +3.0 [-1.6, +13.7] | 2% | yes |
| 1 | 0.3 | 0.0099 [0.0011, 0.0402] | 23% | +3.4 [-1.1, +13.8] | 2% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0062**; boundary excess at the Parkin line **+13.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0055 [0.0006, 0.0337] | 46% | +3.8 [-5.5, +18.2] | 9% |
| 0.03 | 0.1 | 0.0080 [0.0007, 0.1075] | 57% | +4.4 [-4.8, +42.2] | 18% |
| 0.03 | 0.2 | 0.0170 [0.0010, 0.1848] | 72% | +11.5 [-3.9, +62.4] | 44% |
| 0.03 | 0.3 | 0.0372 [0.0011, 0.1714] | 84% | +16.6 [-2.9, +61.4] | 55% |
| 0.1 | 0 | 0.0046 [0.0004, 0.0325] | 38% | +2.3 [-6.1, +17.5] | 6% |
| 0.1 | 0.1 | 0.0058 [0.0005, 0.0984] | 47% | +3.3 [-5.3, +32.4] | 15% |
| 0.1 | 0.2 | 0.0109 [0.0006, 0.1276] | 62% | +5.8 [-4.8, +42.5] | 30% |
| 0.1 | 0.3 | 0.0270 [0.0016, 0.1541] | 83% | +14.6 [-2.9, +48.3] | 52% |
| 0.5 | 0 | 0.0029 [0.0003, 0.0139] | 21% | +1.3 [-5.8, +9.8] | 0.7% |
| 0.5 | 0.1 | 0.0032 [0.0003, 0.0283] | 31% | +1.8 [-6.0, +15.9] | 5% |
| 0.5 | 0.2 | 0.0065 [0.0005, 0.0539] | 52% | +3.9 [-5.2, +25.0] | 15% |
| 0.5 | 0.3 | 0.0137 [0.0007, 0.0616] | 70% | +7.0 [-3.9, +26.8] | 26% |
| 1 | 0 | 0.0019 [0.0002, 0.0110] | 12% | +0.5 [-6.5, +7.5] | 0% |
| 1 | 0.1 | 0.0028 [0.0002, 0.0229] | 20% | +1.0 [-5.9, +13.3] | 2% |
| 1 | 0.2 | 0.0053 [0.0003, 0.0421] | 46% | +2.8 [-4.8, +21.1] | 10% |
| 1 | 0.3 | 0.0086 [0.0005, 0.0387] | 62% | +4.8 [-4.2, +20.5] | 12% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0179 against a median of 0.0029 (95 percent range 0.0007 to 0.0127); 0.7% of runs reach it. Boundary excess at phase lines: observed +14.2 against +1.3 (-2.9 to +6.4); 0% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.92** (prior 0.75).
- P(local innovation, share > 0) = **0.70** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.96 | 0.71 |
| 1% | phase F_ST and boundary excess only | 0.90 | 0.83 |
| 2% | all six summaries | 0.95 | 0.71 |
| 2% | phase F_ST and boundary excess only | 0.89 | 0.84 |
| 5% | all six summaries | 0.92 | 0.70 |
| 5% | phase F_ST and boundary excess only | 0.85 | 0.83 |
| 10% | all six summaries | 0.88 | 0.68 |
| 10% | phase F_ST and boundary excess only | 0.84 | 0.86 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.384 | 0.999 |
| 0.5 | 0.281 | 0.999 |
| 0.1 | 0.126 | 0.995 |
| 0.03 | 0.076 | 0.983 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.39 | | 0 | 0.30 |
| 0.1 | 0.35 | | 0.1 | 0.22 |
| 0.5 | 0.18 | | 0.2 | 0.23 |
| 1 | 0.08 | | 0.3 | 0.26 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.76 | 0.64 |
| 0.1 | 0.2 | 0.77 | 0.59 |
| 0.03 | 0 | 0.83 | 0.72 |
| 0.03 | 0.2 | 0.92 | 0.73 |
| 1 (no boundary) | 0 | 0.66 | 0.52 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
