# Do the published phases behave as groups? (placebo 13)

PLACEBO 13: the groups are NOT the phases but a division with the phases' sizes around random centers (adjusted Rand with the phases 0.346); 'phase' below means this division.

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the placebo division (Kent 11, Parkin 12, Walls 5).

Observed: between-phase cultural F_ST **0.0234**; boundary excess at phase lines **+13.5**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0059 [0.0014, 0.0306] | 5% | +3.4 [-2.8, +14.1] | 3% | no |
| 0.03 | 0.1 | 0.0085 [0.0012, 0.1401] | 17% | +4.9 [-2.4, +43.3] | 17% | yes |
| 0.03 | 0.2 | 0.0189 [0.0020, 0.1873] | 45% | +12.2 [-1.9, +62.6] | 48% | yes |
| 0.03 | 0.3 | 0.0472 [0.0026, 0.1786] | 61% | +26.0 [+1.1, +70.2] | 66% | yes |
| 0.1 | 0 | 0.0046 [0.0009, 0.0225] | 2% | +1.9 [-2.4, +9.2] | 0.3% | yes |
| 0.1 | 0.1 | 0.0066 [0.0011, 0.1002] | 16% | +3.4 [-2.0, +29.1] | 13% | yes |
| 0.1 | 0.2 | 0.0103 [0.0009, 0.1251] | 32% | +5.9 [-1.8, +36.6] | 33% | yes |
| 0.1 | 0.3 | 0.0304 [0.0024, 0.1340] | 55% | +14.0 [+0.2, +43.8] | 52% | yes |
| 0.5 | 0 | 0.0028 [0.0004, 0.0132] | 0% | +1.3 [-2.6, +7.6] | 0% | yes |
| 0.5 | 0.1 | 0.0037 [0.0006, 0.0272] | 4% | +1.7 [-1.9, +11.8] | 1% | yes |
| 0.5 | 0.2 | 0.0064 [0.0009, 0.0501] | 18% | +3.7 [-1.5, +18.6] | 8% | yes |
| 0.5 | 0.3 | 0.0134 [0.0013, 0.0574] | 30% | +7.0 [-1.2, +22.2] | 19% | no |
| 1 | 0 | 0.0023 [0.0005, 0.0113] | 0% | +1.1 [-2.1, +7.0] | 0% | yes |
| 1 | 0.1 | 0.0031 [0.0005, 0.0210] | 2% | +1.7 [-2.2, +10.3] | 0% | yes |
| 1 | 0.2 | 0.0048 [0.0007, 0.0330] | 8% | +2.7 [-2.3, +14.5] | 3% | yes |
| 1 | 0.3 | 0.0072 [0.0010, 0.0377] | 15% | +4.3 [-2.0, +15.3] | 5% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0067**; boundary excess at the Parkin line **+5.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0039 [0.0006, 0.0298] | 29% | -0.1 [-5.5, +8.1] | 6% |
| 0.03 | 0.1 | 0.0053 [0.0007, 0.0913] | 42% | +0.8 [-6.3, +24.2] | 21% |
| 0.03 | 0.2 | 0.0124 [0.0008, 0.1457] | 64% | +4.2 [-5.5, +33.1] | 42% |
| 0.03 | 0.3 | 0.0311 [0.0010, 0.1368] | 77% | +9.1 [-2.8, +37.1] | 60% |
| 0.1 | 0 | 0.0036 [0.0003, 0.0193] | 25% | +0.4 [-5.1, +6.6] | 4% |
| 0.1 | 0.1 | 0.0050 [0.0006, 0.0854] | 41% | +1.3 [-4.2, +24.1] | 20% |
| 0.1 | 0.2 | 0.0083 [0.0006, 0.1073] | 56% | +3.7 [-3.5, +31.3] | 37% |
| 0.1 | 0.3 | 0.0254 [0.0010, 0.1136] | 71% | +8.5 [-3.1, +34.3] | 59% |
| 0.5 | 0 | 0.0021 [0.0002, 0.0124] | 13% | +0.1 [-3.5, +5.7] | 3% |
| 0.5 | 0.1 | 0.0030 [0.0003, 0.0264] | 18% | +0.4 [-4.1, +10.5] | 8% |
| 0.5 | 0.2 | 0.0052 [0.0004, 0.0498] | 45% | +2.1 [-3.1, +21.1] | 27% |
| 0.5 | 0.3 | 0.0116 [0.0005, 0.0557] | 65% | +5.4 [-2.9, +23.0] | 48% |
| 1 | 0 | 0.0017 [0.0002, 0.0112] | 7% | +0.0 [-3.9, +5.9] | 3% |
| 1 | 0.1 | 0.0024 [0.0003, 0.0222] | 17% | +0.5 [-3.8, +9.9] | 9% |
| 1 | 0.2 | 0.0045 [0.0004, 0.0325] | 38% | +1.6 [-3.2, +14.8] | 20% |
| 1 | 0.3 | 0.0067 [0.0005, 0.0387] | 50% | +3.3 [-2.5, +16.1] | 33% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0234 against a median of 0.0023 (95 percent range 0.0005 to 0.0113); 0% of runs reach it. Boundary excess at phase lines: observed +13.5 against +1.1 (-2.1 to +7.0); 0% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.85** (prior 0.75).
- P(local innovation, share > 0) = **0.80** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.94 | 0.90 |
| 1% | phase F_ST and boundary excess only | 0.77 | 1.00 |
| 2% | all six summaries | 0.95 | 0.90 |
| 2% | phase F_ST and boundary excess only | 0.79 | 0.97 |
| 5% | all six summaries | 0.85 | 0.80 |
| 5% | phase F_ST and boundary excess only | 0.81 | 0.95 |
| 10% | all six summaries | 0.83 | 0.78 |
| 10% | phase F_ST and boundary excess only | 0.78 | 0.93 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.354 | 0.979 |
| 0.5 | 0.252 | 0.959 |
| 0.1 | 0.094 | 0.825 |
| 0.03 | 0.042 | 0.585 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.38 | | 0 | 0.20 |
| 0.1 | 0.30 | | 0.1 | 0.31 |
| 0.5 | 0.18 | | 0.2 | 0.24 |
| 1 | 0.15 | | 0.3 | 0.25 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.75 | 0.60 |
| 0.1 | 0.2 | 0.80 | 0.64 |
| 0.03 | 0 | 0.82 | 0.65 |
| 0.03 | 0.2 | 0.93 | 0.71 |
| 1 (no boundary) | 0 | 0.66 | 0.54 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
