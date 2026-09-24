# Do the published phases behave as groups? (placebo 12)

PLACEBO 12: the groups are NOT the phases but a division with the phases' sizes around random centers (adjusted Rand with the phases 0.399); 'phase' below means this division.

Produced by `analyses/84_phases_as_groups.py`, 300 runs per cell, calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 km along rivers. Groups are the placebo division (Kent 8, Parkin 11, Walls 9).

Observed: between-phase cultural F_ST **0.0098**; boundary excess at phase lines **+5.3**; F_ST at the 2 spatial clusters 0.0062.

## 1. Every cell

Leak multiplies copying across phase lines (1 = no copying boundary). Local innovation is the share of classes reordered in each phase's source of new variants (0 = one regional pool).

| leak | local innovation | phase F_ST median [95%] | share reaching observed | boundary excess median [95%] | share reaching observed | diversity matched |
|---|---|---|---|---|---|---|
| 0.03 | 0 | 0.0069 [0.0012, 0.0264] | 37% | +3.2 [-1.8, +12.5] | 31% | no |
| 0.03 | 0.1 | 0.0083 [0.0018, 0.0815] | 44% | +4.7 [-1.7, +36.9] | 44% | yes |
| 0.03 | 0.2 | 0.0203 [0.0026, 0.1025] | 68% | +10.9 [-1.1, +51.1] | 67% | yes |
| 0.03 | 0.3 | 0.0360 [0.0033, 0.1002] | 80% | +21.3 [-0.4, +54.9] | 80% | no |
| 0.1 | 0 | 0.0041 [0.0009, 0.0179] | 13% | +1.8 [-2.3, +7.0] | 10% | yes |
| 0.1 | 0.1 | 0.0055 [0.0009, 0.0408] | 25% | +2.3 [-1.9, +18.2] | 27% | yes |
| 0.1 | 0.2 | 0.0090 [0.0011, 0.0486] | 47% | +5.1 [-1.1, +24.2] | 49% | yes |
| 0.1 | 0.3 | 0.0178 [0.0017, 0.0505] | 67% | +10.1 [-0.7, +26.1] | 68% | no |
| 0.5 | 0 | 0.0024 [0.0004, 0.0106] | 4% | +0.6 [-1.9, +4.7] | 1% | yes |
| 0.5 | 0.1 | 0.0027 [0.0006, 0.0143] | 7% | +1.0 [-1.7, +6.8] | 5% | yes |
| 0.5 | 0.2 | 0.0042 [0.0009, 0.0195] | 17% | +1.9 [-1.6, +8.7] | 12% | yes |
| 0.5 | 0.3 | 0.0058 [0.0011, 0.0213] | 29% | +2.6 [-1.1, +10.0] | 23% | no |
| 1 | 0 | 0.0023 [0.0005, 0.0101] | 3% | +0.7 [-1.6, +4.2] | 1% | yes |
| 1 | 0.1 | 0.0026 [0.0005, 0.0114] | 3% | +0.6 [-1.9, +5.2] | 2% | yes |
| 1 | 0.2 | 0.0038 [0.0006, 0.0186] | 12% | +1.2 [-1.6, +7.1] | 6% | yes |
| 1 | 0.3 | 0.0046 [0.0010, 0.0167] | 17% | +2.0 [-1.6, +7.8] | 12% | no |

## The Parkin phase against the rest of the basin

Observed: Parkin-versus-rest cultural F_ST **0.0062**; boundary excess at the Parkin line **+13.7**. The shares are of 300 runs per cell reaching the observed value.

| copying factor | local innovation | Parkin F_ST median [95%] | share reaching | Parkin-line boundary excess median [95%] | share reaching |
|---|---|---|---|---|---|
| 0.03 | 0 | 0.0030 [0.0004, 0.0186] | 25% | +0.4 [-9.5, +9.2] | 0.7% |
| 0.03 | 0.1 | 0.0040 [0.0003, 0.0427] | 36% | +0.5 [-11.0, +13.4] | 2% |
| 0.03 | 0.2 | 0.0065 [0.0004, 0.0531] | 52% | +0.6 [-11.6, +13.3] | 2% |
| 0.03 | 0.3 | 0.0112 [0.0005, 0.0513] | 65% | +1.1 [-11.4, +14.7] | 3% |
| 0.1 | 0 | 0.0029 [0.0002, 0.0148] | 22% | -0.0 [-7.4, +7.8] | 0% |
| 0.1 | 0.1 | 0.0034 [0.0003, 0.0312] | 29% | +0.4 [-8.0, +11.6] | 1% |
| 0.1 | 0.2 | 0.0048 [0.0004, 0.0404] | 43% | +0.4 [-9.9, +13.9] | 3% |
| 0.1 | 0.3 | 0.0085 [0.0007, 0.0410] | 55% | +0.8 [-10.8, +15.8] | 5% |
| 0.5 | 0 | 0.0019 [0.0002, 0.0126] | 17% | +0.5 [-6.5, +8.8] | 0.3% |
| 0.5 | 0.1 | 0.0025 [0.0003, 0.0176] | 17% | +1.0 [-6.3, +11.5] | 1% |
| 0.5 | 0.2 | 0.0037 [0.0004, 0.0247] | 34% | +1.4 [-7.5, +12.7] | 2% |
| 0.5 | 0.3 | 0.0057 [0.0005, 0.0274] | 46% | +1.5 [-8.1, +15.9] | 5% |
| 1 | 0 | 0.0019 [0.0002, 0.0110] | 12% | +0.5 [-6.5, +7.5] | 0% |
| 1 | 0.1 | 0.0028 [0.0002, 0.0151] | 18% | +0.8 [-5.6, +9.8] | 0.7% |
| 1 | 0.2 | 0.0041 [0.0003, 0.0251] | 33% | +1.7 [-6.3, +16.4] | 3% |
| 1 | 0.3 | 0.0052 [0.0006, 0.0224] | 42% | +2.2 [-6.6, +14.1] | 3% |

## 2. The phases under neutral copying alone (leak 1, regional pool)

Between-phase F_ST: observed 0.0098 against a median of 0.0023 (95 percent range 0.0005 to 0.0101); 3% of runs reach it. Boundary excess at phase lines: observed +5.3 against +0.7 (-1.6 to +4.2); 1% of runs reach it.

## 3. Posterior over the grid (rejection ABC, uniform prior over cells)

Accepted the closest 5% of 4800 runs on standardized phase_fst, spatial_fst, be_phase, hs, rich, ht.

- P(some copying boundary at phase lines, leak < 1) = **0.78** (prior 0.75).
- P(local innovation, share > 0) = **0.68** (prior 0.75).

Sensitivity to the approximation's settings (prior 0.75 for both):

| acceptance | summaries | P(copying boundary) | P(local innovation) |
|---|---|---|---|
| 1% | all six summaries | 0.77 | 0.65 |
| 1% | phase F_ST and boundary excess only | 0.77 | 0.79 |
| 2% | all six summaries | 0.78 | 0.70 |
| 2% | phase F_ST and boundary excess only | 0.78 | 0.78 |
| 5% | all six summaries | 0.78 | 0.68 |
| 5% | phase F_ST and boundary excess only | 0.81 | 0.81 |
| 10% | all six summaries | 0.78 | 0.68 |
| 10% | phase F_ST and boundary excess only | 0.81 | 0.81 |

What the copying factor does to cross-phase copying (share of each site's between-site copying that crosses a phase line, mean and maximum over sites):

| copying factor | mean share | maximum share |
|---|---|---|
| 1 | 0.442 | 0.991 |
| 0.5 | 0.317 | 0.983 |
| 0.1 | 0.125 | 0.918 |
| 0.03 | 0.060 | 0.772 |

| leak | posterior | | local innovation | posterior |
|---|---|---|---|---|
| 0.03 | 0.19 | | 0 | 0.32 |
| 0.1 | 0.39 | | 0.1 | 0.30 |
| 0.5 | 0.20 | | 0.2 | 0.23 |
| 1 | 0.22 | | 0.3 | 0.15 |

## 4. Recovery (rule 20c)

Pseudo-observations drawn from cells WITH a copying boundary, run through the same posterior (each excluded from its own reference table). If the design can see a boundary, P(leak < 1) should rise well above the prior of 0.75.

| true leak | true local innovation | P(leak < 1), median over 40 | 10th percentile |
|---|---|---|---|
| 0.1 | 0 | 0.74 | 0.60 |
| 0.1 | 0.2 | 0.75 | 0.57 |
| 0.03 | 0 | 0.86 | 0.77 |
| 0.03 | 0.2 | 0.93 | 0.76 |
| 1 (no boundary) | 0 | 0.67 | 0.55 |

## Reading

Compare the observed posterior P(leak < 1) with the recovery rows. A value near the prior, or near the no-boundary row, says the data do not ask for a copying boundary at the phase lines once local innovation is available; a value near the boundary rows says they do.
