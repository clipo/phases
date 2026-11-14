# Regional drift at larger amplitudes and other occupation timings

Produced by `.venv-np254/bin/python analyses/per_type/drift_extension.py` at commit 3d67be6ff630 (with uncommitted changes to tracked files). NumPy 2.5.4; lock file `analyses/per_type/drift_extension_requirements.lock`, SHA-256 f51da1f4af069202; data hash 3b5e1cfabcad47f4. Run time: environment check 14 s, primary point 30 s, secondary point 17 s. Design: drift-extension spec of 2026-11-14, draft 3 (`docs/superpowers/specs/2026-11-14-drift-extension-timing-design.md`).

Generator C3 (a neutral copying model whose regional innovation target drifts), unchanged. Per operating point: 3 span-length limits x 3 windows x 7 couplings = 63 observation designs at 8 amplitudes F_W = 0, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9: 504 cells of 50 records. F_W is the regional F accumulated over the cell's own window. The 21 designs that share a window read the same simulated histories (1200 histories per point), so cells that share a window are dependent.

Every comparison is a neutral comparison in the main text's sense: the record's value of a statistic beside the 10th to 90th percentile band of the same statistic over the cell's simulated records (NumPy's linear quantile rule, inclusive endpoints; a band is eligible at 45 finite values of 50). These are descriptive bands, not posterior intervals. Statements are of marginal coverage: each named summary of the record lies inside its own band. No statement here says that any simulated record resembles the record in several summaries at once. Only a value above the band of a per-class F is read as more variation than the generator produces. Losses from the target are reported beside the cells and exclude none.

## The record

Plug-in values on observed shares, computed by the functions applied to the simulated records. They are not the posterior quantities the main text reports for the record. Sign conventions: a positive latitude correlation means the class's share is higher to the north; a positive D_dist means sites farther apart differ more in composition.

| statistic | Parkin_Punctated | Barton/Kent/MPI | Painted |
|---|---|---|---|
| F (site-weighted moment F) | 0.165 | 0.117 | 0.199 |
| latitude rank correlation (signed) | +0.540 | +0.057 | -0.775 |
| A_lat (absolute value) | 0.540 | 0.057 | 0.775 |

D_dist (rank correlation of total-variation distance with distance over 378 pairs): +0.360. Sites present for the seven rare classes: Fortune_Noded 15, Ranch_Incised 22, Walls_Engraved 14, Wallace_Incised 7, Rhodes_Incised 15, Vernon_Paul_Applique 10, Hull_Engraved 8.

## Primary point

Frozen scales: L_d 182, L_R 292, burn-in 5000. Windows and stored lengths: 8LR W = 2336, S = 2600; 2LR W = 584, S = 770; 32LR W = 9344, S = 9600. Span limits: Ld = 182, Ld/4 = 46, 1 = 1.

Environment check: the committed C3 cell at F_W 0.3 was regenerated through `per_class_table.one_record_c3m` and equals the committed findings file to three decimals for the three common classes (Parkin_Punctated 0.035 [0.011, 0.077]; Barton/Kent/MPI 0.032 [0.016, 0.111]; Painted 0.029 [0.014, 0.052]), with the same losses from the target.

### R1. Amplitude at the baseline design d0 (span limit L_d, window 8 L_R, coupling c0)

Cells: median [10th, 90th percentile] finite count status. Lost: records of 50 in which the class's target share is exactly zero at the last stored step.

| F_W | F Parkin_Punctated (record 0.165) | F Barton/Kent/MPI (record 0.117) | F Painted (record 0.199) | lost: Parkin_Punctated, Barton/Kent/MPI, Painted | D_dist (record +0.360) | A_lat Parkin_Punctated (record 0.540) | A_lat Barton/Kent/MPI (record 0.057) | A_lat Painted (record 0.775) |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.014 [0.009, 0.023] 50 above | 0.016 [0.007, 0.022] 50 above | 0.012 [0.006, 0.023] 50 above | 0, 0, 0 | +0.134 [+0.044, +0.292] 50 above | 0.254 [0.045, 0.537] 50 above | 0.285 [0.055, 0.544] 50 inside | 0.253 [0.052, 0.488] 50 above |
| 0.01 | 0.015 [0.009, 0.024] 50 above | 0.014 [0.006, 0.026] 50 above | 0.014 [0.006, 0.022] 50 above | 0, 0, 0 | +0.128 [+0.037, +0.263] 50 above | 0.229 [0.051, 0.488] 50 above | 0.203 [0.084, 0.412] 50 below | 0.193 [0.058, 0.393] 50 above |
| 0.03 | 0.019 [0.006, 0.032] 50 above | 0.016 [0.007, 0.038] 50 above | 0.016 [0.009, 0.033] 50 above | 0, 0, 0 | +0.133 [+0.026, +0.252] 50 above | 0.183 [0.022, 0.408] 50 above | 0.230 [0.045, 0.525] 50 inside | 0.260 [0.063, 0.487] 50 above |
| 0.1 | 0.021 [0.011, 0.036] 50 above | 0.022 [0.011, 0.043] 50 above | 0.024 [0.009, 0.040] 50 above | 0, 0, 12 | +0.107 [-0.012, +0.236] 50 above | 0.210 [0.049, 0.456] 50 above | 0.172 [0.027, 0.471] 50 inside | 0.265 [0.038, 0.478] 50 above |
| 0.3 | 0.036 [0.017, 0.108] 50 above | 0.037 [0.015, 0.116] 50 above | 0.030 [0.013, 0.114] 50 above | 1, 11, 30 | +0.073 [-0.018, +0.162] 50 above | 0.193 [0.024, 0.417] 50 above | 0.195 [0.042, 0.454] 50 inside | 0.159 [0.031, 0.398] 50 above |
| 0.5 | 0.059 [0.019, 0.188] 50 inside | 0.062 [0.026, 0.171] 50 inside | 0.040 [0.019, 0.082] 50 above | 8, 13, 38 | +0.043 [-0.032, +0.162] 50 above | 0.185 [0.026, 0.421] 50 above | 0.184 [0.030, 0.457] 50 inside | 0.193 [0.035, 0.451] 50 above |
| 0.7 | 0.110 [0.031, 0.211] 50 inside | 0.089 [0.030, 0.208] 50 inside | 0.046 [0.026, 0.081] 50 above | 10, 25, 45 | +0.040 [-0.021, +0.174] 50 above | 0.173 [0.033, 0.400] 50 above | 0.130 [0.019, 0.456] 50 inside | 0.168 [0.037, 0.484] 50 above |
| 0.9 | 0.152 [0.075, 0.253] 50 inside | 0.110 [0.059, 0.194] 50 inside | 0.044 [0.018, 0.111] 50 above | 15, 41, 45 | +0.021 [-0.033, +0.120] 50 above | 0.154 [0.049, 0.379] 50 above | 0.186 [0.066, 0.387] 50 below | 0.260 [0.037, 0.507] 50 above |

Parkin_Punctated: inside at F_W 0.5, 0.7, 0.9; above at 0, 0.01, 0.03, 0.1, 0.3; below at none; unavailable at none.
Barton/Kent/MPI: inside at F_W 0.5, 0.7, 0.9; above at 0, 0.01, 0.03, 0.1, 0.3; below at none; unavailable at none.
Painted: inside at F_W none; above at 0, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9; below at none; unavailable at none.

### R2. Level

Cells evaluable for the three common-class F values (E_L): 504 of 504. Cells with all three inside their bands (L): 5, resting on 100 distinct simulated histories. Of the 499 evaluable cells outside L, 499 have at least one value above its band and 0 have none above and at least one below.

Cells of L:

| span limit | window | coupling | F_W | rare classes inside (of 7) | rare classes unavailable | records lost from target: Parkin_Punctated, Barton/Kent/MPI, Painted | finite counts: F x 3, D_dist, A_lat x 3 |
|---|---|---|---|---|---|---|---|
| Ld/4 | 32LR | lat50N | 0.7 | 0 | 0 | 17, 22, 42 | 50, 50, 49, 50, 50, 50, 49 |
| Ld/4 | 32LR | lat100S | 0.9 | 0 | 0 | 16, 32, 45 | 50, 50, 48, 50, 50, 50, 48 |
| 1 | 32LR | c0 | 0.9 | 0 | 0 | 16, 32, 45 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 32LR | lat50N | 0.7 | 0 | 0 | 17, 22, 42 | 50, 50, 49, 50, 50, 50, 49 |
| 1 | 32LR | lat50S | 0.9 | 0 | 0 | 16, 32, 45 | 50, 50, 50, 50, 50, 50, 50 |

### R3. Level and space

Cells of L evaluable for the four spatial readouts (E_S): 5. Cells with the four also inside their bands (S): 0, resting on 0 distinct simulated histories. Cells of L not evaluable for space: 0.

Cells of S:

None.

### R5. By coupling

| coupling | cells | E_L | L | E_S | S |
|---|---|---|---|---|---|
| c0 | 72 | 72 | 1 | 1 | 0 |
| lat50N | 72 | 72 | 2 | 2 | 0 |
| lat100N | 72 | 72 | 0 | 0 | 0 |
| lat50S | 72 | 72 | 1 | 1 | 0 |
| lat100S | 72 | 72 | 1 | 1 | 0 |
| stepL | 72 | 72 | 0 | 0 | 0 |
| stepS | 72 | 72 | 0 | 0 | 0 |

**Statement (form b).** The record's three common-class F values lie inside their bands at 5 cells of the 504 evaluable (of 504 tested), listed below. At none of the 5 of those cells that are evaluable for space do its spatial readouts also lie inside theirs.

A band excludes about one record in five by construction, each cell carries up to fourteen statuses, and cells that share a window share simulated histories. The counts above are counts of grid cells, not of independent confirmations, and a small count can arise from the bands' own sampling.

### Recorded expectations

D_dist under c0 (72 cells): the record's value +0.360 is above the band at 56, inside at 16, below at 0 and unavailable at 0. Cell medians of simulated D_dist under c0 run from -0.007 to +0.319.

Median simulated F by span limit, at window 8 L_R and coupling c0 (median over 50 records, finite count in parentheses):

| F_W | Parkin_Punctated, L = 182 | Parkin_Punctated, L = 46 | Parkin_Punctated, L = 1 | Barton/Kent/MPI, L = 182 | Barton/Kent/MPI, L = 46 | Barton/Kent/MPI, L = 1 | Painted, L = 182 | Painted, L = 46 | Painted, L = 1 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.014 (50) | 0.016 (50) | 0.017 (50) | 0.016 (50) | 0.018 (50) | 0.018 (50) | 0.012 (50) | 0.015 (50) | 0.016 (50) |
| 0.01 | 0.015 (50) | 0.016 (50) | 0.017 (50) | 0.014 (50) | 0.016 (50) | 0.017 (50) | 0.014 (50) | 0.015 (50) | 0.017 (50) |
| 0.03 | 0.019 (50) | 0.019 (50) | 0.021 (50) | 0.016 (50) | 0.020 (50) | 0.022 (50) | 0.016 (50) | 0.018 (50) | 0.017 (50) |
| 0.1 | 0.021 (50) | 0.023 (50) | 0.023 (50) | 0.022 (50) | 0.023 (50) | 0.024 (50) | 0.024 (50) | 0.026 (50) | 0.028 (50) |
| 0.3 | 0.036 (50) | 0.035 (50) | 0.038 (50) | 0.037 (50) | 0.038 (50) | 0.036 (50) | 0.030 (50) | 0.033 (50) | 0.035 (50) |
| 0.5 | 0.059 (50) | 0.059 (50) | 0.063 (50) | 0.062 (50) | 0.063 (50) | 0.066 (50) | 0.040 (50) | 0.040 (50) | 0.041 (50) |
| 0.7 | 0.110 (50) | 0.113 (50) | 0.111 (50) | 0.089 (50) | 0.086 (50) | 0.081 (50) | 0.046 (50) | 0.045 (50) | 0.045 (50) |
| 0.9 | 0.152 (50) | 0.157 (50) | 0.161 (50) | 0.110 (50) | 0.106 (50) | 0.108 (50) | 0.044 (50) | 0.040 (50) | 0.041 (50) |

### Baseline comparison with the committed C3 table

Design d0 at the five earlier amplitudes, this run beside the committed table (`output/findings/per_type_generators_C3_primary.md`). The streams differ by design, so no equality is expected. The committed table reports inside or outside; this run's above and below both correspond to its outside. Agreement does not show that conclusions are robust to the seed set.

| F_W | class | this run: median [band] finite, lost, status | committed: median [band] finite, lost, status | status agrees |
|---|---|---|---|---|
| 0 | Parkin_Punctated | 0.014 [0.009, 0.023] 50, 0, above | 0.013 [0.006, 0.025] 50, 0, outside | yes |
| 0 | Barton/Kent/MPI | 0.016 [0.007, 0.022] 50, 0, above | 0.013 [0.006, 0.025] 50, 0, outside | yes |
| 0 | Painted | 0.012 [0.006, 0.023] 50, 0, above | 0.014 [0.007, 0.021] 50, 0, outside | yes |
| 0.01 | Parkin_Punctated | 0.015 [0.009, 0.024] 50, 0, above | 0.016 [0.008, 0.029] 50, 0, outside | yes |
| 0.01 | Barton/Kent/MPI | 0.014 [0.006, 0.026] 50, 0, above | 0.016 [0.007, 0.025] 50, 0, outside | yes |
| 0.01 | Painted | 0.014 [0.006, 0.022] 50, 0, above | 0.016 [0.009, 0.024] 50, 0, outside | yes |
| 0.03 | Parkin_Punctated | 0.019 [0.006, 0.032] 50, 0, above | 0.016 [0.008, 0.027] 50, 0, outside | yes |
| 0.03 | Barton/Kent/MPI | 0.016 [0.007, 0.038] 50, 0, above | 0.016 [0.008, 0.026] 50, 0, outside | yes |
| 0.03 | Painted | 0.016 [0.009, 0.033] 50, 0, above | 0.014 [0.007, 0.024] 50, 0, outside | yes |
| 0.1 | Parkin_Punctated | 0.021 [0.011, 0.036] 50, 0, above | 0.025 [0.011, 0.051] 50, 0, outside | yes |
| 0.1 | Barton/Kent/MPI | 0.022 [0.011, 0.043] 50, 0, above | 0.024 [0.012, 0.045] 50, 0, outside | yes |
| 0.1 | Painted | 0.024 [0.009, 0.040] 50, 12, above | 0.023 [0.011, 0.041] 50, 8, outside | yes |
| 0.3 | Parkin_Punctated | 0.036 [0.017, 0.108] 50, 1, above | 0.035 [0.011, 0.077] 50, 2, outside | yes |
| 0.3 | Barton/Kent/MPI | 0.037 [0.015, 0.116] 50, 11, above | 0.032 [0.016, 0.111] 50, 6, outside | yes |
| 0.3 | Painted | 0.030 [0.013, 0.114] 50, 30, above | 0.029 [0.014, 0.052] 50, 32, outside | yes |

Statuses that differ between the two runs: none.

## Secondary point

Frozen scales: L_d 12, L_R 40, burn-in 5000. Windows and stored lengths: 8LR W = 320, S = 340; 2LR W = 80, S = 100; 32LR W = 1280, S = 1300. Span limits: Ld = 12, Ld/4 = 3, 1 = 1.

Environment check: the committed C3 cell at F_W 0.3 was regenerated through `per_class_table.one_record_c3m` and equals the committed findings file to three decimals for the three common classes (Parkin_Punctated 0.033 [0.014, 0.092]; Barton/Kent/MPI 0.035 [0.013, 0.069]; Painted 0.041 [0.018, 0.082]), with the same losses from the target.

### R1. Amplitude at the baseline design d0 (span limit L_d, window 8 L_R, coupling c0)

Cells: median [10th, 90th percentile] finite count status. Lost: records of 50 in which the class's target share is exactly zero at the last stored step.

| F_W | F Parkin_Punctated (record 0.165) | F Barton/Kent/MPI (record 0.117) | F Painted (record 0.199) | lost: Parkin_Punctated, Barton/Kent/MPI, Painted | D_dist (record +0.360) | A_lat Parkin_Punctated (record 0.540) | A_lat Barton/Kent/MPI (record 0.057) | A_lat Painted (record 0.775) |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.017 [0.008, 0.025] 50 above | 0.016 [0.009, 0.026] 50 above | 0.018 [0.010, 0.025] 50 above | 0, 0, 0 | +0.074 [+0.004, +0.179] 50 above | 0.189 [0.053, 0.444] 50 above | 0.162 [0.051, 0.415] 50 inside | 0.145 [0.018, 0.500] 50 above |
| 0.01 | 0.017 [0.010, 0.025] 50 above | 0.018 [0.012, 0.030] 50 above | 0.015 [0.009, 0.025] 50 above | 0, 0, 0 | +0.053 [+0.001, +0.172] 50 above | 0.204 [0.079, 0.397] 50 above | 0.199 [0.050, 0.423] 50 inside | 0.255 [0.068, 0.429] 50 above |
| 0.03 | 0.020 [0.012, 0.034] 50 above | 0.020 [0.010, 0.034] 50 above | 0.019 [0.010, 0.037] 50 above | 0, 0, 0 | +0.052 [+0.007, +0.161] 50 above | 0.149 [0.035, 0.340] 50 above | 0.127 [0.036, 0.357] 50 inside | 0.191 [0.025, 0.462] 50 above |
| 0.1 | 0.023 [0.012, 0.045] 50 above | 0.024 [0.011, 0.058] 50 above | 0.023 [0.014, 0.056] 50 above | 0, 0, 9 | +0.054 [-0.039, +0.133] 50 above | 0.154 [0.036, 0.377] 50 above | 0.166 [0.031, 0.404] 50 inside | 0.169 [0.051, 0.327] 50 above |
| 0.3 | 0.050 [0.025, 0.099] 50 above | 0.046 [0.020, 0.123] 50 inside | 0.035 [0.018, 0.074] 50 above | 1, 6, 24 | +0.018 [-0.049, +0.130] 50 above | 0.160 [0.023, 0.412] 50 above | 0.182 [0.027, 0.432] 50 inside | 0.140 [0.037, 0.441] 50 above |
| 0.5 | 0.080 [0.027, 0.218] 50 inside | 0.106 [0.036, 0.190] 50 inside | 0.042 [0.021, 0.076] 50 above | 4, 24, 39 | +0.003 [-0.058, +0.079] 50 above | 0.138 [0.013, 0.280] 50 above | 0.153 [0.037, 0.288] 50 inside | 0.127 [0.015, 0.346] 50 above |
| 0.7 | 0.156 [0.028, 0.224] 50 inside | 0.120 [0.042, 0.204] 50 inside | 0.053 [0.027, 0.185] 50 above | 13, 25, 41 | +0.017 [-0.066, +0.083] 50 above | 0.137 [0.013, 0.332] 50 above | 0.131 [0.026, 0.302] 50 inside | 0.174 [0.023, 0.345] 50 above |
| 0.9 | 0.185 [0.058, 0.342] 50 inside | 0.141 [0.057, 0.311] 50 inside | 0.072 [0.016, 0.353] 50 inside | 25, 30, 40 | +0.007 [-0.058, +0.093] 50 above | 0.140 [0.031, 0.278] 50 above | 0.134 [0.031, 0.340] 50 inside | 0.166 [0.043, 0.333] 50 above |

Parkin_Punctated: inside at F_W 0.5, 0.7, 0.9; above at 0, 0.01, 0.03, 0.1, 0.3; below at none; unavailable at none.
Barton/Kent/MPI: inside at F_W 0.3, 0.5, 0.7, 0.9; above at 0, 0.01, 0.03, 0.1; below at none; unavailable at none.
Painted: inside at F_W 0.9; above at 0, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7; below at none; unavailable at none.

### R2. Level

Cells evaluable for the three common-class F values (E_L): 496 of 504. Cells with all three inside their bands (L): 45, resting on 150 distinct simulated histories. Of the 451 evaluable cells outside L, 446 have at least one value above its band and 5 have none above and at least one below.

Cells of L:

| span limit | window | coupling | F_W | rare classes inside (of 7) | rare classes unavailable | records lost from target: Parkin_Punctated, Barton/Kent/MPI, Painted | finite counts: F x 3, D_dist, A_lat x 3 |
|---|---|---|---|---|---|---|---|
| Ld | 8LR | c0 | 0.9 | 1 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 8LR | lat50N | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 8LR | lat100N | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 8LR | lat50S | 0.9 | 1 | 0 | 25, 30, 40 | 50, 50, 49, 50, 50, 50, 49 |
| Ld | 8LR | lat100S | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 8LR | stepL | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 48, 50, 50, 50, 48 |
| Ld | 8LR | stepS | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 32LR | c0 | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 32LR | c0 | 0.9 | 0 | 0 | 20, 30, 45 | 50, 50, 46, 50, 50, 50, 46 |
| Ld | 32LR | lat50N | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 32LR | lat50N | 0.9 | 0 | 0 | 20, 30, 45 | 50, 50, 46, 50, 50, 50, 46 |
| Ld | 32LR | lat100N | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 32LR | lat100N | 0.9 | 0 | 0 | 20, 30, 45 | 50, 50, 45, 50, 50, 50, 45 |
| Ld | 32LR | lat50S | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 32LR | lat100S | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld | 32LR | stepL | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 49, 50, 50, 50, 49 |
| Ld | 32LR | stepS | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 8LR | c0 | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 8LR | lat50N | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 8LR | lat100N | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 8LR | lat50S | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 49, 50, 50, 50, 49 |
| Ld/4 | 8LR | lat100S | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 49, 50, 50, 50, 49 |
| Ld/4 | 8LR | stepL | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 48, 50, 50, 50, 48 |
| Ld/4 | 8LR | stepS | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 49, 50, 50, 50, 49 |
| Ld/4 | 32LR | c0 | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 32LR | lat50N | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 32LR | lat100N | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 32LR | lat50S | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 32LR | lat100S | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| Ld/4 | 32LR | stepL | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 48, 50, 50, 50, 48 |
| Ld/4 | 32LR | stepS | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 49, 50, 50, 50, 49 |
| 1 | 8LR | c0 | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 8LR | lat50N | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 48, 50, 50, 50, 48 |
| 1 | 8LR | lat100N | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 8LR | lat50S | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 49, 50, 50, 50, 49 |
| 1 | 8LR | lat100S | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 49, 50, 50, 50, 49 |
| 1 | 8LR | stepL | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 48, 50, 50, 50, 48 |
| 1 | 8LR | stepS | 0.9 | 0 | 0 | 25, 30, 40 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 32LR | c0 | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 32LR | lat50N | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 32LR | lat100N | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 32LR | lat50S | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 32LR | lat100S | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 50, 50, 50, 50, 50 |
| 1 | 32LR | stepL | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 48, 50, 50, 50, 48 |
| 1 | 32LR | stepS | 0.7 | 0 | 0 | 17, 25, 39 | 50, 50, 49, 50, 50, 50, 49 |

### R3. Level and space

Cells of L evaluable for the four spatial readouts (E_S): 45. Cells with the four also inside their bands (S): 0, resting on 0 distinct simulated histories. Cells of L not evaluable for space: 0.

Cells of S:

None.

### R5. By coupling

| coupling | cells | E_L | L | E_S | S |
|---|---|---|---|---|---|
| c0 | 72 | 70 | 7 | 7 | 0 |
| lat50N | 72 | 72 | 7 | 7 | 0 |
| lat100N | 72 | 70 | 7 | 7 | 0 |
| lat50S | 72 | 70 | 6 | 6 | 0 |
| lat100S | 72 | 71 | 6 | 6 | 0 |
| stepL | 72 | 71 | 6 | 6 | 0 |
| stepS | 72 | 72 | 6 | 6 | 0 |

**Statement (form b).** The record's three common-class F values lie inside their bands at 45 cells of the 496 evaluable (of 504 tested), listed below. At none of the 45 of those cells that are evaluable for space do its spatial readouts also lie inside theirs.

A band excludes about one record in five by construction, each cell carries up to fourteen statuses, and cells that share a window share simulated histories. The counts above are counts of grid cells, not of independent confirmations, and a small count can arise from the bands' own sampling.

### Recorded expectations

D_dist under c0 (72 cells): the record's value +0.360 is above the band at 72, inside at 0, below at 0 and unavailable at 0. Cell medians of simulated D_dist under c0 run from -0.017 to +0.148.

Median simulated F by span limit, at window 8 L_R and coupling c0 (median over 50 records, finite count in parentheses):

| F_W | Parkin_Punctated, L = 12 | Parkin_Punctated, L = 3 | Parkin_Punctated, L = 1 | Barton/Kent/MPI, L = 12 | Barton/Kent/MPI, L = 3 | Barton/Kent/MPI, L = 1 | Painted, L = 12 | Painted, L = 3 | Painted, L = 1 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.017 (50) | 0.021 (50) | 0.022 (50) | 0.016 (50) | 0.020 (50) | 0.021 (50) | 0.018 (50) | 0.022 (50) | 0.024 (50) |
| 0.01 | 0.017 (50) | 0.021 (50) | 0.024 (50) | 0.018 (50) | 0.023 (50) | 0.024 (50) | 0.015 (50) | 0.020 (50) | 0.023 (50) |
| 0.03 | 0.020 (50) | 0.025 (50) | 0.029 (50) | 0.020 (50) | 0.026 (50) | 0.028 (50) | 0.019 (50) | 0.025 (50) | 0.029 (50) |
| 0.1 | 0.023 (50) | 0.027 (50) | 0.030 (50) | 0.024 (50) | 0.027 (50) | 0.028 (50) | 0.023 (50) | 0.029 (50) | 0.030 (50) |
| 0.3 | 0.050 (50) | 0.053 (50) | 0.051 (50) | 0.046 (50) | 0.053 (50) | 0.053 (50) | 0.035 (50) | 0.042 (50) | 0.042 (50) |
| 0.5 | 0.080 (50) | 0.086 (50) | 0.086 (50) | 0.106 (50) | 0.113 (50) | 0.111 (50) | 0.042 (50) | 0.047 (50) | 0.052 (50) |
| 0.7 | 0.156 (50) | 0.158 (50) | 0.162 (50) | 0.120 (50) | 0.120 (50) | 0.118 (50) | 0.053 (50) | 0.055 (50) | 0.057 (50) |
| 0.9 | 0.185 (50) | 0.175 (50) | 0.172 (50) | 0.141 (50) | 0.148 (50) | 0.145 (50) | 0.072 (50) | 0.070 (50) | 0.073 (50) |

### Baseline comparison with the committed C3 table

Design d0 at the five earlier amplitudes, this run beside the committed table (`output/findings/per_type_generators_C3_secondary.md`). The streams differ by design, so no equality is expected. The committed table reports inside or outside; this run's above and below both correspond to its outside. Agreement does not show that conclusions are robust to the seed set.

| F_W | class | this run: median [band] finite, lost, status | committed: median [band] finite, lost, status | status agrees |
|---|---|---|---|---|
| 0 | Parkin_Punctated | 0.017 [0.008, 0.025] 50, 0, above | 0.017 [0.010, 0.025] 50, 0, outside | yes |
| 0 | Barton/Kent/MPI | 0.016 [0.009, 0.026] 50, 0, above | 0.017 [0.010, 0.026] 50, 0, outside | yes |
| 0 | Painted | 0.018 [0.010, 0.025] 50, 0, above | 0.017 [0.009, 0.025] 50, 0, outside | yes |
| 0.01 | Parkin_Punctated | 0.017 [0.010, 0.025] 50, 0, above | 0.017 [0.011, 0.031] 50, 0, outside | yes |
| 0.01 | Barton/Kent/MPI | 0.018 [0.012, 0.030] 50, 0, above | 0.017 [0.011, 0.025] 50, 0, outside | yes |
| 0.01 | Painted | 0.015 [0.009, 0.025] 50, 0, above | 0.018 [0.009, 0.032] 50, 0, outside | yes |
| 0.03 | Parkin_Punctated | 0.020 [0.012, 0.034] 50, 0, above | 0.020 [0.010, 0.042] 50, 0, outside | yes |
| 0.03 | Barton/Kent/MPI | 0.020 [0.010, 0.034] 50, 0, above | 0.020 [0.010, 0.035] 50, 0, outside | yes |
| 0.03 | Painted | 0.019 [0.010, 0.037] 50, 0, above | 0.016 [0.009, 0.031] 50, 0, outside | yes |
| 0.1 | Parkin_Punctated | 0.023 [0.012, 0.045] 50, 0, above | 0.025 [0.015, 0.048] 50, 0, outside | yes |
| 0.1 | Barton/Kent/MPI | 0.024 [0.011, 0.058] 50, 0, above | 0.026 [0.014, 0.049] 50, 0, outside | yes |
| 0.1 | Painted | 0.023 [0.014, 0.056] 50, 9, above | 0.022 [0.014, 0.041] 50, 6, outside | yes |
| 0.3 | Parkin_Punctated | 0.050 [0.025, 0.099] 50, 1, above | 0.033 [0.014, 0.092] 50, 0, outside | yes |
| 0.3 | Barton/Kent/MPI | 0.046 [0.020, 0.123] 50, 6, inside | 0.035 [0.013, 0.069] 50, 3, outside | NO |
| 0.3 | Painted | 0.035 [0.018, 0.074] 50, 24, above | 0.041 [0.018, 0.082] 50, 26, outside | yes |

Statuses that differ between the two runs: F_W 0.3, Barton/Kent/MPI.

## What these tables cannot give

Joint reproduction of the record by any simulated record; the sufficiency or necessity of a timing arrangement; whether the sites differ in date; type origination, which this generator does not have; identified parameters, since results hold at two operating points; robustness to the seed set, since each cell is 50 records from one set of streams; a test of boundaries; and the fitted spatial share of the A0 model, which is not computed on simulated records (spec section 10).

Per-cell values: `output/findings/drift_extension_cells.csv`. Per-record values: `output/probe_latent_groups/drift_extension_<point>.csv` (not committed).
