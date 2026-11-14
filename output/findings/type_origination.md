# Type origination during the occupation (generator C4): exploratory

Produced by `.venv-np254/bin/python analyses/per_type/type_origination.py` at commit 3db81f4f9947; NumPy 2.5.4; data hash 3b5e1cfabcad47f4; run time primary point 189 s, secondary point 48 s. Design: `docs/superpowers/specs/2026-11-14-type-origination-generator-design.md`, draft 1. Exploratory by the author's decision of 2026-11-14. The spec was not audited. The build and the first run were audited by Codex on 2026-11-14, after the run (`docs/audits/codex_audit_2026-11-14_build_ledger.md`). Each cell is 50 records from one seed set.

Generator C4 is the neutral copying model of C3 with the regional innovation target given by a finite pool, resampled at every stored step, in which each draw is a new type with probability theta / (2 N_R). theta = 0 has no origination. The Ewens estimate from the record's richness (ten classes in 14,101 sherds) is theta = 0.98, a reference point only, since archaeological types are not neutral variants. Grid per point: 2 span limits x 2 windows x 3 date arrangements x 5 amplitudes x 5 values of theta = 300 cells.

Readouts are by rank of pooled count, because a simulated new type has no counterpart among the record's named classes. Every comparison places the record's value beside the 10th to 90th percentile band of the cell's simulated records. These are descriptive bands, not posterior intervals, and every statement is of marginal coverage: each named summary lies inside its own band. Nothing here shows that types originated during the occupation at any intensity.

## The record, by rank

F of ranks 1 to 3: 0.165, 0.117, 0.199. Sites present for ranks 4 to 10: 22, 15, 14, 15, 7, 10, 8. Richness 10; share beyond rank 10: 0.000. D_dist +0.360; latitude rank correlations of ranks 1 to 3: +0.540, +0.057, -0.775.

## Primary point

Pool sizes N_R (realized F over the window in parentheses): 8LR F_W 0.1: 22172 (0.100), F_W 0.3: 6550 (0.300), F_W 0.5: 3371 (0.500), F_W 0.7: 1941 (0.700), F_W 0.9: 1015 (0.900); 32LR F_W 0.1: 88686 (0.100), F_W 0.3: 26198 (0.300), F_W 0.5: 13481 (0.500), F_W 0.7: 7761 (0.700), F_W 0.9: 4059 (0.900).

### Starting classes in the pool

The pool starts from the record's pooled shares rounded to N_R draws by largest remainder. A class allocated zero draws is absent from the start and, since a lost type cannot return, is absent from every record of that cell at every theta. This is absence from the start, distinct from loss by drift during the window. Where a class is listed, a theta = 0 record has at most the stated number of classes, so the occupancy of the lowest ranks is zero by construction.

| window | F_W | N_R | starting classes present | classes rounded to zero draws |
|---|---|---|---|---|
| 8LR | 0.1 | 22172 | 10 | none |
| 8LR | 0.3 | 6550 | 10 | none |
| 8LR | 0.5 | 3371 | 10 | none |
| 8LR | 0.7 | 1941 | 10 | none |
| 8LR | 0.9 | 1015 | 10 | none |
| 32LR | 0.1 | 88686 | 10 | none |
| 32LR | 0.3 | 26198 | 10 | none |
| 32LR | 0.5 | 13481 | 10 | none |
| 32LR | 0.7 | 7761 | 10 | none |
| 32LR | 0.9 | 4059 | 10 | none |

### Sets by origination intensity

Cells of 60 per theta at which the record's summaries each lie inside their own band. L: F of ranks 1 to 3. R: sites present for all of ranks 4 to 10. R5: at least five of those seven (a threshold chosen after the drift-extension run). K: richness and the share beyond rank 10. S3: D_dist and the latitude trends of ranks 1 and 3.

| theta | L | R | R5 | K | S3 | L and R | all of L, R, K, S3 | L and R5 | L and R5 and K | L and R5 and S3 | all of L, R5, K, S3 | histories behind L and R5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 10 | 0 | 0 | 0 | 31 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0.3 | 17 | 0 | 1 | 33 | 26 | 0 | 0 | 1 | 0 | 0 | 0 | 50 |
| 1 | 22 | 20 | 30 | 0 | 17 | 3 | 0 | 8 | 0 | 0 | 0 | 100 |
| 3 | 22 | 29 | 60 | 0 | 10 | 4 | 0 | 22 | 0 | 0 | 0 | 250 |
| 10 | 8 | 3 | 45 | 0 | 1 | 1 | 0 | 5 | 0 | 0 | 0 | 100 |

**Statement (form b).** The level of variation and the occupancy of all seven rare ranks both lie inside their bands at 8 cells. None of those cells also has richness and the three spatial readouts inside their bands: richness is outside at 8 of them and space at 8.

With the threshold relaxed to at least five of the seven rare ranks (chosen after the drift-extension run, so not part of the statement form): level and rare ranks together at 36 cells, of which 0 also have richness inside, 0 also have the three spatial readouts inside, and 0 have both. The cells are listed below.

Counts are of grid cells from one seed set. The six designs of a window, amplitude and theta share simulated histories, and a band excludes about one record in five by construction, so a small count can arise from the bands' own sampling.

### Cells with level and at least five rare ranks, or level and space

| span | window | dates | F_W | theta | rare ranks inside | richness: median [band] (record 10) | beyond rank 10: median | K | S3 | new types among top 3 / top 10 (median) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 32LR | c0 | 0.7 | 3 | 6 | 108 [97, 119] | 0.149 | no | no | 1 / 7 |
| 1 | 32LR | c0 | 0.9 | 3 | 5 | 146 [127, 158] | 0.229 | no | no | 2 / 8 |
| 1 | 32LR | lat100N | 0.5 | 1 | 3 | 31 [25, 37] | 0.011 | no | yes | 1 / 7 |
| 1 | 32LR | lat100N | 0.5 | 3 | 7 | 76 [67, 93] | 0.087 | no | no | 1 / 7 |
| 1 | 32LR | lat100N | 0.5 | 10 | 6 | 234 [217, 262] | 0.402 | no | no | 2 / 8 |
| 1 | 32LR | lat100N | 0.7 | 0 | 0 | 5 [4, 6] | 0.000 | no | yes | 0 / 0 |
| 1 | 32LR | lat100N | 0.7 | 0.3 | 0 | 14 [11, 20] | 0.001 | no | yes | 0 / 6 |
| 1 | 32LR | lat100N | 0.7 | 1 | 5 | 38 [32, 45] | 0.019 | no | no | 1 / 7 |
| 1 | 32LR | lat100N | 0.7 | 3 | 6 | 102 [91, 113] | 0.153 | no | no | 2 / 8 |
| 1 | 32LR | lat100N | 0.9 | 0 | 0 | 5 [3, 6] | 0.000 | no | yes | 0 / 0 |
| 1 | 32LR | lat100N | 0.9 | 0.3 | 0 | 18 [13, 22] | 0.002 | no | yes | 1 / 6 |
| 1 | 32LR | lat100N | 0.9 | 3 | 5 | 134 [119, 152] | 0.223 | no | no | 2 / 8 |
| 1 | 32LR | lat100S | 0.7 | 3 | 5 | 101 [92, 109] | 0.136 | no | no | 1 / 7 |
| 1 | 32LR | lat100S | 0.7 | 10 | 5 | 309 [286, 328] | 0.487 | no | no | 2 / 8 |
| 1 | 32LR | lat100S | 0.9 | 0 | 0 | 5 [4, 6] | 0.000 | no | yes | 0 / 0 |
| 1 | 32LR | lat100S | 0.9 | 3 | 5 | 135 [118, 150] | 0.204 | no | no | 2 / 8 |
| 1 | 8LR | c0 | 0.9 | 1 | 7 | 52 [44, 62] | 0.045 | no | no | 1 / 7 |
| 1 | 8LR | c0 | 0.9 | 3 | 6 | 132 [119, 151] | 0.224 | no | no | 1 / 7 |
| 1 | 8LR | lat100N | 0.7 | 3 | 6 | 88 [74, 98] | 0.128 | no | no | 1 / 7 |
| 1 | 8LR | lat100N | 0.9 | 1 | 7 | 46 [39, 53] | 0.039 | no | no | 1 / 7 |
| 1 | 8LR | lat100N | 0.9 | 3 | 6 | 122 [104, 136] | 0.222 | no | no | 2 / 8 |
| 1 | 8LR | lat100S | 0.7 | 3 | 6 | 88 [77, 101] | 0.124 | no | no | 1 / 7 |
| 1 | 8LR | lat100S | 0.9 | 1 | 6 | 45 [38, 55] | 0.035 | no | no | 1 / 7 |
| Ld | 32LR | c0 | 0.5 | 10 | 7 | 276 [257, 299] | 0.390 | no | no | 1 / 7 |
| Ld | 32LR | c0 | 0.7 | 3 | 7 | 116 [106, 131] | 0.151 | no | no | 1 / 7 |
| Ld | 32LR | c0 | 0.9 | 3 | 5 | 158 [138, 176] | 0.228 | no | no | 2 / 8 |
| Ld | 32LR | lat100N | 0.5 | 3 | 7 | 84 [74, 100] | 0.084 | no | no | 1 / 7 |
| Ld | 32LR | lat100N | 0.5 | 10 | 5 | 262 [238, 284] | 0.403 | no | no | 1 / 8 |
| Ld | 32LR | lat100N | 0.7 | 0 | 0 | 5 [4, 6] | 0.000 | no | yes | 0 / 0 |
| Ld | 32LR | lat100N | 0.7 | 0.3 | 0 | 16 [12, 20] | 0.001 | no | yes | 0 / 6 |
| Ld | 32LR | lat100N | 0.7 | 3 | 7 | 109 [98, 125] | 0.151 | no | no | 2 / 8 |
| Ld | 32LR | lat100N | 0.9 | 0 | 0 | 5 [3, 6] | 0.000 | no | yes | 0 / 0 |
| Ld | 32LR | lat100N | 0.9 | 0.3 | 0 | 20 [16, 24] | 0.002 | no | yes | 1 / 7 |
| Ld | 32LR | lat100N | 0.9 | 3 | 5 | 148 [136, 160] | 0.222 | no | no | 2 / 8 |
| Ld | 32LR | lat100S | 0.7 | 0 | 0 | 5 [4, 7] | 0.000 | no | yes | 0 / 0 |
| Ld | 32LR | lat100S | 0.7 | 1 | 5 | 41 [34, 49] | 0.023 | no | no | 1 / 7 |
| Ld | 32LR | lat100S | 0.7 | 3 | 6 | 111 [99, 122] | 0.132 | no | no | 1 / 7 |
| Ld | 32LR | lat100S | 0.7 | 10 | 6 | 334 [308, 359] | 0.477 | no | no | 2 / 8 |
| Ld | 32LR | lat100S | 0.9 | 0 | 0 | 5 [4, 6] | 0.000 | no | yes | 0 / 0 |
| Ld | 32LR | lat100S | 0.9 | 3 | 5 | 148 [129, 163] | 0.206 | no | no | 2 / 8 |
| Ld | 8LR | c0 | 0.7 | 3 | 6 | 106 [94, 120] | 0.137 | no | no | 1 / 7 |
| Ld | 8LR | c0 | 0.9 | 1 | 7 | 56 [45, 65] | 0.041 | no | no | 1 / 7 |
| Ld | 8LR | c0 | 0.9 | 3 | 5 | 147 [132, 161] | 0.219 | no | no | 1 / 7 |
| Ld | 8LR | lat100N | 0.7 | 3 | 5 | 97 [83, 109] | 0.124 | no | no | 1 / 7 |
| Ld | 8LR | lat100N | 0.9 | 1 | 6 | 50 [41, 60] | 0.039 | no | no | 1 / 7 |
| Ld | 8LR | lat100N | 0.9 | 3 | 5 | 136 [118, 151] | 0.219 | no | no | 2 / 8 |
| Ld | 8LR | lat100S | 0.9 | 0.3 | 5 | 18 [15, 24] | 0.002 | no | no | 1 / 6 |
| Ld | 8LR | lat100S | 0.9 | 1 | 6 | 51 [41, 59] | 0.035 | no | no | 1 / 7 |

### Baseline design (span limit L_d, window 8 L_R, dates independent of location) by amplitude and theta

Cells: median [10th, 90th percentile] status. Rare ranks: number of ranks 4 to 10 whose record occupancy is inside its band.

| F_W | theta | F rank 1 (record 0.165) | F rank 2 (record 0.117) | F rank 3 (record 0.199) | rare ranks inside | sites, rank 4 (record 22) | sites, rank 7 (record 15) | sites, rank 10 (record 8) | richness (record 10) | beyond rank 10 | types originated | new among top 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0 | 0.020 [0.009, 0.045] above | 0.024 [0.014, 0.040] above | 0.020 [0.009, 0.053] above | 2 | 23 [12, 28] inside | 2 [0, 7] above | 0 [0, 0] above | 7 [6, 9] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.1 | 0.3 | 0.018 [0.007, 0.033] above | 0.021 [0.011, 0.046] above | 0.025 [0.014, 0.047] above | 3 | 22 [9, 27] inside | 4 [2, 8] above | 1 [0, 3] above | 10 [9, 13] inside | 0.000 [0.000, 0.000] inside | 390 | 2 |
| 0.1 | 1 | 0.024 [0.011, 0.042] above | 0.022 [0.009, 0.045] above | 0.027 [0.011, 0.058] above | 4 | 21 [14, 27] inside | 6 [3, 12] above | 3 [1, 7] above | 18 [13, 23] below | 0.001 [0.000, 0.006] below | 1300 | 4 |
| 0.1 | 3 | 0.021 [0.011, 0.042] above | 0.024 [0.008, 0.061] above | 0.019 [0.010, 0.032] above | 7 | 24 [14, 27] inside | 10 [6, 20] inside | 7 [4, 14] inside | 37 [30, 44] below | 0.016 [0.008, 0.029] below | 3889 | 4 |
| 0.1 | 10 | 0.023 [0.015, 0.061] above | 0.041 [0.015, 0.079] above | 0.024 [0.011, 0.050] above | 6 | 22 [15, 27] inside | 16 [10, 23] inside | 12 [7, 18] inside | 106 [95, 117] below | 0.106 [0.073, 0.144] below | 12998 | 6 |
| 0.3 | 0 | 0.046 [0.016, 0.123] above | 0.030 [0.011, 0.098] above | 0.036 [0.021, 0.080] above | 2 | 12 [3, 26] inside | 1 [0, 4] above | 0 [0, 0] above | 7 [6, 8] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.3 | 0.3 | 0.034 [0.014, 0.081] above | 0.036 [0.018, 0.081] above | 0.034 [0.012, 0.086] above | 2 | 15 [7, 25] inside | 4 [1, 8] above | 2 [1, 3] above | 13 [10, 16] inside | 0.000 [0.000, 0.002] inside | 390 | 3 |
| 0.3 | 1 | 0.028 [0.014, 0.078] above | 0.039 [0.013, 0.097] above | 0.038 [0.016, 0.094] above | 7 | 22 [15, 27] inside | 8 [3, 15] inside | 6 [2, 9] inside | 26 [20, 32] below | 0.006 [0.002, 0.016] below | 1310 | 5 |
| 0.3 | 3 | 0.026 [0.010, 0.082] above | 0.044 [0.014, 0.108] above | 0.045 [0.026, 0.104] above | 6 | 22 [15, 28] inside | 16 [9, 22] inside | 10 [5, 17] inside | 62 [56, 74] below | 0.051 [0.032, 0.077] below | 3904 | 6 |
| 0.3 | 10 | 0.087 [0.034, 0.169] inside | 0.066 [0.020, 0.117] above | 0.037 [0.018, 0.092] above | 5 | 22 [17, 27] inside | 18 [11, 23] inside | 15 [9, 20] below | 190 [172, 208] below | 0.236 [0.199, 0.303] below | 13028 | 7 |
| 0.5 | 0 | 0.058 [0.016, 0.167] inside | 0.072 [0.022, 0.161] inside | 0.044 [0.018, 0.090] above | 0 | 9 [3, 20] above | 1 [0, 4] above | 0 [0, 0] above | 7 [5, 8] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.5 | 0.3 | 0.053 [0.015, 0.154] above | 0.068 [0.021, 0.132] inside | 0.054 [0.020, 0.145] above | 4 | 14 [6, 24] inside | 6 [2, 11] above | 2 [1, 5] above | 14 [12, 18] below | 0.001 [0.000, 0.003] below | 396 | 4 |
| 0.5 | 1 | 0.056 [0.021, 0.129] above | 0.075 [0.027, 0.145] inside | 0.048 [0.020, 0.102] above | 7 | 19 [12, 28] inside | 10 [4, 17] inside | 6 [3, 11] inside | 34 [29, 42] below | 0.015 [0.007, 0.030] below | 1296 | 5 |
| 0.5 | 3 | 0.050 [0.014, 0.128] above | 0.083 [0.025, 0.175] inside | 0.069 [0.030, 0.138] above | 6 | 22 [13, 27] inside | 16 [8, 23] inside | 11 [5, 17] inside | 83 [73, 91] below | 0.078 [0.051, 0.111] below | 3902 | 7 |
| 0.5 | 10 | 0.132 [0.059, 0.227] inside | 0.093 [0.037, 0.172] inside | 0.050 [0.021, 0.122] above | 4 | 22 [15, 26] inside | 18 [13, 23] inside | 14 [10, 21] below | 260 [236, 281] below | 0.372 [0.296, 0.461] below | 13000 | 7 |
| 0.7 | 0 | 0.137 [0.042, 0.227] inside | 0.078 [0.036, 0.176] inside | 0.046 [0.025, 0.106] above | 0 | 7 [1, 20] above | 1 [0, 4] above | 0 [0, 0] above | 7 [4, 8] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.7 | 0.3 | 0.089 [0.028, 0.229] inside | 0.099 [0.021, 0.224] inside | 0.055 [0.026, 0.157] above | 3 | 14 [7, 22] inside | 5 [3, 11] above | 2 [1, 5] above | 18 [14, 21] below | 0.002 [0.000, 0.005] below | 389 | 5 |
| 0.7 | 1 | 0.058 [0.020, 0.177] inside | 0.094 [0.040, 0.247] inside | 0.068 [0.033, 0.155] above | 7 | 18 [11, 27] inside | 10 [5, 18] inside | 7 [3, 12] inside | 39 [32, 48] below | 0.017 [0.008, 0.033] below | 1302 | 6 |
| 0.7 | 3 | 0.085 [0.024, 0.192] inside | 0.079 [0.040, 0.207] inside | 0.086 [0.044, 0.211] inside | 6 | 21 [14, 28] inside | 15 [9, 22] inside | 13 [8, 17] inside | 106 [94, 120] below | 0.137 [0.084, 0.207] below | 3912 | 7 |
| 0.7 | 10 | 0.148 [0.060, 0.233] inside | 0.096 [0.048, 0.163] inside | 0.065 [0.027, 0.129] above | 5 | 20 [13, 24] inside | 16 [11, 21] inside | 15 [10, 19] below | 334 [304, 359] below | 0.466 [0.389, 0.543] below | 13020 | 7 |
| 0.9 | 0 | 0.194 [0.065, 0.293] inside | 0.149 [0.048, 0.220] inside | 0.041 [0.022, 0.100] above | 0 | 7 [4, 14] above | 0 [0, 3] above | 0 [0, 0] above | 6 [5, 8] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.9 | 0.3 | 0.125 [0.038, 0.285] inside | 0.106 [0.039, 0.286] inside | 0.087 [0.033, 0.219] inside | 3 | 13 [7, 21] above | 7 [3, 9] above | 3 [1, 7] above | 21 [17, 25] below | 0.003 [0.001, 0.009] below | 391 | 5 |
| 0.9 | 1 | 0.091 [0.035, 0.316] inside | 0.139 [0.070, 0.255] inside | 0.122 [0.050, 0.247] inside | 7 | 18 [10, 24] inside | 11 [7, 17] inside | 9 [4, 13] inside | 56 [45, 65] below | 0.041 [0.025, 0.070] below | 1303 | 7 |
| 0.9 | 3 | 0.146 [0.064, 0.244] inside | 0.116 [0.061, 0.216] inside | 0.109 [0.042, 0.216] inside | 5 | 20 [14, 26] inside | 17 [10, 21] inside | 14 [10, 18] below | 147 [132, 161] below | 0.219 [0.158, 0.290] below | 3894 | 7 |
| 0.9 | 10 | 0.157 [0.073, 0.238] inside | 0.100 [0.059, 0.200] inside | 0.053 [0.034, 0.106] above | 5 | 17 [14, 22] inside | 16 [12, 20] inside | 14 [9, 20] below | 443 [415, 468] below | 0.577 [0.502, 0.642] below | 13012 | 8 |

### Calibration: theta = 0 beside generator C3

The comparison of C4 at theta = 0 with C3 on identical ranked readouts, with medians and bands for every statistic, is in `output/findings/type_origination_calibration.md` (`analyses/per_type/type_origination_calibration.py`).

## Secondary point

Pool sizes N_R (realized F over the window in parentheses): 8LR F_W 0.1: 3038 (0.100), F_W 0.3: 898 (0.300), F_W 0.5: 462 (0.500), F_W 0.7: 266 (0.700), F_W 0.9: 139 (0.901); 32LR F_W 0.1: 12149 (0.100), F_W 0.3: 3589 (0.300), F_W 0.5: 1847 (0.500), F_W 0.7: 1064 (0.700), F_W 0.9: 556 (0.900).

### Starting classes in the pool

The pool starts from the record's pooled shares rounded to N_R draws by largest remainder. A class allocated zero draws is absent from the start and, since a lost type cannot return, is absent from every record of that cell at every theta. This is absence from the start, distinct from loss by drift during the window. Where a class is listed, a theta = 0 record has at most the stated number of classes, so the occupancy of the lowest ranks is zero by construction.

| window | F_W | N_R | starting classes present | classes rounded to zero draws |
|---|---|---|---|---|
| 8LR | 0.1 | 3038 | 10 | none |
| 8LR | 0.3 | 898 | 10 | none |
| 8LR | 0.5 | 462 | 10 | none |
| 8LR | 0.7 | 266 | 9 | Hull_Engraved |
| 8LR | 0.9 | 139 | 8 | Vernon_Paul_Applique, Hull_Engraved |
| 32LR | 0.1 | 12149 | 10 | none |
| 32LR | 0.3 | 3589 | 10 | none |
| 32LR | 0.5 | 1847 | 10 | none |
| 32LR | 0.7 | 1064 | 10 | none |
| 32LR | 0.9 | 556 | 10 | none |

### Sets by origination intensity

Cells of 60 per theta at which the record's summaries each lie inside their own band. L: F of ranks 1 to 3. R: sites present for all of ranks 4 to 10. R5: at least five of those seven (a threshold chosen after the drift-extension run). K: richness and the share beyond rank 10. S3: D_dist and the latitude trends of ranks 1 and 3.

| theta | L | R | R5 | K | S3 | L and R | all of L, R, K, S3 | L and R5 | L and R5 and K | L and R5 and S3 | all of L, R5, K, S3 | histories behind L and R5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 11 | 0 | 0 | 0 | 33 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0.3 | 24 | 0 | 0 | 48 | 37 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | 23 | 3 | 18 | 4 | 25 | 3 | 0 | 9 | 0 | 2 | 0 | 100 |
| 3 | 24 | 27 | 57 | 0 | 8 | 2 | 0 | 22 | 0 | 0 | 0 | 250 |
| 10 | 7 | 14 | 49 | 0 | 3 | 0 | 0 | 4 | 0 | 0 | 0 | 100 |

**Statement (form b).** The level of variation and the occupancy of all seven rare ranks both lie inside their bands at 5 cells. None of those cells also has richness and the three spatial readouts inside their bands: richness is outside at 5 of them and space at 4.

With the threshold relaxed to at least five of the seven rare ranks (chosen after the drift-extension run, so not part of the statement form): level and rare ranks together at 35 cells, of which 0 also have richness inside, 2 also have the three spatial readouts inside, and 0 have both. The cells are listed below.

Counts are of grid cells from one seed set. The six designs of a window, amplitude and theta share simulated histories, and a band excludes about one record in five by construction, so a small count can arise from the bands' own sampling.

### Cells with level and at least five rare ranks, or level and space

| span | window | dates | F_W | theta | rare ranks inside | richness: median [band] (record 10) | beyond rank 10: median | K | S3 | new types among top 3 / top 10 (median) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 32LR | c0 | 0.5 | 3 | 6 | 62 [57, 73] | 0.084 | no | no | 1 / 7 |
| 1 | 32LR | c0 | 0.7 | 3 | 6 | 82 [70, 96] | 0.155 | no | no | 2 / 7 |
| 1 | 32LR | c0 | 0.7 | 10 | 5 | 242 [224, 258] | 0.513 | no | no | 1 / 8 |
| 1 | 32LR | lat100N | 0.5 | 0.3 | 0 | 10 [8, 15] | 0.000 | yes | yes | 0 / 4 |
| 1 | 32LR | lat100N | 0.5 | 1 | 3 | 23 [19, 30] | 0.009 | no | yes | 1 / 7 |
| 1 | 32LR | lat100N | 0.5 | 3 | 5 | 63 [51, 72] | 0.080 | no | no | 1 / 7 |
| 1 | 32LR | lat100N | 0.7 | 0 | 0 | 5 [3, 7] | 0.000 | no | yes | 0 / 0 |
| 1 | 32LR | lat100N | 0.7 | 0.3 | 0 | 14 [10, 17] | 0.001 | yes | yes | 1 / 6 |
| 1 | 32LR | lat100N | 0.7 | 1 | 2 | 31 [24, 37] | 0.019 | no | yes | 1 / 7 |
| 1 | 32LR | lat100N | 0.7 | 3 | 5 | 82 [71, 96] | 0.143 | no | no | 2 / 8 |
| 1 | 32LR | lat100N | 0.9 | 0.3 | 0 | 14 [10, 18] | 0.001 | yes | yes | 1 / 6 |
| 1 | 32LR | lat100S | 0.5 | 0 | 0 | 5 [4, 7] | 0.000 | no | yes | 0 / 0 |
| 1 | 32LR | lat100S | 0.5 | 0.3 | 0 | 10 [7, 15] | 0.000 | yes | yes | 0 / 4 |
| 1 | 32LR | lat100S | 0.7 | 0.3 | 0 | 13 [10, 17] | 0.001 | yes | yes | 0 / 6 |
| 1 | 32LR | lat100S | 0.7 | 1 | 3 | 32 [24, 37] | 0.019 | no | yes | 1 / 7 |
| 1 | 32LR | lat100S | 0.7 | 3 | 5 | 80 [68, 93] | 0.139 | no | no | 1 / 7 |
| 1 | 32LR | lat100S | 0.9 | 0.3 | 0 | 15 [11, 19] | 0.002 | no | yes | 1 / 6 |
| 1 | 8LR | c0 | 0.7 | 3 | 7 | 72 [62, 83] | 0.134 | no | no | 1 / 7 |
| 1 | 8LR | c0 | 0.9 | 1 | 6 | 36 [29, 40] | 0.037 | no | no | 1 / 7 |
| 1 | 8LR | c0 | 0.9 | 3 | 6 | 98 [83, 108] | 0.210 | no | no | 2 / 8 |
| 1 | 8LR | lat100N | 0.7 | 0 | 0 | 5 [4, 7] | 0.000 | no | yes | 0 / 0 |
| 1 | 8LR | lat100N | 0.7 | 0.3 | 1 | 12 [8, 16] | 0.000 | yes | yes | 1 / 5 |
| 1 | 8LR | lat100N | 0.7 | 3 | 7 | 68 [57, 79] | 0.134 | no | no | 1 / 8 |
| 1 | 8LR | lat100N | 0.9 | 0.3 | 1 | 13 [9, 18] | 0.001 | yes | yes | 1 / 6 |
| 1 | 8LR | lat100N | 0.9 | 1 | 7 | 34 [27, 39] | 0.033 | no | no | 1 / 7 |
| 1 | 8LR | lat100N | 0.9 | 3 | 5 | 90 [80, 103] | 0.211 | no | no | 2 / 8 |
| 1 | 8LR | lat100S | 0.7 | 1 | 5 | 26 [21, 33] | 0.017 | no | yes | 1 / 6 |
| 1 | 8LR | lat100S | 0.7 | 3 | 5 | 67 [59, 79] | 0.133 | no | no | 1 / 7 |
| 1 | 8LR | lat100S | 0.9 | 0.3 | 1 | 14 [10, 19] | 0.001 | yes | yes | 1 / 5 |
| 1 | 8LR | lat100S | 0.9 | 1 | 5 | 34 [27, 39] | 0.032 | no | no | 1 / 7 |
| 1 | 8LR | lat100S | 0.9 | 3 | 5 | 92 [81, 103] | 0.200 | no | no | 2 / 8 |
| Ld | 32LR | c0 | 0.5 | 3 | 6 | 78 [69, 88] | 0.085 | no | no | 1 / 7 |
| Ld | 32LR | c0 | 0.7 | 3 | 6 | 102 [87, 116] | 0.151 | no | no | 2 / 7 |
| Ld | 32LR | c0 | 0.7 | 10 | 5 | 302 [281, 331] | 0.515 | no | no | 1 / 8 |
| Ld | 32LR | lat100N | 0.5 | 1 | 3 | 28 [23, 35] | 0.010 | no | yes | 1 / 7 |
| Ld | 32LR | lat100N | 0.5 | 10 | 6 | 227 [208, 255] | 0.413 | no | no | 2 / 8 |
| Ld | 32LR | lat100N | 0.7 | 0.3 | 0 | 15 [12, 20] | 0.001 | no | yes | 1 / 6 |
| Ld | 32LR | lat100N | 0.7 | 1 | 2 | 38 [28, 45] | 0.021 | no | yes | 1 / 7 |
| Ld | 32LR | lat100N | 0.7 | 3 | 5 | 100 [86, 110] | 0.145 | no | no | 2 / 8 |
| Ld | 32LR | lat100N | 0.9 | 0.3 | 0 | 18 [13, 22] | 0.002 | no | yes | 1 / 6 |
| Ld | 32LR | lat100N | 0.9 | 3 | 5 | 130 [112, 144] | 0.222 | no | no | 2 / 8 |
| Ld | 32LR | lat100S | 0.5 | 0.3 | 0 | 12 [7, 18] | 0.000 | yes | yes | 0 / 5 |
| Ld | 32LR | lat100S | 0.7 | 0.3 | 0 | 16 [11, 20] | 0.001 | no | yes | 0 / 6 |
| Ld | 32LR | lat100S | 0.7 | 1 | 3 | 39 [29, 45] | 0.020 | no | yes | 1 / 7 |
| Ld | 32LR | lat100S | 0.7 | 3 | 6 | 99 [83, 113] | 0.137 | no | no | 1 / 7 |
| Ld | 32LR | lat100S | 0.7 | 10 | 5 | 293 [277, 318] | 0.518 | no | no | 2 / 8 |
| Ld | 32LR | lat100S | 0.9 | 0.3 | 0 | 17 [13, 22] | 0.002 | no | yes | 1 / 6 |
| Ld | 8LR | c0 | 0.7 | 3 | 6 | 88 [75, 98] | 0.141 | no | no | 1 / 7 |
| Ld | 8LR | c0 | 0.9 | 1 | 6 | 42 [36, 49] | 0.037 | no | no | 1 / 7 |
| Ld | 8LR | c0 | 0.9 | 3 | 5 | 116 [105, 133] | 0.210 | no | no | 2 / 8 |
| Ld | 8LR | lat100N | 0.7 | 0 | 0 | 5 [4, 7] | 0.000 | no | yes | 0 / 0 |
| Ld | 8LR | lat100N | 0.7 | 0.3 | 2 | 13 [10, 18] | 0.000 | yes | yes | 1 / 5 |
| Ld | 8LR | lat100N | 0.7 | 1 | 6 | 32 [23, 38] | 0.016 | no | no | 1 / 7 |
| Ld | 8LR | lat100N | 0.7 | 3 | 6 | 84 [72, 92] | 0.138 | no | no | 1 / 8 |
| Ld | 8LR | lat100N | 0.9 | 0.3 | 2 | 16 [12, 21] | 0.002 | no | yes | 1 / 5 |
| Ld | 8LR | lat100N | 0.9 | 1 | 7 | 40 [33, 47] | 0.034 | no | no | 1 / 7 |
| Ld | 8LR | lat100N | 0.9 | 3 | 6 | 110 [97, 124] | 0.208 | no | no | 2 / 8 |
| Ld | 8LR | lat100S | 0.7 | 1 | 7 | 31 [27, 39] | 0.018 | no | yes | 1 / 6 |
| Ld | 8LR | lat100S | 0.7 | 3 | 6 | 83 [70, 94] | 0.130 | no | no | 1 / 7 |
| Ld | 8LR | lat100S | 0.9 | 0.3 | 2 | 16 [12, 20] | 0.001 | no | yes | 1 / 5 |
| Ld | 8LR | lat100S | 0.9 | 1 | 6 | 41 [31, 47] | 0.033 | no | no | 1 / 7 |
| Ld | 8LR | lat100S | 0.9 | 3 | 6 | 112 [97, 125] | 0.208 | no | no | 2 / 8 |

### Baseline design (span limit L_d, window 8 L_R, dates independent of location) by amplitude and theta

Cells: median [10th, 90th percentile] status. Rare ranks: number of ranks 4 to 10 whose record occupancy is inside its band.

| F_W | theta | F rank 1 (record 0.165) | F rank 2 (record 0.117) | F rank 3 (record 0.199) | rare ranks inside | sites, rank 4 (record 22) | sites, rank 7 (record 15) | sites, rank 10 (record 8) | richness (record 10) | beyond rank 10 | types originated | new among top 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0 | 0.027 [0.016, 0.051] above | 0.025 [0.015, 0.045] above | 0.026 [0.013, 0.049] above | 2 | 18 [7, 26] inside | 1 [0, 6] above | 0 [0, 0] above | 7 [6, 9] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.1 | 0.3 | 0.024 [0.011, 0.049] above | 0.026 [0.014, 0.045] above | 0.024 [0.011, 0.045] above | 2 | 22 [7, 27] inside | 3 [1, 7] above | 1 [0, 3] above | 10 [8, 14] inside | 0.000 [0.000, 0.001] inside | 52 | 2 |
| 0.1 | 1 | 0.025 [0.014, 0.039] above | 0.027 [0.011, 0.051] above | 0.022 [0.011, 0.052] above | 4 | 21 [10, 27] inside | 5 [2, 12] above | 3 [1, 6] above | 14 [11, 20] below | 0.002 [0.000, 0.006] below | 168 | 4 |
| 0.1 | 3 | 0.030 [0.016, 0.053] above | 0.024 [0.012, 0.052] above | 0.025 [0.017, 0.055] above | 7 | 20 [11, 27] inside | 11 [5, 17] inside | 7 [3, 10] inside | 31 [27, 37] below | 0.017 [0.009, 0.033] below | 510 | 5 |
| 0.1 | 10 | 0.032 [0.016, 0.065] above | 0.034 [0.017, 0.067] above | 0.026 [0.013, 0.062] above | 7 | 20 [13, 27] inside | 14 [9, 22] inside | 11 [7, 18] inside | 89 [77, 100] below | 0.108 [0.067, 0.146] below | 1703 | 6 |
| 0.3 | 0 | 0.036 [0.019, 0.130] above | 0.036 [0.014, 0.083] above | 0.041 [0.021, 0.075] above | 1 | 10 [3, 27] inside | 0 [0, 3] above | 0 [0, 0] above | 6 [5, 8] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.3 | 0.3 | 0.047 [0.020, 0.159] above | 0.044 [0.019, 0.133] inside | 0.041 [0.016, 0.085] above | 1 | 11 [5, 26] inside | 3 [1, 5] above | 1 [0, 2] above | 11 [9, 14] inside | 0.000 [0.000, 0.001] inside | 50 | 4 |
| 0.3 | 1 | 0.036 [0.016, 0.080] above | 0.046 [0.021, 0.134] inside | 0.044 [0.019, 0.093] above | 5 | 18 [8, 26] inside | 8 [2, 14] above | 4 [2, 8] inside | 23 [18, 28] below | 0.006 [0.002, 0.015] below | 174 | 5 |
| 0.3 | 3 | 0.051 [0.019, 0.101] above | 0.050 [0.019, 0.134] inside | 0.050 [0.027, 0.112] above | 7 | 18 [11, 27] inside | 11 [6, 19] inside | 8 [3, 15] inside | 52 [45, 64] below | 0.048 [0.029, 0.089] below | 509 | 6 |
| 0.3 | 10 | 0.079 [0.028, 0.177] inside | 0.094 [0.035, 0.197] inside | 0.052 [0.024, 0.107] above | 5 | 20 [15, 25] inside | 16 [11, 22] inside | 13 [9, 19] below | 159 [141, 168] below | 0.268 [0.209, 0.319] below | 1698 | 7 |
| 0.5 | 0 | 0.083 [0.029, 0.193] inside | 0.069 [0.034, 0.174] inside | 0.061 [0.025, 0.143] above | 0 | 8 [2, 17] above | 1 [0, 3] above | 0 [0, 0] above | 7 [5, 8] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.5 | 0.3 | 0.055 [0.025, 0.129] above | 0.060 [0.029, 0.179] inside | 0.061 [0.028, 0.176] above | 2 | 11 [4, 24] inside | 3 [1, 7] above | 1 [1, 3] above | 12 [10, 16] inside | 0.000 [0.000, 0.002] inside | 50 | 4 |
| 0.5 | 1 | 0.050 [0.026, 0.164] above | 0.098 [0.032, 0.173] inside | 0.064 [0.031, 0.232] inside | 6 | 16 [8, 27] inside | 8 [3, 14] above | 4 [2, 9] inside | 30 [22, 35] below | 0.013 [0.004, 0.029] below | 170 | 6 |
| 0.5 | 3 | 0.058 [0.020, 0.165] above | 0.072 [0.033, 0.167] inside | 0.084 [0.035, 0.151] above | 7 | 18 [12, 26] inside | 13 [8, 20] inside | 10 [5, 16] inside | 68 [58, 78] below | 0.090 [0.056, 0.112] below | 500 | 7 |
| 0.5 | 10 | 0.135 [0.053, 0.273] inside | 0.113 [0.037, 0.225] inside | 0.064 [0.033, 0.165] above | 5 | 17 [12, 22] inside | 16 [10, 22] inside | 12 [9, 18] below | 212 [193, 236] below | 0.364 [0.306, 0.445] below | 1696 | 7 |
| 0.7 | 0 | 0.127 [0.040, 0.297] inside | 0.094 [0.030, 0.264] inside | 0.056 [0.025, 0.246] inside | 0 | 5 [2, 12] above | 0 [0, 2] above | 0 [0, 0] above | 6 [5, 8] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.7 | 0.3 | 0.093 [0.036, 0.222] inside | 0.095 [0.037, 0.226] inside | 0.073 [0.032, 0.177] above | 1 | 10 [5, 20] above | 5 [1, 9] above | 2 [1, 5] above | 14 [10, 19] inside | 0.001 [0.000, 0.005] inside | 50 | 5 |
| 0.7 | 1 | 0.069 [0.024, 0.216] inside | 0.111 [0.053, 0.216] inside | 0.095 [0.039, 0.178] above | 6 | 14 [8, 22] inside | 9 [3, 14] above | 5 [3, 10] inside | 33 [28, 40] below | 0.020 [0.008, 0.040] below | 170 | 6 |
| 0.7 | 3 | 0.121 [0.035, 0.292] inside | 0.117 [0.054, 0.272] inside | 0.105 [0.062, 0.281] inside | 6 | 18 [10, 23] inside | 14 [8, 18] inside | 11 [7, 15] inside | 88 [75, 98] below | 0.141 [0.091, 0.204] below | 510 | 7 |
| 0.7 | 10 | 0.207 [0.096, 0.288] inside | 0.113 [0.064, 0.234] inside | 0.089 [0.047, 0.179] above | 6 | 17 [13, 22] inside | 14 [9, 19] inside | 12 [8, 17] inside | 262 [242, 290] below | 0.484 [0.402, 0.566] below | 1706 | 8 |
| 0.9 | 0 | 0.206 [0.043, 0.374] inside | 0.163 [0.063, 0.316] inside | 0.057 [0.021, 0.134] above | 0 | 6 [2, 15] above | 0 [0, 2] above | 0 [0, 0] above | 6 [4, 7] above | 0.000 [0.000, 0.000] inside | 0 | 0 |
| 0.9 | 0.3 | 0.117 [0.042, 0.300] inside | 0.151 [0.039, 0.317] inside | 0.100 [0.042, 0.230] inside | 2 | 11 [5, 19] above | 4 [2, 7] above | 2 [1, 5] above | 17 [14, 21] below | 0.002 [0.001, 0.008] below | 50 | 5 |
| 0.9 | 1 | 0.114 [0.044, 0.336] inside | 0.170 [0.073, 0.350] inside | 0.128 [0.043, 0.305] inside | 6 | 16 [7, 21] above | 10 [5, 17] inside | 7 [3, 11] inside | 42 [36, 49] below | 0.037 [0.020, 0.059] below | 171 | 7 |
| 0.9 | 3 | 0.197 [0.061, 0.337] inside | 0.171 [0.085, 0.289] inside | 0.138 [0.088, 0.243] inside | 5 | 16 [9, 21] above | 13 [8, 19] inside | 11 [7, 16] inside | 116 [105, 133] below | 0.210 [0.156, 0.306] below | 512 | 8 |
| 0.9 | 10 | 0.196 [0.071, 0.313] inside | 0.124 [0.072, 0.221] inside | 0.085 [0.049, 0.169] above | 6 | 14 [10, 20] above | 13 [8, 17] inside | 12 [8, 15] inside | 339 [313, 365] below | 0.607 [0.527, 0.655] below | 1699 | 8 |

### Calibration: theta = 0 beside generator C3

The comparison of C4 at theta = 0 with C3 on identical ranked readouts, with medians and bands for every statistic, is in `output/findings/type_origination_calibration.md` (`analyses/per_type/type_origination_calibration.py`).
