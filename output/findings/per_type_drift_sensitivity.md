# Regional-drift sensitivity: coverage and separation readings

Produced by `analyses/per_type/drift_sensitivity_figure.py` from the C3 and M per-record CSVs of `per_class_table.py` (regional-drift spec draft 3, section 5). Statuses: inside, outside or unavailable (band ineligible: fewer than 45 finite values). Separation: planted records with a defined statistic outside the neutral band of the same class and amplitude, with the undefined count beside them; the neutral cell's own count outside its band is a description of the band, not a calibration. Descriptive summaries only; no threshold other than the band, no test statistic, no success level; nothing here identifies the cause of the record's pattern.

Planting one class rescales the other classes' shares in the marked rows and redraws the sherds, so the realized occupancies of unplanted rare classes in a marker cell differ from the paired neutral record's; only the planted class is read. The paired direction (up / down / tied) is marked uninformative where the class was lost from the target in at least half of the neutral records.

## primary operating point

W = 8 L_R = 2336, S = 2600, L_d = 182, L_R = 292, burn-in = 5000, records per cell 50.

### Coverage status, record against the neutral bands

| class | record | F_W 0 | F_W 0.01 | F_W 0.03 | F_W 0.1 | F_W 0.3 |
|---|---|---|---|---|---|---|
| Parkin_Punctated | 0.165 | 0.013 [0.006, 0.025] n 50: **outside** | 0.016 [0.008, 0.029] n 50: **outside** | 0.016 [0.008, 0.027] n 50: **outside** | 0.025 [0.011, 0.051] n 50: **outside** | 0.035 [0.011, 0.077] n 50: **outside**; lost 2 |
| Barton/Kent/MPI | 0.117 | 0.013 [0.006, 0.025] n 50: **outside** | 0.016 [0.007, 0.025] n 50: **outside** | 0.016 [0.008, 0.026] n 50: **outside** | 0.024 [0.012, 0.045] n 50: **outside** | 0.032 [0.016, 0.111] n 50: **outside**; lost 6 |
| Painted | 0.199 | 0.014 [0.007, 0.021] n 50: **outside** | 0.016 [0.009, 0.024] n 50: **outside** | 0.014 [0.007, 0.024] n 50: **outside** | 0.023 [0.011, 0.041] n 50: **outside**; lost 8 | 0.029 [0.014, 0.052] n 50: **outside**; lost 32 |
| Fortune_Noded | 15 | 19 [16, 23] n 50: **outside** | 20 [10, 25] n 50: **inside**; lost 8 | 12 [2, 25] n 50: **inside**; lost 34 | 9 [2, 26] n 50: **inside**; lost 41; absent 4 | 4 [0, 21] n 50: **inside**; lost 46; absent 8 |
| Ranch_Incised | 22 | 23 [19, 25] n 50: **inside** | 22 [13, 26] n 50: **inside**; lost 1 | 20 [7, 26] n 50: **inside**; lost 20 | 8 [2, 27] n 50: **inside**; lost 39; absent 1 | 6 [1, 24] n 50: **inside**; lost 46; absent 4 |
| Walls_Engraved | 14 | 18 [14, 21] n 50: **inside** | 16 [9, 22] n 50: **inside**; lost 5 | 14 [3, 24] n 50: **inside**; lost 29 | 8 [0, 25] n 50: **inside**; lost 37; absent 8 | 3 [0, 28] n 50: **inside**; lost 45; absent 11 |
| Wallace_Incised | 7 | 8 [5, 14] n 50: **inside**; absent 1 | 7 [1, 16] n 50: **inside**; lost 29; absent 5 | 5 [0, 21] n 50: **inside**; lost 38; absent 10 | 1 [0, 7] n 50: **inside**; lost 50; absent 22 | 0 [0, 6] n 50: **outside**; lost 50; absent 30 |
| Rhodes_Incised | 15 | 11 [7, 18] n 50: **inside** | 10 [1, 18] n 50: **inside**; lost 22; absent 2 | 4 [0, 23] n 50: **inside**; lost 38; absent 8 | 1 [0, 19] n 50: **inside**; lost 46; absent 23 | 0 [0, 13] n 50: **outside**; lost 49; absent 26 |
| Vernon_Paul_Applique | 10 | 7 [3, 12] n 50: **inside** | 4 [0, 14] n 50: **inside**; lost 32; absent 8 | 2 [0, 16] n 50: **inside**; lost 41; absent 15 | 0 [0, 11] n 50: **inside**; lost 48; absent 27 | 0 [0, 3] n 50: **outside**; lost 50; absent 34 |
| Hull_Engraved | 8 | 4 [1, 6] n 50: **outside**; absent 2 | 0 [0, 11] n 50: **inside**; lost 44; absent 25 | 0 [0, 8] n 50: **inside**; lost 48; absent 33 | 0 [0, 3] n 50: **outside**; lost 50; absent 39 | 0 [0, 3] n 50: **outside**; lost 50; absent 37 |

Sets over the tested grid: C = {}; R = {0.01, 0.03}; H = {}; E_C = {0, 0.01, 0.03, 0.1, 0.3}; E_R = {0, 0.01, 0.03, 0.1, 0.3}; E_H = {0, 0.01, 0.03, 0.1, 0.3}.

Reading: not inside the bands at any evaluable tested amplitude (common classes evaluable at {0, 0.01, 0.03, 0.1, 0.3}, inside at {}; rare classes evaluable at {0, 0.01, 0.03, 0.1, 0.3}, inside at {0.01, 0.03}).

### Separation of planted markers

### Barton/Kent/MPI, large cluster (17 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 25 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.025 (eligible) | [0.006, 0.025] n 50 | 47 / 3 / 0 / 0 |
| 0 | 2 | 45 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.044 (eligible) | [0.006, 0.025] n 50 | 49 / 1 / 0 / 0 |
| 0 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.085 (eligible) | [0.006, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0.01 | 1.5 | 26 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.026 (eligible) | [0.007, 0.025] n 50 | 43 / 7 / 0 / 0 |
| 0.01 | 2 | 43 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.045 (eligible) | [0.007, 0.025] n 50 | 48 / 2 / 0 / 0 |
| 0.01 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.087 (eligible) | [0.007, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0.03 | 1.5 | 23 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.024 (eligible) | [0.008, 0.026] n 50 | 43 / 7 / 0 / 0 |
| 0.03 | 2 | 43 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.043 (eligible) | [0.008, 0.026] n 50 | 48 / 2 / 0 / 0 |
| 0.03 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.087 (eligible) | [0.008, 0.026] n 50 | 50 / 0 / 0 / 0 |
| 0.1 | 1.5 | 18 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.032 (eligible) | [0.012, 0.045] n 50 | 44 / 6 / 0 / 0 |
| 0.1 | 2 | 30 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.052 (eligible) | [0.012, 0.045] n 50 | 49 / 1 / 0 / 0 |
| 0.1 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.089 (eligible) | [0.012, 0.045] n 50 | 50 / 0 / 0 / 0 |
| 0.3 | 1.5 | 9 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.044 (eligible) | [0.016, 0.111] n 50 | 40 / 10 / 0 / 0 |
| 0.3 | 2 | 9 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.060 (eligible) | [0.016, 0.111] n 50 | 45 / 5 / 0 / 0 |
| 0.3 | 3 | 18 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.094 (eligible) | [0.016, 0.111] n 50 | 46 / 4 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: 1.5; F_W 0.01: 1.5; F_W 0.03: 2; F_W 0.1: 2; F_W 0.3: none on this grid.

### Barton/Kent/MPI, small cluster (11 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 25 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.023 (eligible) | [0.006, 0.025] n 50 | 46 / 4 / 0 / 0 |
| 0 | 2 | 44 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.040 (eligible) | [0.006, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.080 (eligible) | [0.006, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0.01 | 1.5 | 22 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.024 (eligible) | [0.007, 0.025] n 50 | 37 / 13 / 0 / 0 |
| 0.01 | 2 | 39 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.042 (eligible) | [0.007, 0.025] n 50 | 48 / 2 / 0 / 0 |
| 0.01 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.083 (eligible) | [0.007, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0.03 | 1.5 | 24 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.023 (eligible) | [0.008, 0.026] n 50 | 44 / 6 / 0 / 0 |
| 0.03 | 2 | 42 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.042 (eligible) | [0.008, 0.026] n 50 | 49 / 1 / 0 / 0 |
| 0.03 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.078 (eligible) | [0.008, 0.026] n 50 | 50 / 0 / 0 / 0 |
| 0.1 | 1.5 | 8 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.033 (eligible) | [0.012, 0.045] n 50 | 41 / 9 / 0 / 0 |
| 0.1 | 2 | 32 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.049 (eligible) | [0.012, 0.045] n 50 | 46 / 4 / 0 / 0 |
| 0.1 | 3 | 49 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.089 (eligible) | [0.012, 0.045] n 50 | 49 / 1 / 0 / 0 |
| 0.3 | 1.5 | 6 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.041 (eligible) | [0.016, 0.111] n 50 | 43 / 7 / 0 / 0 |
| 0.3 | 2 | 10 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.063 (eligible) | [0.016, 0.111] n 50 | 47 / 3 / 0 / 0 |
| 0.3 | 3 | 24 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.103 (eligible) | [0.016, 0.111] n 50 | 49 / 1 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: 2; F_W 0.01: 2; F_W 0.03: 2; F_W 0.1: 2; F_W 0.3: none on this grid.

### Painted, large cluster (17 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 20 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.017 (eligible) | [0.007, 0.021] n 50 | 37 / 13 / 0 / 0 |
| 0 | 2 | 32 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.026 (eligible) | [0.007, 0.021] n 50 | 45 / 5 / 0 / 0 |
| 0 | 3 | 47 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.051 (eligible) | [0.007, 0.021] n 50 | 48 / 2 / 0 / 0 |
| 0.01 | 1.5 | 24 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.021 (eligible) | [0.009, 0.024] n 50 | 40 / 10 / 0 / 0 |
| 0.01 | 2 | 34 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.031 (eligible) | [0.009, 0.024] n 50 | 46 / 4 / 0 / 0 |
| 0.01 | 3 | 49 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.058 (eligible) | [0.009, 0.024] n 50 | 50 / 0 / 0 / 0 |
| 0.03 | 1.5 | 25 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.024 (eligible) | [0.007, 0.024] n 50 | 42 / 8 / 0 / 0 |
| 0.03 | 2 | 41 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.034 (eligible) | [0.007, 0.024] n 50 | 47 / 3 / 0 / 0 |
| 0.03 | 3 | 49 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.057 (eligible) | [0.007, 0.024] n 50 | 50 / 0 / 0 / 0 |
| 0.1 | 1.5 | 21 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.032 (eligible) | [0.011, 0.041] n 50 | 44 / 6 / 0 / 0 |
| 0.1 | 2 | 32 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.047 (eligible) | [0.011, 0.041] n 50 | 46 / 4 / 0 / 0 |
| 0.1 | 3 | 45 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.076 (eligible) | [0.011, 0.041] n 50 | 50 / 0 / 0 / 0 |
| 0.3 | 1.5 | 21 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.045 (eligible) | [0.014, 0.052] n 50 | 47 / 3 / 0 / 0 (class lost from the target in 32 of 50 neutral records; paired direction uninformative) |
| 0.3 | 2 | 28 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.055 (eligible) | [0.014, 0.052] n 50 | 48 / 2 / 0 / 0 (class lost from the target in 32 of 50 neutral records; paired direction uninformative) |
| 0.3 | 3 | 42 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.092 (eligible) | [0.014, 0.052] n 50 | 49 / 1 / 0 / 0 (class lost from the target in 32 of 50 neutral records; paired direction uninformative) |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: 2; F_W 0.01: 2; F_W 0.03: 2; F_W 0.1: 2; F_W 0.3: 2.

### Painted, small cluster (11 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 25 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.020 (eligible) | [0.007, 0.021] n 50 | 41 / 9 / 0 / 0 |
| 0 | 2 | 35 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.031 (eligible) | [0.007, 0.021] n 50 | 44 / 6 / 0 / 0 |
| 0 | 3 | 46 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.059 (eligible) | [0.007, 0.021] n 50 | 47 / 3 / 0 / 0 |
| 0.01 | 1.5 | 25 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.021 (eligible) | [0.009, 0.024] n 50 | 39 / 11 / 0 / 0 |
| 0.01 | 2 | 34 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.034 (eligible) | [0.009, 0.024] n 50 | 46 / 4 / 0 / 0 |
| 0.01 | 3 | 46 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.061 (eligible) | [0.009, 0.024] n 50 | 49 / 1 / 0 / 0 |
| 0.03 | 1.5 | 19 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.020 (eligible) | [0.007, 0.024] n 50 | 41 / 9 / 0 / 0 |
| 0.03 | 2 | 34 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.030 (eligible) | [0.007, 0.024] n 50 | 45 / 5 / 0 / 0 |
| 0.03 | 3 | 47 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.055 (eligible) | [0.007, 0.024] n 50 | 49 / 1 / 0 / 0 |
| 0.1 | 1.5 | 22 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.030 (eligible) | [0.011, 0.041] n 50 | 39 / 11 / 0 / 0 |
| 0.1 | 2 | 28 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.042 (eligible) | [0.011, 0.041] n 50 | 42 / 8 / 0 / 0 |
| 0.1 | 3 | 37 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.065 (eligible) | [0.011, 0.041] n 50 | 46 / 4 / 0 / 0 |
| 0.3 | 1.5 | 13 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.035 (eligible) | [0.014, 0.052] n 50 | 41 / 8 / 1 / 0 (class lost from the target in 32 of 50 neutral records; paired direction uninformative) |
| 0.3 | 2 | 22 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.048 (eligible) | [0.014, 0.052] n 50 | 44 / 5 / 1 / 0 (class lost from the target in 32 of 50 neutral records; paired direction uninformative) |
| 0.3 | 3 | 38 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.075 (eligible) | [0.014, 0.052] n 50 | 46 / 3 / 1 / 0 (class lost from the target in 32 of 50 neutral records; paired direction uninformative) |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: 2; F_W 0.01: 2; F_W 0.03: 2; F_W 0.1: 2; F_W 0.3: 3.

### Parkin_Punctated, large cluster (17 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 21 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.022 (eligible) | [0.006, 0.025] n 50 | 43 / 7 / 0 / 0 |
| 0 | 2 | 43 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.040 (eligible) | [0.006, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.071 (eligible) | [0.006, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0.01 | 1.5 | 21 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.027 (eligible) | [0.008, 0.029] n 50 | 37 / 13 / 0 / 0 |
| 0.01 | 2 | 37 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.042 (eligible) | [0.008, 0.029] n 50 | 47 / 3 / 0 / 0 |
| 0.01 | 3 | 49 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.076 (eligible) | [0.008, 0.029] n 50 | 50 / 0 / 0 / 0 |
| 0.03 | 1.5 | 23 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.025 (eligible) | [0.008, 0.027] n 50 | 39 / 11 / 0 / 0 |
| 0.03 | 2 | 33 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.038 (eligible) | [0.008, 0.027] n 50 | 48 / 2 / 0 / 0 |
| 0.03 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.066 (eligible) | [0.008, 0.027] n 50 | 50 / 0 / 0 / 0 |
| 0.1 | 1.5 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.033 (eligible) | [0.011, 0.051] n 50 | 34 / 16 / 0 / 0 |
| 0.1 | 2 | 23 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.047 (eligible) | [0.011, 0.051] n 50 | 41 / 9 / 0 / 0 |
| 0.1 | 3 | 42 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.077 (eligible) | [0.011, 0.051] n 50 | 48 / 2 / 0 / 0 |
| 0.3 | 1.5 | 7 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.037 (eligible) | [0.011, 0.077] n 50 | 37 / 13 / 0 / 0 |
| 0.3 | 2 | 14 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.051 (eligible) | [0.011, 0.077] n 50 | 42 / 8 / 0 / 0 |
| 0.3 | 3 | 30 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.083 (eligible) | [0.011, 0.077] n 50 | 46 / 4 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: 2; F_W 0.01: 2; F_W 0.03: 2; F_W 0.1: 3; F_W 0.3: 3.

### Parkin_Punctated, small cluster (11 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 21 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.022 (eligible) | [0.006, 0.025] n 50 | 43 / 7 / 0 / 0 |
| 0 | 2 | 41 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.039 (eligible) | [0.006, 0.025] n 50 | 47 / 3 / 0 / 0 |
| 0 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.074 (eligible) | [0.006, 0.025] n 50 | 50 / 0 / 0 / 0 |
| 0.01 | 1.5 | 21 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.025 (eligible) | [0.008, 0.029] n 50 | 39 / 11 / 0 / 0 |
| 0.01 | 2 | 35 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.040 (eligible) | [0.008, 0.029] n 50 | 48 / 2 / 0 / 0 |
| 0.01 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.073 (eligible) | [0.008, 0.029] n 50 | 49 / 1 / 0 / 0 |
| 0.03 | 1.5 | 21 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.024 (eligible) | [0.008, 0.027] n 50 | 38 / 12 / 0 / 0 |
| 0.03 | 2 | 41 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.040 (eligible) | [0.008, 0.027] n 50 | 44 / 6 / 0 / 0 |
| 0.03 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.075 (eligible) | [0.008, 0.027] n 50 | 49 / 1 / 0 / 0 |
| 0.1 | 1.5 | 12 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.036 (eligible) | [0.011, 0.051] n 50 | 47 / 3 / 0 / 0 |
| 0.1 | 2 | 23 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.050 (eligible) | [0.011, 0.051] n 50 | 49 / 1 / 0 / 0 |
| 0.1 | 3 | 46 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.079 (eligible) | [0.011, 0.051] n 50 | 50 / 0 / 0 / 0 |
| 0.3 | 1.5 | 5 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.039 (eligible) | [0.011, 0.077] n 50 | 38 / 12 / 0 / 0 |
| 0.3 | 2 | 13 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.050 (eligible) | [0.011, 0.077] n 50 | 42 / 8 / 0 / 0 |
| 0.3 | 3 | 25 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.079 (eligible) | [0.011, 0.077] n 50 | 45 / 5 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: 2; F_W 0.01: 2; F_W 0.03: 2; F_W 0.1: 3; F_W 0.3: 3.

### Walls_Engraved, large cluster (17 sites), statistic: sites present

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 16 of 50 (0 undefined) | 10 of 50 (0 undefined) | 19 (eligible) | [14, 21] n 50 | 31 / 5 / 14 / 0 |
| 0 | 2 | 19 of 50 (0 undefined) | 10 of 50 (0 undefined) | 20 (eligible) | [14, 21] n 50 | 37 / 1 / 12 / 0 |
| 0 | 3 | 22 of 50 (0 undefined) | 10 of 50 (0 undefined) | 21 (eligible) | [14, 21] n 50 | 42 / 0 / 8 / 0 |
| 0.01 | 1.5 | 13 of 50 (0 undefined) | 10 of 50 (0 undefined) | 17 (eligible) | [9, 22] n 50 | 34 / 1 / 15 / 0 |
| 0.01 | 2 | 16 of 50 (0 undefined) | 10 of 50 (0 undefined) | 19 (eligible) | [9, 22] n 50 | 38 / 0 / 12 / 0 |
| 0.01 | 3 | 16 of 50 (0 undefined) | 10 of 50 (0 undefined) | 19 (eligible) | [9, 22] n 50 | 45 / 0 / 5 / 0 |
| 0.03 | 1.5 | 11 of 50 (0 undefined) | 7 of 50 (0 undefined) | 14 (eligible) | [3, 24] n 50 | 23 / 4 / 23 / 0 (class lost from the target in 29 of 50 neutral records; paired direction uninformative) |
| 0.03 | 2 | 11 of 50 (0 undefined) | 7 of 50 (0 undefined) | 14 (eligible) | [3, 24] n 50 | 24 / 2 / 24 / 0 (class lost from the target in 29 of 50 neutral records; paired direction uninformative) |
| 0.03 | 3 | 13 of 50 (0 undefined) | 7 of 50 (0 undefined) | 15 (eligible) | [3, 24] n 50 | 34 / 1 / 15 / 0 (class lost from the target in 29 of 50 neutral records; paired direction uninformative) |
| 0.1 | 1.5 | 6 of 50 (0 undefined) | 4 of 50 (0 undefined) | 8 (eligible) | [0, 25] n 50 | 16 / 1 / 33 / 0 (class lost from the target in 37 of 50 neutral records; paired direction uninformative) |
| 0.1 | 2 | 6 of 50 (0 undefined) | 4 of 50 (0 undefined) | 9 (eligible) | [0, 25] n 50 | 19 / 0 / 31 / 0 (class lost from the target in 37 of 50 neutral records; paired direction uninformative) |
| 0.1 | 3 | 8 of 50 (0 undefined) | 4 of 50 (0 undefined) | 9 (eligible) | [0, 25] n 50 | 25 / 0 / 25 / 0 (class lost from the target in 37 of 50 neutral records; paired direction uninformative) |
| 0.3 | 1.5 | 0 of 50 (0 undefined) | 0 of 50 (0 undefined) | 4 (eligible) | [0, 28] n 50 | 10 / 0 / 40 / 0 (class lost from the target in 45 of 50 neutral records; paired direction uninformative) |
| 0.3 | 2 | 0 of 50 (0 undefined) | 0 of 50 (0 undefined) | 4 (eligible) | [0, 28] n 50 | 16 / 0 / 34 / 0 (class lost from the target in 45 of 50 neutral records; paired direction uninformative) |
| 0.3 | 3 | 0 of 50 (0 undefined) | 0 of 50 (0 undefined) | 4 (eligible) | [0, 28] n 50 | 19 / 1 / 30 / 0 (class lost from the target in 45 of 50 neutral records; paired direction uninformative) |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: none on this grid; F_W 0.01: none on this grid; F_W 0.03: none on this grid; F_W 0.1: none on this grid; F_W 0.3: none on this grid.

### Walls_Engraved, small cluster (11 sites), statistic: sites present

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0 | 1.5 | 12 of 50 (0 undefined) | 10 of 50 (0 undefined) | 19 (eligible) | [14, 21] n 50 | 27 / 3 / 20 / 0 |
| 0 | 2 | 14 of 50 (0 undefined) | 10 of 50 (0 undefined) | 20 (eligible) | [14, 21] n 50 | 35 / 1 / 14 / 0 |
| 0 | 3 | 15 of 50 (0 undefined) | 10 of 50 (0 undefined) | 20 (eligible) | [14, 21] n 50 | 41 / 0 / 9 / 0 |
| 0.01 | 1.5 | 12 of 50 (0 undefined) | 10 of 50 (0 undefined) | 17 (eligible) | [9, 22] n 50 | 19 / 4 / 27 / 0 |
| 0.01 | 2 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 17 (eligible) | [9, 22] n 50 | 26 / 2 / 22 / 0 |
| 0.01 | 3 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 18 (eligible) | [9, 22] n 50 | 29 / 1 / 20 / 0 |
| 0.03 | 1.5 | 8 of 50 (0 undefined) | 7 of 50 (0 undefined) | 14 (eligible) | [3, 24] n 50 | 18 / 2 / 30 / 0 (class lost from the target in 29 of 50 neutral records; paired direction uninformative) |
| 0.03 | 2 | 9 of 50 (0 undefined) | 7 of 50 (0 undefined) | 15 (eligible) | [3, 24] n 50 | 26 / 1 / 23 / 0 (class lost from the target in 29 of 50 neutral records; paired direction uninformative) |
| 0.03 | 3 | 9 of 50 (0 undefined) | 7 of 50 (0 undefined) | 16 (eligible) | [3, 24] n 50 | 31 / 0 / 19 / 0 (class lost from the target in 29 of 50 neutral records; paired direction uninformative) |
| 0.1 | 1.5 | 6 of 50 (0 undefined) | 4 of 50 (0 undefined) | 9 (eligible) | [0, 25] n 50 | 12 / 0 / 38 / 0 (class lost from the target in 37 of 50 neutral records; paired direction uninformative) |
| 0.1 | 2 | 6 of 50 (0 undefined) | 4 of 50 (0 undefined) | 9 (eligible) | [0, 25] n 50 | 16 / 0 / 34 / 0 (class lost from the target in 37 of 50 neutral records; paired direction uninformative) |
| 0.1 | 3 | 7 of 50 (0 undefined) | 4 of 50 (0 undefined) | 9 (eligible) | [0, 25] n 50 | 16 / 0 / 34 / 0 (class lost from the target in 37 of 50 neutral records; paired direction uninformative) |
| 0.3 | 1.5 | 0 of 50 (0 undefined) | 0 of 50 (0 undefined) | 3 (eligible) | [0, 28] n 50 | 3 / 1 / 46 / 0 (class lost from the target in 45 of 50 neutral records; paired direction uninformative) |
| 0.3 | 2 | 0 of 50 (0 undefined) | 0 of 50 (0 undefined) | 3 (eligible) | [0, 28] n 50 | 5 / 0 / 45 / 0 (class lost from the target in 45 of 50 neutral records; paired direction uninformative) |
| 0.3 | 3 | 0 of 50 (0 undefined) | 0 of 50 (0 undefined) | 3 (eligible) | [0, 28] n 50 | 7 / 0 / 43 / 0 (class lost from the target in 45 of 50 neutral records; paired direction uninformative) |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: none on this grid; F_W 0.01: none on this grid; F_W 0.03: none on this grid; F_W 0.1: none on this grid; F_W 0.3: none on this grid.

## secondary operating point

W = 8 L_R = 320, S = 340, L_d = 12, L_R = 40, burn-in = 5000, records per cell 50.

Marker cells at this point were run at F_W 0.1 only, the fallback of spec 3.4: H was empty on the primary grid.

### Coverage status, record against the neutral bands

| class | record | F_W 0 | F_W 0.01 | F_W 0.03 | F_W 0.1 | F_W 0.3 |
|---|---|---|---|---|---|---|
| Parkin_Punctated | 0.165 | 0.017 [0.010, 0.025] n 50: **outside** | 0.017 [0.011, 0.031] n 50: **outside** | 0.020 [0.010, 0.042] n 50: **outside** | 0.025 [0.015, 0.048] n 50: **outside** | 0.033 [0.014, 0.092] n 50: **outside** |
| Barton/Kent/MPI | 0.117 | 0.017 [0.010, 0.026] n 50: **outside** | 0.017 [0.011, 0.025] n 50: **outside** | 0.020 [0.010, 0.035] n 50: **outside** | 0.026 [0.014, 0.049] n 50: **outside** | 0.035 [0.013, 0.069] n 50: **outside**; lost 3 |
| Painted | 0.199 | 0.017 [0.009, 0.025] n 50: **outside** | 0.018 [0.009, 0.032] n 50: **outside** | 0.016 [0.009, 0.031] n 50: **outside** | 0.022 [0.014, 0.041] n 50: **outside**; lost 6 | 0.041 [0.018, 0.082] n 50: **outside**; lost 26 |
| Fortune_Noded | 15 | 18 [13, 22] n 50: **inside** | 14 [6, 23] n 50: **inside**; lost 7 | 12 [2, 24] n 50: **inside**; lost 27; absent 2 | 8 [1, 27] n 50: **inside**; lost 37; absent 5 | 2 [0, 11] n 50: **outside**; lost 50; absent 14 |
| Ranch_Incised | 22 | 20 [17, 24] n 50: **inside** | 20 [11, 24] n 50: **inside**; lost 1 | 18 [4, 25] n 50: **inside**; lost 19 | 6 [2, 26] n 50: **inside**; lost 36; absent 3 | 5 [0, 27] n 50: **inside**; lost 42; absent 7 |
| Walls_Engraved | 14 | 16 [11, 20] n 50: **inside** | 16 [6, 23] n 50: **inside**; lost 5 | 9 [1, 23] n 50: **inside**; lost 33; absent 4 | 6 [1, 24] n 50: **inside**; lost 38; absent 5 | 3 [0, 18] n 50: **inside**; lost 49; absent 16 |
| Wallace_Incised | 7 | 6 [3, 13] n 50: **inside** | 3 [0, 12] n 50: **inside**; lost 30; absent 8 | 3 [0, 17] n 50: **inside**; lost 35; absent 12 | 1 [0, 14] n 50: **inside**; lost 47; absent 19 | 0 [0, 7] n 50: **inside**; lost 49; absent 29 |
| Rhodes_Incised | 15 | 10 [5, 14] n 50: **outside** | 8 [1, 17] n 50: **inside**; lost 23; absent 3 | 10 [0, 21] n 50: **inside**; lost 28; absent 7 | 1 [0, 14] n 50: **outside**; lost 46; absent 19 | 0 [0, 6] n 50: **outside**; lost 49; absent 27 |
| Vernon_Paul_Applique | 10 | 6 [3, 10] n 50: **inside** | 2 [0, 11] n 50: **inside**; lost 34; absent 14 | 1 [0, 17] n 50: **inside**; lost 42; absent 21 | 1 [0, 6] n 50: **outside**; lost 49; absent 24 | 0 [0, 3] n 50: **outside**; lost 49; absent 33 |
| Hull_Engraved | 8 | 3 [1, 5] n 50: **outside**; absent 4 | 1 [0, 7] n 50: **outside**; lost 40; absent 20 | 0 [0, 6] n 50: **outside**; lost 45; absent 27 | 0 [0, 4] n 50: **outside**; lost 49; absent 39 | 0 [0, 4] n 50: **outside**; lost 50; absent 35 |

Sets over the tested grid: C = {}; R = {}; H = {}; E_C = {0, 0.01, 0.03, 0.1, 0.3}; E_R = {0, 0.01, 0.03, 0.1, 0.3}; E_H = {0, 0.01, 0.03, 0.1, 0.3}.

Reading: not inside the bands at any evaluable tested amplitude (common classes evaluable at {0, 0.01, 0.03, 0.1, 0.3}, inside at {}; rare classes evaluable at {0, 0.01, 0.03, 0.1, 0.3}, inside at {}).

### Separation of planted markers

### Barton/Kent/MPI, large cluster (17 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.036 (eligible) | [0.014, 0.049] n 50 | 45 / 5 / 0 / 0 |
| 0.1 | 2 | 26 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.050 (eligible) | [0.014, 0.049] n 50 | 49 / 1 / 0 / 0 |
| 0.1 | 3 | 50 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.092 (eligible) | [0.014, 0.049] n 50 | 50 / 0 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: 2; F_W 0.3: not run.

### Barton/Kent/MPI, small cluster (11 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 14 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.036 (eligible) | [0.014, 0.049] n 50 | 43 / 7 / 0 / 0 |
| 0.1 | 2 | 31 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.055 (eligible) | [0.014, 0.049] n 50 | 47 / 3 / 0 / 0 |
| 0.1 | 3 | 49 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.094 (eligible) | [0.014, 0.049] n 50 | 49 / 1 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: 2; F_W 0.3: not run.

### Painted, large cluster (17 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 18 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.030 (eligible) | [0.014, 0.041] n 50 | 45 / 5 / 0 / 0 |
| 0.1 | 2 | 30 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.042 (eligible) | [0.014, 0.041] n 50 | 49 / 1 / 0 / 0 |
| 0.1 | 3 | 46 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.071 (eligible) | [0.014, 0.041] n 50 | 50 / 0 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: 2; F_W 0.3: not run.

### Painted, small cluster (11 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 18 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.031 (eligible) | [0.014, 0.041] n 50 | 41 / 9 / 0 / 0 |
| 0.1 | 2 | 29 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.043 (eligible) | [0.014, 0.041] n 50 | 49 / 1 / 0 / 0 |
| 0.1 | 3 | 41 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.072 (eligible) | [0.014, 0.041] n 50 | 50 / 0 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: 2; F_W 0.3: not run.

### Parkin_Punctated, large cluster (17 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 12 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.033 (eligible) | [0.015, 0.048] n 50 | 43 / 7 / 0 / 0 |
| 0.1 | 2 | 26 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.049 (eligible) | [0.015, 0.048] n 50 | 47 / 3 / 0 / 0 |
| 0.1 | 3 | 49 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.077 (eligible) | [0.015, 0.048] n 50 | 50 / 0 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: 2; F_W 0.3: not run.

### Parkin_Punctated, small cluster (11 sites), statistic: F

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 7 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.033 (eligible) | [0.015, 0.048] n 50 | 42 / 8 / 0 / 0 |
| 0.1 | 2 | 26 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.048 (eligible) | [0.015, 0.048] n 50 | 48 / 2 / 0 / 0 |
| 0.1 | 3 | 49 of 50 (0 undefined) | 10 of 50 (0 undefined) | 0.076 (eligible) | [0.015, 0.048] n 50 | 50 / 0 / 0 / 0 |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: 2; F_W 0.3: not run.

### Walls_Engraved, large cluster (17 sites), statistic: sites present

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 7 (eligible) | [1, 24] n 50 | 11 / 3 / 36 / 0 (class lost from the target in 38 of 50 neutral records; paired direction uninformative) |
| 0.1 | 2 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 7 (eligible) | [1, 24] n 50 | 16 / 4 / 30 / 0 (class lost from the target in 38 of 50 neutral records; paired direction uninformative) |
| 0.1 | 3 | 10 of 50 (0 undefined) | 10 of 50 (0 undefined) | 8 (eligible) | [1, 24] n 50 | 23 / 0 / 27 / 0 (class lost from the target in 38 of 50 neutral records; paired direction uninformative) |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: none on this grid; F_W 0.3: not run.

### Walls_Engraved, small cluster (11 sites), statistic: sites present

| F_W | m | planted outside the neutral band | neutral outside its own band | planted median (eligible) | neutral band | paired: up / down / tied / unavailable |
|---|---|---|---|---|---|---|
| 0.1 | 1.5 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 7 (eligible) | [1, 24] n 50 | 8 / 0 / 42 / 0 (class lost from the target in 38 of 50 neutral records; paired direction uninformative) |
| 0.1 | 2 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 7 (eligible) | [1, 24] n 50 | 10 / 0 / 40 / 0 (class lost from the target in 38 of 50 neutral records; paired direction uninformative) |
| 0.1 | 3 | 11 of 50 (0 undefined) | 10 of 50 (0 undefined) | 7 (eligible) | [1, 24] n 50 | 13 / 0 / 37 / 0 (class lost from the target in 38 of 50 neutral records; paired direction uninformative) |

First tested m whose eligible planted median lies outside the neutral band: F_W 0: not run; F_W 0.01: not run; F_W 0.03: not run; F_W 0.1: none on this grid; F_W 0.3: not run.

## Figure

`figures/per_type_drift_sensitivity.png` (primary operating point; W = 8 L_R = 2336, S = 2600, 50 records per cell, record-key seed base 2026101950 + 100000 point + 1000 a + r). (a) Per-class moment F of the three common classes against F_W: neutral median and 10th to 90th percentile band (NumPy linear quantile rule, endpoints inclusive), the record's value as a dashed line, and the planted medians for each cluster (colors) and multiplier (solid 1.5, dashed 2, dotted 3). (b) Fraction of planted records with a defined statistic outside the neutral band of the same class and amplitude. (c) Rare-class sites present: neutral median and band per class, the record's value as a short horizontal bar, and the planted Walls_Engraved occupancy at m = 3 in each cluster.

## The sentence of spec section 1, filled in

> Under neutral copying with regional drift of the innovation target (Dirichlet recursion), at the primary operating point, window W = 8 L_R = 2336 steps, 50 records per cell (every common-class band over 50 finite values), the record's common-class F values and rare-class occupancies are not inside the bands at any evaluable tested amplitude (common classes evaluable at {0, 0.01, 0.03, 0.1, 0.3}, inside at {}; rare classes evaluable at {0, 0.01, 0.03, 0.1, 0.3}, inside at {0.01, 0.03}). For Barton/Kent/MPI multiplied on the odds scale in the large cluster (17 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 25 of 50 (0 undefined); m 2: 45 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 1.5); F_W 0.01 (m 1.5: 26 of 50 (0 undefined); m 2: 43 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 1.5); F_W 0.03 (m 1.5: 23 of 50 (0 undefined); m 2: 43 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.1 (m 1.5: 18 of 50 (0 undefined); m 2: 30 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.3 (m 1.5: 9 of 50 (0 undefined); m 2: 9 of 50 (0 undefined); m 3: 18 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid). For Barton/Kent/MPI multiplied on the odds scale in the small cluster (11 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 25 of 50 (0 undefined); m 2: 44 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.01 (m 1.5: 22 of 50 (0 undefined); m 2: 39 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.03 (m 1.5: 24 of 50 (0 undefined); m 2: 42 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.1 (m 1.5: 8 of 50 (0 undefined); m 2: 32 of 50 (0 undefined); m 3: 49 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.3 (m 1.5: 6 of 50 (0 undefined); m 2: 10 of 50 (0 undefined); m 3: 24 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid). For Painted multiplied on the odds scale in the large cluster (17 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 20 of 50 (0 undefined); m 2: 32 of 50 (0 undefined); m 3: 47 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.01 (m 1.5: 24 of 50 (0 undefined); m 2: 34 of 50 (0 undefined); m 3: 49 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.03 (m 1.5: 25 of 50 (0 undefined); m 2: 41 of 50 (0 undefined); m 3: 49 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.1 (m 1.5: 21 of 50 (0 undefined); m 2: 32 of 50 (0 undefined); m 3: 45 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.3 (m 1.5: 21 of 50 (0 undefined); m 2: 28 of 50 (0 undefined); m 3: 42 of 50 (0 undefined); first tested m with the eligible median outside: 2). For Painted multiplied on the odds scale in the small cluster (11 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 25 of 50 (0 undefined); m 2: 35 of 50 (0 undefined); m 3: 46 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.01 (m 1.5: 25 of 50 (0 undefined); m 2: 34 of 50 (0 undefined); m 3: 46 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.03 (m 1.5: 19 of 50 (0 undefined); m 2: 34 of 50 (0 undefined); m 3: 47 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.1 (m 1.5: 22 of 50 (0 undefined); m 2: 28 of 50 (0 undefined); m 3: 37 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.3 (m 1.5: 13 of 50 (0 undefined); m 2: 22 of 50 (0 undefined); m 3: 38 of 50 (0 undefined); first tested m with the eligible median outside: 3). For Parkin_Punctated multiplied on the odds scale in the large cluster (17 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 21 of 50 (0 undefined); m 2: 43 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.01 (m 1.5: 21 of 50 (0 undefined); m 2: 37 of 50 (0 undefined); m 3: 49 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.03 (m 1.5: 23 of 50 (0 undefined); m 2: 33 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.1 (m 1.5: 11 of 50 (0 undefined); m 2: 23 of 50 (0 undefined); m 3: 42 of 50 (0 undefined); first tested m with the eligible median outside: 3); F_W 0.3 (m 1.5: 7 of 50 (0 undefined); m 2: 14 of 50 (0 undefined); m 3: 30 of 50 (0 undefined); first tested m with the eligible median outside: 3). For Parkin_Punctated multiplied on the odds scale in the small cluster (11 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 21 of 50 (0 undefined); m 2: 41 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.01 (m 1.5: 21 of 50 (0 undefined); m 2: 35 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.03 (m 1.5: 21 of 50 (0 undefined); m 2: 41 of 50 (0 undefined); m 3: 50 of 50 (0 undefined); first tested m with the eligible median outside: 2); F_W 0.1 (m 1.5: 12 of 50 (0 undefined); m 2: 23 of 50 (0 undefined); m 3: 46 of 50 (0 undefined); first tested m with the eligible median outside: 3); F_W 0.3 (m 1.5: 5 of 50 (0 undefined); m 2: 13 of 50 (0 undefined); m 3: 25 of 50 (0 undefined); first tested m with the eligible median outside: 3). For Walls_Engraved multiplied on the odds scale in the large cluster (17 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 16 of 50 (0 undefined); m 2: 19 of 50 (0 undefined); m 3: 22 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.01 (m 1.5: 13 of 50 (0 undefined); m 2: 16 of 50 (0 undefined); m 3: 16 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.03 (m 1.5: 11 of 50 (0 undefined); m 2: 11 of 50 (0 undefined); m 3: 13 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.1 (m 1.5: 6 of 50 (0 undefined); m 2: 6 of 50 (0 undefined); m 3: 8 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.3 (m 1.5: 0 of 50 (0 undefined); m 2: 0 of 50 (0 undefined); m 3: 0 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid). For Walls_Engraved multiplied on the odds scale in the small cluster (11 sites), planted records with a defined statistic falling outside the unmarked band: F_W 0 (m 1.5: 12 of 50 (0 undefined); m 2: 14 of 50 (0 undefined); m 3: 15 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.01 (m 1.5: 12 of 50 (0 undefined); m 2: 11 of 50 (0 undefined); m 3: 11 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.03 (m 1.5: 8 of 50 (0 undefined); m 2: 9 of 50 (0 undefined); m 3: 9 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.1 (m 1.5: 6 of 50 (0 undefined); m 2: 6 of 50 (0 undefined); m 3: 7 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid); F_W 0.3 (m 1.5: 0 of 50 (0 undefined); m 2: 0 of 50 (0 undefined); m 3: 0 of 50 (0 undefined); first tested m with the eligible median outside: none on this grid). These comparisons do not identify the cause of the observed pattern and do not exclude marker patterns outside this planting design.
