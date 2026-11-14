# Type origination: calibration of C4 at theta = 0 against C3, by rank (exploratory)

Produced by `.venv-np254/bin/python analyses/per_type/type_origination_calibration.py`; NumPy 2.5.4; data hash 3b5e1cfabcad47f4. Both arms are read with the same function (`ranked_readouts`) at the baseline design (span limit L_d, dates independent of location). C3: the drift-extension records regenerated from their registered seeds. C4: the theta = 0 cells of `type_origination_cells.csv`. Each cell is median [10th, 90th percentile] over 50 records. The arms are independent sets of records, so medians differ by sampling as well as by construction. Descriptive bands from one seed set per arm, not posterior intervals.

By construction C3 loses a class from the regional target only at floating-point zero, and C4 loses it when its pool count reaches zero. C4 also starts from pooled shares rounded to the pool size, so a class allocated zero draws is absent from the start (listed per cell below).

## Primary point

### Window 8LR

| F_W | pool size | starting classes present | classes rounded to zero draws |
|---|---|---|---|
| 0.1 | 22172 | 10 | none |
| 0.3 | 6550 | 10 | none |
| 0.5 | 3371 | 10 | none |
| 0.7 | 1941 | 10 | none |
| 0.9 | 1015 | 10 | none |

| F_W | statistic | C3 by rank | C4 theta 0 by rank | C4 median in C3 band | C3 median in C4 band |
|---|---|---|---|---|---|
| 0.1 | F rank 1 | 0.021 [0.011, 0.037] | 0.020 [0.009, 0.045] | yes | yes |
| 0.1 | F rank 2 | 0.021 [0.011, 0.044] | 0.024 [0.014, 0.040] | yes | yes |
| 0.1 | F rank 3 | 0.026 [0.009, 0.040] | 0.020 [0.009, 0.053] | yes | yes |
| 0.1 | sites, rank 4 | 24 [11, 28] | 23 [12, 28] | yes | yes |
| 0.1 | sites, rank 5 | 12 [5, 24] | 13 [5, 23] | yes | yes |
| 0.1 | sites, rank 6 | 5 [2, 21] | 5 [1, 13] | yes | yes |
| 0.1 | sites, rank 7 | 2 [0, 7] | 2 [0, 7] | yes | yes |
| 0.1 | sites, rank 8 | 1 [0, 2] | 0 [0, 4] | yes | yes |
| 0.1 | sites, rank 9 | 0 [0, 1] | 0 [0, 1] | yes | yes |
| 0.1 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.1 | richness | 8 [6, 9] | 7 [6, 9] | yes | yes |
| 0.1 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.1 | D_dist | 0.107 [-0.012, 0.236] | 0.112 [0.034, 0.234] | yes | yes |
| 0.1 | |lat| rank 1 | 0.213 [0.045, 0.471] | 0.231 [0.048, 0.559] | yes | yes |
| 0.1 | |lat| rank 3 | 0.248 [0.024, 0.463] | 0.210 [0.044, 0.589] | yes | yes |
| 0.1 | |lat| rank 2 | 0.172 [0.035, 0.443] | 0.245 [0.050, 0.542] | yes | yes |
| 0.3 | F rank 1 | 0.036 [0.019, 0.108] | 0.046 [0.016, 0.123] | yes | yes |
| 0.3 | F rank 2 | 0.037 [0.014, 0.118] | 0.030 [0.011, 0.098] | yes | yes |
| 0.3 | F rank 3 | 0.032 [0.014, 0.117] | 0.036 [0.021, 0.080] | yes | yes |
| 0.3 | sites, rank 4 | 15 [7, 28] | 12 [3, 26] | yes | yes |
| 0.3 | sites, rank 5 | 7 [3, 19] | 6 [2, 15] | yes | yes |
| 0.3 | sites, rank 6 | 3 [1, 9] | 3 [1, 8] | yes | yes |
| 0.3 | sites, rank 7 | 1 [0, 6] | 1 [0, 4] | yes | yes |
| 0.3 | sites, rank 8 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.3 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | richness | 7 [6, 8] | 7 [6, 8] | yes | yes |
| 0.3 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.3 | D_dist | 0.073 [-0.018, 0.162] | 0.081 [0.020, 0.189] | yes | yes |
| 0.3 | |lat| rank 1 | 0.190 [0.036, 0.425] | 0.186 [0.058, 0.434] | yes | yes |
| 0.3 | |lat| rank 3 | 0.188 [0.036, 0.420] | 0.221 [0.046, 0.467] | yes | yes |
| 0.3 | |lat| rank 2 | 0.190 [0.032, 0.391] | 0.223 [0.064, 0.439] | yes | yes |
| 0.5 | F rank 1 | 0.070 [0.019, 0.180] | 0.058 [0.016, 0.167] | yes | yes |
| 0.5 | F rank 2 | 0.044 [0.021, 0.132] | 0.072 [0.022, 0.161] | yes | yes |
| 0.5 | F rank 3 | 0.041 [0.020, 0.107] | 0.044 [0.018, 0.090] | yes | yes |
| 0.5 | sites, rank 4 | 10 [2, 24] | 9 [3, 20] | yes | yes |
| 0.5 | sites, rank 5 | 4 [1, 10] | 4 [1, 9] | yes | yes |
| 0.5 | sites, rank 6 | 2 [1, 7] | 3 [0, 7] | yes | yes |
| 0.5 | sites, rank 7 | 1 [0, 4] | 1 [0, 4] | yes | yes |
| 0.5 | sites, rank 8 | 0 [0, 1] | 0 [0, 2] | yes | yes |
| 0.5 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | richness | 7 [6, 8] | 7 [5, 8] | yes | yes |
| 0.5 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.5 | D_dist | 0.043 [-0.032, 0.162] | 0.034 [-0.014, 0.152] | yes | yes |
| 0.5 | |lat| rank 1 | 0.166 [0.016, 0.446] | 0.194 [0.050, 0.434] | yes | yes |
| 0.5 | |lat| rank 3 | 0.218 [0.036, 0.439] | 0.199 [0.026, 0.439] | yes | yes |
| 0.5 | |lat| rank 2 | 0.180 [0.027, 0.411] | 0.197 [0.064, 0.386] | yes | yes |
| 0.7 | F rank 1 | 0.104 [0.025, 0.239] | 0.137 [0.042, 0.227] | yes | yes |
| 0.7 | F rank 2 | 0.090 [0.036, 0.168] | 0.078 [0.036, 0.176] | yes | yes |
| 0.7 | F rank 3 | 0.044 [0.025, 0.088] | 0.046 [0.025, 0.106] | yes | yes |
| 0.7 | sites, rank 4 | 8 [3, 17] | 7 [1, 20] | yes | yes |
| 0.7 | sites, rank 5 | 4 [1, 10] | 2 [0, 11] | yes | yes |
| 0.7 | sites, rank 6 | 2 [0, 7] | 2 [0, 7] | yes | yes |
| 0.7 | sites, rank 7 | 0 [0, 3] | 1 [0, 4] | yes | yes |
| 0.7 | sites, rank 8 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.7 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | richness | 6 [5, 8] | 7 [4, 8] | yes | yes |
| 0.7 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.7 | D_dist | 0.040 [-0.021, 0.174] | 0.035 [-0.028, 0.169] | yes | yes |
| 0.7 | |lat| rank 1 | 0.184 [0.032, 0.400] | 0.188 [0.041, 0.391] | yes | yes |
| 0.7 | |lat| rank 3 | 0.191 [0.037, 0.539] | 0.181 [0.045, 0.387] | yes | yes |
| 0.7 | |lat| rank 2 | 0.142 [0.019, 0.334] | 0.216 [0.082, 0.469] | yes | yes |
| 0.9 | F rank 1 | 0.159 [0.072, 0.367] | 0.194 [0.065, 0.293] | yes | yes |
| 0.9 | F rank 2 | 0.121 [0.049, 0.248] | 0.149 [0.048, 0.220] | yes | yes |
| 0.9 | F rank 3 | 0.049 [0.024, 0.124] | 0.041 [0.022, 0.100] | yes | yes |
| 0.9 | sites, rank 4 | 9 [4, 16] | 7 [4, 14] | yes | yes |
| 0.9 | sites, rank 5 | 5 [1, 12] | 4 [1, 9] | yes | yes |
| 0.9 | sites, rank 6 | 2 [0, 4] | 2 [0, 4] | yes | yes |
| 0.9 | sites, rank 7 | 1 [0, 4] | 0 [0, 3] | yes | yes |
| 0.9 | sites, rank 8 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.9 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | richness | 7 [5, 8] | 6 [5, 8] | yes | yes |
| 0.9 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.9 | D_dist | 0.021 [-0.033, 0.120] | 0.022 [-0.051, 0.101] | yes | yes |
| 0.9 | |lat| rank 1 | 0.179 [0.052, 0.368] | 0.173 [0.047, 0.362] | yes | yes |
| 0.9 | |lat| rank 3 | 0.271 [0.027, 0.514] | 0.170 [0.041, 0.430] | yes | yes |
| 0.9 | |lat| rank 2 | 0.169 [0.050, 0.404] | 0.182 [0.071, 0.362] | yes | yes |

### Window 32LR

| F_W | pool size | starting classes present | classes rounded to zero draws |
|---|---|---|---|
| 0.1 | 88686 | 10 | none |
| 0.3 | 26198 | 10 | none |
| 0.5 | 13481 | 10 | none |
| 0.7 | 7761 | 10 | none |
| 0.9 | 4059 | 10 | none |

| F_W | statistic | C3 by rank | C4 theta 0 by rank | C4 median in C3 band | C3 median in C4 band |
|---|---|---|---|---|---|
| 0.1 | F rank 1 | 0.031 [0.015, 0.049] | 0.027 [0.014, 0.042] | yes | yes |
| 0.1 | F rank 2 | 0.026 [0.015, 0.053] | 0.025 [0.015, 0.049] | yes | yes |
| 0.1 | F rank 3 | 0.029 [0.013, 0.048] | 0.029 [0.017, 0.042] | yes | yes |
| 0.1 | sites, rank 4 | 20 [5, 27] | 21 [7, 27] | yes | yes |
| 0.1 | sites, rank 5 | 9 [2, 24] | 9 [2, 23] | yes | yes |
| 0.1 | sites, rank 6 | 4 [1, 12] | 4 [1, 9] | yes | yes |
| 0.1 | sites, rank 7 | 2 [0, 5] | 2 [1, 6] | yes | yes |
| 0.1 | sites, rank 8 | 1 [0, 3] | 0 [0, 2] | yes | yes |
| 0.1 | sites, rank 9 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.1 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.1 | richness | 8 [6, 9] | 7 [7, 9] | yes | yes |
| 0.1 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.1 | D_dist | 0.032 [-0.040, 0.142] | 0.022 [-0.044, 0.118] | yes | yes |
| 0.1 | |lat| rank 1 | 0.188 [0.057, 0.380] | 0.131 [0.030, 0.295] | yes | yes |
| 0.1 | |lat| rank 3 | 0.157 [0.051, 0.381] | 0.150 [0.040, 0.320] | yes | yes |
| 0.1 | |lat| rank 2 | 0.204 [0.041, 0.359] | 0.120 [0.022, 0.322] | yes | yes |
| 0.3 | F rank 1 | 0.060 [0.027, 0.145] | 0.052 [0.020, 0.146] | yes | yes |
| 0.3 | F rank 2 | 0.058 [0.031, 0.150] | 0.058 [0.027, 0.148] | yes | yes |
| 0.3 | F rank 3 | 0.046 [0.028, 0.078] | 0.051 [0.021, 0.111] | yes | yes |
| 0.3 | sites, rank 4 | 9 [2, 27] | 9 [2, 22] | yes | yes |
| 0.3 | sites, rank 5 | 4 [1, 10] | 3 [1, 11] | yes | yes |
| 0.3 | sites, rank 6 | 2 [0, 5] | 2 [0, 4] | yes | yes |
| 0.3 | sites, rank 7 | 1 [0, 4] | 0 [0, 2] | yes | yes |
| 0.3 | sites, rank 8 | 0 [0, 1] | 0 [0, 1] | yes | yes |
| 0.3 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | richness | 7 [5, 8] | 6 [5, 8] | yes | yes |
| 0.3 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.3 | D_dist | 0.037 [-0.031, 0.150] | 0.011 [-0.045, 0.151] | yes | yes |
| 0.3 | |lat| rank 1 | 0.165 [0.042, 0.377] | 0.136 [0.037, 0.373] | yes | yes |
| 0.3 | |lat| rank 3 | 0.181 [0.039, 0.405] | 0.166 [0.042, 0.384] | yes | yes |
| 0.3 | |lat| rank 2 | 0.158 [0.020, 0.387] | 0.187 [0.020, 0.360] | yes | yes |
| 0.5 | F rank 1 | 0.099 [0.034, 0.238] | 0.100 [0.038, 0.230] | yes | yes |
| 0.5 | F rank 2 | 0.065 [0.028, 0.195] | 0.078 [0.037, 0.200] | yes | yes |
| 0.5 | F rank 3 | 0.077 [0.036, 0.165] | 0.051 [0.033, 0.109] | yes | yes |
| 0.5 | sites, rank 4 | 6 [2, 15] | 4 [1, 13] | yes | yes |
| 0.5 | sites, rank 5 | 3 [1, 5] | 2 [0, 6] | yes | yes |
| 0.5 | sites, rank 6 | 1 [0, 2] | 1 [0, 3] | yes | yes |
| 0.5 | sites, rank 7 | 0 [0, 2] | 0 [0, 2] | yes | yes |
| 0.5 | sites, rank 8 | 0 [0, 1] | 0 [0, 1] | yes | yes |
| 0.5 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | richness | 6 [5, 8] | 6 [4, 8] | yes | yes |
| 0.5 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.5 | D_dist | 0.007 [-0.045, 0.078] | 0.001 [-0.044, 0.077] | yes | yes |
| 0.5 | |lat| rank 1 | 0.164 [0.010, 0.309] | 0.109 [0.018, 0.297] | yes | yes |
| 0.5 | |lat| rank 3 | 0.131 [0.037, 0.322] | 0.140 [0.020, 0.327] | yes | yes |
| 0.5 | |lat| rank 2 | 0.159 [0.035, 0.330] | 0.106 [0.022, 0.333] | yes | yes |
| 0.7 | F rank 1 | 0.162 [0.039, 0.326] | 0.151 [0.056, 0.353] | yes | yes |
| 0.7 | F rank 2 | 0.135 [0.045, 0.254] | 0.102 [0.040, 0.256] | yes | yes |
| 0.7 | F rank 3 | 0.071 [0.022, 0.205] | 0.083 [0.041, 0.202] | yes | yes |
| 0.7 | sites, rank 4 | 4 [1, 12] | 4 [1, 11] | yes | yes |
| 0.7 | sites, rank 5 | 1 [0, 4] | 2 [0, 6] | yes | yes |
| 0.7 | sites, rank 6 | 0 [0, 2] | 1 [0, 2] | yes | yes |
| 0.7 | sites, rank 7 | 0 [0, 1] | 0 [0, 1] | yes | yes |
| 0.7 | sites, rank 8 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | richness | 5 [4, 7] | 6 [4, 7] | yes | yes |
| 0.7 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.7 | D_dist | 0.011 [-0.058, 0.081] | 0.011 [-0.054, 0.102] | yes | yes |
| 0.7 | |lat| rank 1 | 0.120 [0.014, 0.356] | 0.158 [0.028, 0.275] | yes | yes |
| 0.7 | |lat| rank 3 | 0.139 [0.039, 0.311] | 0.145 [0.034, 0.317] | yes | yes |
| 0.7 | |lat| rank 2 | 0.123 [0.016, 0.359] | 0.131 [0.010, 0.270] | yes | yes |
| 0.9 | F rank 1 | 0.286 [0.100, 0.436] | 0.270 [0.074, 0.450] | yes | yes |
| 0.9 | F rank 2 | 0.192 [0.099, 0.383] | 0.176 [0.072, 0.330] | yes | yes |
| 0.9 | F rank 3 | 0.092 [0.033, 0.234] | 0.101 [0.030, 0.289] | yes | yes |
| 0.9 | sites, rank 4 | 4 [1, 9] | 4 [1, 9] | yes | yes |
| 0.9 | sites, rank 5 | 2 [0, 5] | 2 [0, 4] | yes | yes |
| 0.9 | sites, rank 6 | 1 [0, 2] | 0 [0, 2] | yes | yes |
| 0.9 | sites, rank 7 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.9 | sites, rank 8 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | richness | 6 [4, 7] | 5 [4, 7] | yes | yes |
| 0.9 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.9 | D_dist | -0.003 [-0.067, 0.085] | -0.012 [-0.066, 0.072] | yes | yes |
| 0.9 | |lat| rank 1 | 0.109 [0.016, 0.306] | 0.126 [0.027, 0.341] | yes | yes |
| 0.9 | |lat| rank 3 | 0.130 [0.027, 0.303] | 0.144 [0.030, 0.289] | yes | yes |
| 0.9 | |lat| rank 2 | 0.101 [0.022, 0.327] | 0.140 [0.032, 0.340] | yes | yes |

## Secondary point

### Window 8LR

| F_W | pool size | starting classes present | classes rounded to zero draws |
|---|---|---|---|
| 0.1 | 3038 | 10 | none |
| 0.3 | 898 | 10 | none |
| 0.5 | 462 | 10 | none |
| 0.7 | 266 | 9 | Hull_Engraved |
| 0.9 | 139 | 8 | Vernon_Paul_Applique, Hull_Engraved |

| F_W | statistic | C3 by rank | C4 theta 0 by rank | C4 median in C3 band | C3 median in C4 band |
|---|---|---|---|---|---|
| 0.1 | F rank 1 | 0.023 [0.012, 0.054] | 0.027 [0.016, 0.051] | yes | yes |
| 0.1 | F rank 2 | 0.024 [0.012, 0.052] | 0.025 [0.015, 0.045] | yes | yes |
| 0.1 | F rank 3 | 0.021 [0.013, 0.057] | 0.026 [0.013, 0.049] | yes | yes |
| 0.1 | sites, rank 4 | 20 [7, 26] | 18 [7, 26] | yes | yes |
| 0.1 | sites, rank 5 | 8 [3, 23] | 8 [2, 22] | yes | yes |
| 0.1 | sites, rank 6 | 4 [1, 15] | 4 [1, 12] | yes | yes |
| 0.1 | sites, rank 7 | 2 [0, 9] | 1 [0, 6] | yes | yes |
| 0.1 | sites, rank 8 | 0 [0, 5] | 0 [0, 3] | yes | yes |
| 0.1 | sites, rank 9 | 0 [0, 1] | 0 [0, 1] | yes | yes |
| 0.1 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.1 | richness | 8 [6, 9] | 7 [6, 9] | yes | yes |
| 0.1 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.1 | D_dist | 0.054 [-0.039, 0.133] | 0.004 [-0.041, 0.131] | yes | yes |
| 0.1 | |lat| rank 1 | 0.162 [0.036, 0.377] | 0.178 [0.056, 0.348] | yes | yes |
| 0.1 | |lat| rank 3 | 0.150 [0.048, 0.329] | 0.207 [0.071, 0.409] | yes | yes |
| 0.1 | |lat| rank 2 | 0.160 [0.031, 0.404] | 0.143 [0.017, 0.378] | yes | yes |
| 0.3 | F rank 1 | 0.052 [0.023, 0.116] | 0.036 [0.019, 0.130] | yes | yes |
| 0.3 | F rank 2 | 0.045 [0.023, 0.114] | 0.036 [0.014, 0.083] | yes | yes |
| 0.3 | F rank 3 | 0.040 [0.018, 0.079] | 0.041 [0.021, 0.075] | yes | yes |
| 0.3 | sites, rank 4 | 12 [5, 27] | 10 [3, 27] | yes | yes |
| 0.3 | sites, rank 5 | 4 [1, 11] | 4 [1, 13] | yes | yes |
| 0.3 | sites, rank 6 | 2 [0, 7] | 2 [0, 5] | yes | yes |
| 0.3 | sites, rank 7 | 1 [0, 4] | 0 [0, 3] | yes | yes |
| 0.3 | sites, rank 8 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.3 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | richness | 7 [5, 8] | 6 [5, 8] | yes | yes |
| 0.3 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.3 | D_dist | 0.018 [-0.049, 0.130] | 0.034 [-0.042, 0.107] | yes | yes |
| 0.3 | |lat| rank 1 | 0.148 [0.033, 0.412] | 0.139 [0.049, 0.420] | yes | yes |
| 0.3 | |lat| rank 3 | 0.192 [0.051, 0.441] | 0.185 [0.015, 0.491] | yes | yes |
| 0.3 | |lat| rank 2 | 0.192 [0.019, 0.433] | 0.242 [0.028, 0.379] | yes | yes |
| 0.5 | F rank 1 | 0.081 [0.027, 0.222] | 0.083 [0.029, 0.193] | yes | yes |
| 0.5 | F rank 2 | 0.087 [0.035, 0.181] | 0.069 [0.034, 0.174] | yes | yes |
| 0.5 | F rank 3 | 0.048 [0.027, 0.102] | 0.061 [0.025, 0.143] | yes | yes |
| 0.5 | sites, rank 4 | 6 [2, 19] | 8 [2, 17] | yes | yes |
| 0.5 | sites, rank 5 | 3 [1, 11] | 4 [1, 13] | yes | yes |
| 0.5 | sites, rank 6 | 1 [0, 4] | 2 [0, 5] | yes | yes |
| 0.5 | sites, rank 7 | 0 [0, 2] | 1 [0, 3] | yes | yes |
| 0.5 | sites, rank 8 | 0 [0, 0] | 0 [0, 1] | yes | yes |
| 0.5 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | richness | 6 [5, 7] | 7 [5, 8] | yes | yes |
| 0.5 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.5 | D_dist | 0.003 [-0.058, 0.079] | 0.010 [-0.030, 0.106] | yes | yes |
| 0.5 | |lat| rank 1 | 0.135 [0.014, 0.280] | 0.160 [0.051, 0.323] | yes | yes |
| 0.5 | |lat| rank 3 | 0.150 [0.027, 0.354] | 0.149 [0.036, 0.450] | yes | yes |
| 0.5 | |lat| rank 2 | 0.153 [0.013, 0.282] | 0.147 [0.021, 0.330] | yes | yes |
| 0.7 | F rank 1 | 0.161 [0.028, 0.319] | 0.127 [0.040, 0.297] | yes | yes |
| 0.7 | F rank 2 | 0.124 [0.041, 0.219] | 0.094 [0.030, 0.264] | yes | yes |
| 0.7 | F rank 3 | 0.051 [0.026, 0.157] | 0.056 [0.025, 0.246] | yes | yes |
| 0.7 | sites, rank 4 | 5 [2, 16] | 5 [2, 12] | yes | yes |
| 0.7 | sites, rank 5 | 3 [0, 9] | 3 [1, 8] | yes | yes |
| 0.7 | sites, rank 6 | 1 [0, 4] | 1 [0, 4] | yes | yes |
| 0.7 | sites, rank 7 | 0 [0, 2] | 0 [0, 2] | yes | yes |
| 0.7 | sites, rank 8 | 0 [0, 0] | 0 [0, 1] | yes | yes |
| 0.7 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | richness | 6 [4, 7] | 6 [5, 8] | yes | yes |
| 0.7 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.7 | D_dist | 0.017 [-0.066, 0.083] | 0.003 [-0.059, 0.095] | yes | yes |
| 0.7 | |lat| rank 1 | 0.127 [0.007, 0.318] | 0.111 [0.022, 0.252] | yes | yes |
| 0.7 | |lat| rank 3 | 0.169 [0.037, 0.340] | 0.123 [0.022, 0.357] | yes | yes |
| 0.7 | |lat| rank 2 | 0.141 [0.023, 0.332] | 0.119 [0.023, 0.297] | yes | yes |
| 0.9 | F rank 1 | 0.217 [0.058, 0.485] | 0.206 [0.043, 0.374] | yes | yes |
| 0.9 | F rank 2 | 0.128 [0.058, 0.283] | 0.163 [0.063, 0.316] | yes | yes |
| 0.9 | F rank 3 | 0.075 [0.030, 0.240] | 0.057 [0.021, 0.134] | yes | yes |
| 0.9 | sites, rank 4 | 4 [1, 12] | 6 [2, 15] | yes | yes |
| 0.9 | sites, rank 5 | 3 [0, 6] | 3 [0, 7] | yes | yes |
| 0.9 | sites, rank 6 | 1 [0, 5] | 1 [0, 3] | yes | yes |
| 0.9 | sites, rank 7 | 0 [0, 3] | 0 [0, 2] | yes | yes |
| 0.9 | sites, rank 8 | 0 [0, 1] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | richness | 6 [4, 8] | 6 [4, 7] | yes | yes |
| 0.9 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.9 | D_dist | 0.007 [-0.058, 0.093] | 0.007 [-0.049, 0.098] | yes | yes |
| 0.9 | |lat| rank 1 | 0.136 [0.038, 0.274] | 0.127 [0.011, 0.350] | yes | yes |
| 0.9 | |lat| rank 3 | 0.175 [0.032, 0.343] | 0.159 [0.025, 0.316] | yes | yes |
| 0.9 | |lat| rank 2 | 0.120 [0.014, 0.320] | 0.117 [0.028, 0.351] | yes | yes |

### Window 32LR

| F_W | pool size | starting classes present | classes rounded to zero draws |
|---|---|---|---|
| 0.1 | 12149 | 10 | none |
| 0.3 | 3589 | 10 | none |
| 0.5 | 1847 | 10 | none |
| 0.7 | 1064 | 10 | none |
| 0.9 | 556 | 10 | none |

| F_W | statistic | C3 by rank | C4 theta 0 by rank | C4 median in C3 band | C3 median in C4 band |
|---|---|---|---|---|---|
| 0.1 | F rank 1 | 0.033 [0.018, 0.066] | 0.031 [0.016, 0.050] | yes | yes |
| 0.1 | F rank 2 | 0.036 [0.021, 0.056] | 0.032 [0.018, 0.056] | yes | yes |
| 0.1 | F rank 3 | 0.027 [0.018, 0.050] | 0.030 [0.016, 0.054] | yes | yes |
| 0.1 | sites, rank 4 | 20 [6, 27] | 13 [4, 26] | yes | yes |
| 0.1 | sites, rank 5 | 8 [2, 21] | 8 [2, 20] | yes | yes |
| 0.1 | sites, rank 6 | 4 [1, 11] | 4 [1, 16] | yes | yes |
| 0.1 | sites, rank 7 | 2 [0, 6] | 2 [0, 6] | yes | yes |
| 0.1 | sites, rank 8 | 1 [0, 3] | 1 [0, 4] | yes | yes |
| 0.1 | sites, rank 9 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.1 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.1 | richness | 8 [6, 9] | 8 [6, 9] | yes | yes |
| 0.1 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.1 | D_dist | 0.039 [-0.045, 0.099] | 0.035 [-0.026, 0.105] | yes | yes |
| 0.1 | |lat| rank 1 | 0.127 [0.027, 0.357] | 0.157 [0.052, 0.326] | yes | yes |
| 0.1 | |lat| rank 3 | 0.170 [0.013, 0.377] | 0.174 [0.068, 0.353] | yes | yes |
| 0.1 | |lat| rank 2 | 0.158 [0.008, 0.392] | 0.173 [0.078, 0.281] | yes | yes |
| 0.3 | F rank 1 | 0.058 [0.035, 0.158] | 0.060 [0.030, 0.160] | yes | yes |
| 0.3 | F rank 2 | 0.048 [0.027, 0.124] | 0.057 [0.023, 0.151] | yes | yes |
| 0.3 | F rank 3 | 0.054 [0.024, 0.102] | 0.052 [0.026, 0.097] | yes | yes |
| 0.3 | sites, rank 4 | 8 [2, 26] | 6 [2, 26] | yes | yes |
| 0.3 | sites, rank 5 | 3 [0, 10] | 3 [1, 13] | yes | yes |
| 0.3 | sites, rank 6 | 1 [0, 3] | 2 [0, 4] | yes | yes |
| 0.3 | sites, rank 7 | 0 [0, 2] | 1 [0, 4] | yes | yes |
| 0.3 | sites, rank 8 | 0 [0, 1] | 0 [0, 1] | yes | yes |
| 0.3 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.3 | richness | 6 [4, 8] | 7 [5, 8] | yes | yes |
| 0.3 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.3 | D_dist | 0.006 [-0.055, 0.077] | 0.010 [-0.044, 0.117] | yes | yes |
| 0.3 | |lat| rank 1 | 0.108 [0.025, 0.322] | 0.153 [0.060, 0.303] | yes | yes |
| 0.3 | |lat| rank 3 | 0.111 [0.021, 0.311] | 0.145 [0.032, 0.294] | yes | yes |
| 0.3 | |lat| rank 2 | 0.110 [0.017, 0.324] | 0.122 [0.019, 0.323] | yes | yes |
| 0.5 | F rank 1 | 0.111 [0.046, 0.239] | 0.097 [0.037, 0.276] | yes | yes |
| 0.5 | F rank 2 | 0.094 [0.039, 0.195] | 0.095 [0.045, 0.216] | yes | yes |
| 0.5 | F rank 3 | 0.069 [0.032, 0.132] | 0.070 [0.034, 0.165] | yes | yes |
| 0.5 | sites, rank 4 | 5 [1, 15] | 6 [1, 19] | yes | yes |
| 0.5 | sites, rank 5 | 2 [0, 6] | 1 [0, 7] | yes | yes |
| 0.5 | sites, rank 6 | 1 [0, 3] | 0 [0, 3] | yes | yes |
| 0.5 | sites, rank 7 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.5 | sites, rank 8 | 0 [0, 1] | 0 [0, 0] | yes | yes |
| 0.5 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.5 | richness | 6 [4, 8] | 5 [4, 7] | yes | yes |
| 0.5 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.5 | D_dist | 0.012 [-0.037, 0.089] | 0.012 [-0.062, 0.111] | yes | yes |
| 0.5 | |lat| rank 1 | 0.134 [0.040, 0.338] | 0.163 [0.034, 0.416] | yes | yes |
| 0.5 | |lat| rank 3 | 0.134 [0.025, 0.294] | 0.140 [0.027, 0.387] | yes | yes |
| 0.5 | |lat| rank 2 | 0.136 [0.038, 0.277] | 0.146 [0.024, 0.436] | yes | yes |
| 0.7 | F rank 1 | 0.229 [0.071, 0.407] | 0.208 [0.076, 0.397] | yes | yes |
| 0.7 | F rank 2 | 0.165 [0.091, 0.294] | 0.171 [0.061, 0.302] | yes | yes |
| 0.7 | F rank 3 | 0.096 [0.042, 0.196] | 0.091 [0.038, 0.158] | yes | yes |
| 0.7 | sites, rank 4 | 3 [1, 8] | 4 [1, 10] | yes | yes |
| 0.7 | sites, rank 5 | 1 [0, 4] | 1 [0, 5] | yes | yes |
| 0.7 | sites, rank 6 | 0 [0, 2] | 0 [0, 3] | yes | yes |
| 0.7 | sites, rank 7 | 0 [0, 1] | 0 [0, 1] | yes | yes |
| 0.7 | sites, rank 8 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.7 | richness | 5 [4, 7] | 5 [4, 7] | yes | yes |
| 0.7 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.7 | D_dist | -0.011 [-0.051, 0.047] | -0.001 [-0.064, 0.093] | yes | yes |
| 0.7 | |lat| rank 1 | 0.109 [0.027, 0.349] | 0.105 [0.012, 0.326] | yes | yes |
| 0.7 | |lat| rank 3 | 0.146 [0.023, 0.336] | 0.149 [0.032, 0.299] | yes | yes |
| 0.7 | |lat| rank 2 | 0.133 [0.045, 0.315] | 0.102 [0.016, 0.363] | yes | yes |
| 0.9 | F rank 1 | 0.329 [0.092, 0.548] | 0.295 [0.132, 0.453] | yes | yes |
| 0.9 | F rank 2 | 0.220 [0.074, 0.411] | 0.224 [0.099, 0.397] | yes | yes |
| 0.9 | F rank 3 | 0.090 [0.022, 0.262] | 0.117 [0.037, 0.265] | yes | yes |
| 0.9 | sites, rank 4 | 2 [0, 6] | 2 [0, 5] | yes | yes |
| 0.9 | sites, rank 5 | 1 [0, 4] | 1 [0, 2] | yes | yes |
| 0.9 | sites, rank 6 | 0 [0, 2] | 0 [0, 1] | yes | yes |
| 0.9 | sites, rank 7 | 0 [0, 1] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 8 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 9 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | sites, rank 10 | 0 [0, 0] | 0 [0, 0] | yes | yes |
| 0.9 | richness | 5 [3, 7] | 5 [3, 6] | yes | yes |
| 0.9 | share beyond rank 10 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | yes | yes |
| 0.9 | D_dist | -0.003 [-0.061, 0.085] | -0.003 [-0.054, 0.070] | yes | yes |
| 0.9 | |lat| rank 1 | 0.157 [0.013, 0.313] | 0.153 [0.047, 0.299] | yes | yes |
| 0.9 | |lat| rank 3 | 0.123 [0.032, 0.270] | 0.134 [0.029, 0.339] | yes | yes |
| 0.9 | |lat| rank 2 | 0.139 [0.029, 0.338] | 0.168 [0.053, 0.308] | yes | yes |

## Overlap by statistic

Cells are the 20 combinations of point, window and amplitude. Each check asks whether one arm's median lies inside the other arm's 10th to 90th percentile band. The check is coarse: two bands shifted by less than half their width pass it, so passing shows overlap and does not show that the generators agree to any stated tolerance (Codex repair check 2026-11-14, RC-F2). No tolerance for a material difference has been set. The last two columns give the size of what was observed: the largest absolute difference between the two medians over the 20 cells, and the largest ratio of that difference to the C3 band's width, over the cells where the width is positive.

| statistic | C4 median inside C3 band | C3 median inside C4 band | cells where neither holds | largest median difference | as a fraction of the C3 band width |
|---|---|---|---|---|---|
| F rank 1 | 20 of 20 | 20 of 20 | none | 0.036 | 0.18 |
| F rank 2 | 20 of 20 | 20 of 20 | none | 0.035 | 0.25 |
| F rank 3 | 20 of 20 | 20 of 20 | none | 0.027 | 0.21 |
| sites, rank 4 | 20 of 20 | 20 of 20 | none | 6 | 0.31 |
| sites, rank 5 | 20 of 20 | 20 of 20 | none | 2 | 0.25 |
| sites, rank 6 | 20 of 20 | 20 of 20 | none | 1 | 0.50 |
| sites, rank 7 | 20 of 20 | 20 of 20 | none | 1 | 0.50 |
| sites, rank 8 | 20 of 20 | 20 of 20 | none | 1 | 0.50 (5 cells have a C3 band of zero width and equal medians) |
| sites, rank 9 | 20 of 20 | 20 of 20 | none | 0 | 0.00 (16 cells have a C3 band of zero width and equal medians) |
| sites, rank 10 | 20 of 20 | 20 of 20 | none | 0 | 0.00 (19 cells have a C3 band of zero width and equal medians) |
| richness | 20 of 20 | 20 of 20 | none | 1 | 0.50 |
| share beyond rank 10 | 20 of 20 | 20 of 20 | none | 0.000 | not defined (20 cells have a C3 band of zero width and equal medians) |
| D_dist | 20 of 20 | 20 of 20 | none | 0.050 | 0.29 |
| |lat| rank 1 | 20 of 20 | 20 of 20 | none | 0.057 | 0.19 |
| |lat| rank 3 | 20 of 20 | 20 of 20 | none | 0.101 | 0.21 |
| |lat| rank 2 | 20 of 20 | 20 of 20 | none | 0.083 | 0.26 |

## Status of the record under each generator

The author's definition of a material difference (2026-11-14): one that changes the record's status (inside, above or below the band) for a statistic in a cell, since status is the reading rule of both designs. The table counts, over the 20 cells, how often the record has the same status against the C3 band and against the C4 band, and lists the cells where it does not. Each band is the 10th to 90th percentile of 50 records from one seed set, and the two arms are independent, so a record value near a band edge can change status by sampling alone. The count is a description, not a test.

| statistic | same status | cells where the status differs (C3 status, C4 status) |
|---|---|---|
| F rank 1 | 20 of 20 | none |
| F rank 2 | 19 of 20 | primary 8LR F_W 0.3 (inside, above) |
| F rank 3 | 18 of 20 | secondary 8LR F_W 0.7 (above, inside); secondary 8LR F_W 0.9 (inside, above) |
| sites, rank 4 | 19 of 20 | primary 8LR F_W 0.5 (inside, above) |
| sites, rank 5 | 20 of 20 | none |
| sites, rank 6 | 17 of 20 | primary 8LR F_W 0.1 (inside, above); secondary 32LR F_W 0.1 (above, inside); secondary 8LR F_W 0.1 (inside, above) |
| sites, rank 7 | 20 of 20 | none |
| sites, rank 8 | 20 of 20 | none |
| sites, rank 9 | 20 of 20 | none |
| sites, rank 10 | 20 of 20 | none |
| richness | 20 of 20 | none |
| share beyond rank 10 | 20 of 20 | none |
| D_dist | 20 of 20 | none |
| |lat| rank 1 | 19 of 20 | primary 8LR F_W 0.1 (above, inside) |
| |lat| rank 3 | 20 of 20 | none |
| |lat| rank 2 | 15 of 20 | primary 8LR F_W 0.3 (inside, below); primary 8LR F_W 0.5 (inside, below); primary 8LR F_W 0.7 (inside, below); primary 8LR F_W 0.9 (inside, below); secondary 32LR F_W 0.1 (inside, below) |
| all statistics | 307 of 320 | |

Common-class F (ranks 1 to 3): same status in 57 of 60. Rare-class occupancy (ranks 4 to 10): same status in 136 of 140.
