# A3 identification test with the placement supplied (step 6b; design section 4)

Produced by `analyses/probe_latent_groups/validate_a3.R --alpha-sd=2` on 2026-11-23. Records from `a3_contracts.py --alpha-sd 2`, fitted under the Normal(0, 2) class-intercept prior. Exploratory design validation, not evidence about the record.

## Code equivalence

A0 against A3 with omega fixed at 0 on `a3_own_G1_om0_r1`: largest |z| over 41 shared posterior means 2.64 at `log_lik[10]` (threshold 3).
- A0: R-hat max 1.0023 (tau[1]) | bulk ESS min 1451 (lp__) | tail ESS min 2647 (lp__) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.721 | stuck 0
- A3 fixed 0: R-hat max 1.0057 (lp__) | bulk ESS min 1047 (lp__) | tail ESS min 1791 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.725 | stuck 0

## omega against the truth, placement supplied

Prior half-Normal(0, 2): median 1.35, 90 percent below 3.3. A posterior that stays near that, whatever the planted value, has learned little about omega.

| record | kind | placement | truth omega / leak / rho | between-region F_ST | omega post | sigma | rho | tau | share | gate | rung | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| own_G1_om0_r1 | own-generator A3 | G1 | omega 0.00, rho 203 | 0.1348 | 1.126 [0.187, 2.548] | 3.741 [2.559, 4.795] | 183.307 [80.958, 283.025] | 0.254 [0.030, 0.520] | 0.995 [0.975, 1.000] | pass | 2 | 558 |
| own_G1_om0_r2 | own-generator A3 | G1 | omega 0.00, rho 203 | 0.0315 | 0.530 [0.063, 1.557] | 4.130 [3.089, 5.120] | 196.518 [106.832, 286.050] | 0.422 [0.142, 0.654] | 0.990 [0.969, 0.999] | pass | 1 | 146 |
| own_G1_om0_r3 | own-generator A3 | G1 | omega 0.00, rho 203 | 0.0095 | 0.607 [0.059, 1.691] | 3.538 [2.655, 4.474] | 221.023 [120.013, 291.604] | 0.510 [0.265, 0.709] | 0.979 [0.946, 0.996] | pass | 1 | 147 |
| own_G1_om0.25_r1 | own-generator A3 | G1 | omega 0.25, rho 203 | 0.0738 | 0.746 [0.081, 1.960] | 3.961 [2.836, 5.180] | 148.545 [72.487, 261.674] | 0.251 [0.024, 0.582] | 0.996 [0.975, 1.000] | pass | 1 | 119 |
| own_G1_om0.25_r2 | own-generator A3 | G1 | omega 0.25, rho 203 | 0.0122 | 0.841 [0.089, 2.086] | 4.037 [3.001, 5.135] | 186.622 [95.205, 283.480] | 0.412 [0.101, 0.668] | 0.990 [0.965, 0.999] | pass | 1 | 161 |
| own_G1_om0.25_r3 | own-generator A3 | G1 | omega 0.25, rho 203 | 0.1327 | 0.837 [0.094, 2.242] | 3.849 [2.673, 5.005] | 177.264 [79.911, 280.221] | 0.346 [0.041, 0.660] | 0.992 [0.961, 1.000] | pass | 2 | 556 |
| own_G1_om0.5_r1 | own-generator A3 | G1 | omega 0.50, rho 203 | 0.0948 | 0.963 [0.135, 2.347] | 3.310 [2.332, 4.301] | 231.727 [123.287, 293.743] | 0.699 [0.490, 0.894] | 0.957 [0.892, 0.985] | pass | 1 | 134 |
| own_G1_om0.5_r2 | own-generator A3 | G1 | omega 0.50, rho 203 | 0.3625 | 1.374 [0.402, 2.685] | 3.841 [2.722, 4.945] | 182.845 [83.363, 280.851] | 0.516 [0.184, 0.755] | 0.983 [0.948, 0.998] | pass | 1 | 108 |
| own_G1_om0.5_r3 | own-generator A3 | G1 | omega 0.50, rho 203 | 0.0237 | 0.796 [0.090, 2.077] | 3.711 [2.664, 4.716] | 190.278 [95.745, 283.858] | 0.492 [0.292, 0.688] | 0.983 [0.952, 0.995] | pass | 1 | 153 |
| own_G1_om1_r1 | own-generator A3 | G1 | omega 1.00, rho 203 | 0.0191 | 1.362 [0.387, 2.828] | 3.308 [1.998, 4.445] | 170.904 [56.767, 283.013] | 0.587 [0.336, 0.827] | 0.969 [0.900, 0.992] | pass | 1 | 155 |
| own_G1_om1_r2 | own-generator A3 | G1 | omega 1.00, rho 203 | 0.0188 | 1.265 [0.336, 2.494] | 3.824 [2.777, 4.892] | 197.945 [96.641, 287.984] | 0.571 [0.284, 0.804] | 0.978 [0.940, 0.995] | pass | 2 | 557 |
| own_G1_om1_r3 | own-generator A3 | G1 | omega 1.00, rho 203 | 0.2004 | 2.585 [1.473, 4.370] | 3.607 [2.249, 4.844] | 195.182 [71.111, 289.039] | 0.673 [0.442, 0.890] | 0.966 [0.901, 0.989] | pass | 1 | 113 |
| own_L_om0_r1 | own-generator A3 | L | omega 0.00, rho 203 | 0.0856 | 1.051 [0.180, 2.317] | 3.754 [2.709, 4.823] | 182.141 [88.428, 282.493] | 0.378 [0.097, 0.600] | 0.990 [0.968, 0.999] | pass | 1 | 105 |
| own_L_om0_r2 | own-generator A3 | L | omega 0.00, rho 203 | 0.1457 | 0.615 [0.066, 1.731] | 4.025 [2.974, 5.222] | 153.776 [80.095, 266.458] | 0.260 [0.034, 0.539] | 0.996 [0.980, 1.000] | pass | 1 | 156 |
| own_L_om0_r3 | own-generator A3 | L | omega 0.00, rho 203 | 0.2496 | 0.707 [0.082, 1.815] | 3.615 [2.626, 4.550] | 205.195 [105.849, 287.867] | 0.292 [0.040, 0.548] | 0.994 [0.970, 1.000] | pass | 1 | 106 |
| own_L_om0.25_r1 | own-generator A3 | L | omega 0.25, rho 203 | 0.0413 | 0.697 [0.078, 1.768] | 3.559 [2.614, 4.584] | 209.548 [105.522, 289.504] | 0.542 [0.242, 0.775] | 0.977 [0.937, 0.996] | pass | 1 | 163 |
| own_L_om0.25_r2 | own-generator A3 | L | omega 0.25, rho 203 | 0.0462 | 0.445 [0.035, 1.411] | 3.789 [2.764, 4.878] | 166.010 [84.085, 273.210] | 0.491 [0.226, 0.716] | 0.984 [0.954, 0.997] | pass | 1 | 146 |
| own_L_om0.25_r3 | own-generator A3 | L | omega 0.25, rho 203 | 0.1016 | 0.821 [0.080, 2.218] | 3.901 [2.709, 5.173] | 154.603 [63.927, 272.833] | 0.682 [0.315, 0.990] | 0.971 [0.915, 0.994] | pass | 1 | 115 |
| own_L_om0.5_r1 | own-generator A3 | L | omega 0.50, rho 203 | 0.1211 | 2.478 [1.202, 4.466] | 3.222 [1.355, 4.739] | 123.404 [15.141, 270.702] | 0.292 [0.032, 0.623] | 0.991 [0.935, 1.000] | pass | 1 | 477 |
| own_L_om0.5_r2 | own-generator A3 | L | omega 0.50, rho 203 | 0.0763 | 0.654 [0.075, 1.897] | 4.329 [3.261, 5.492] | 179.234 [97.830, 279.645] | 0.276 [0.032, 0.617] | 0.996 [0.976, 1.000] | pass | 1 | 184 |
| own_L_om0.5_r3 | own-generator A3 | L | omega 0.50, rho 203 | 0.2408 | 0.695 [0.082, 1.814] | 3.913 [2.854, 5.073] | 160.932 [82.022, 270.222] | 0.197 [0.018, 0.504] | 0.997 [0.981, 1.000] | pass | 1 | 117 |
| own_L_om1_r1 | own-generator A3 | L | omega 1.00, rho 203 | 0.2099 | 2.137 [0.883, 3.837] | 3.540 [2.128, 5.044] | 106.460 [31.206, 228.708] | 0.366 [0.047, 0.686] | 0.990 [0.950, 1.000] | pass | 1 | 149 |
| own_L_om1_r2 | own-generator A3 | L | omega 1.00, rho 203 | 0.2136 | 1.014 [0.136, 2.437] | 4.106 [3.062, 5.258] | 174.366 [92.194, 278.820] | 0.157 [0.016, 0.418] | 0.999 [0.988, 1.000] | pass | 1 | 137 |
| own_L_om1_r3 | own-generator A3 | L | omega 1.00, rho 203 | 0.0074 | 0.797 [0.090, 2.179] | 4.294 [3.143, 5.504] | 168.244 [78.376, 275.086] | 0.595 [0.152, 0.886] | 0.981 [0.946, 0.999] | pass | 1 | 154 |
| sim_S0_r1 | copying simulator S0, leak 1.0, T1, primary operating point | G1 | leak 1.00 | 0.0005 | 2.694 [0.543, 4.665] | 0.787 [0.332, 1.394] | 163.525 [21.901, 290.433] | 0.078 [0.008, 0.177] | 0.990 [0.888, 1.000] | pass | 2 | 491 |
| sim_S0_r2 | copying simulator S0, leak 1.0, T1, primary operating point | G1 | leak 1.00 | 0.0073 | 1.502 [0.552, 3.038] | 2.734 [2.084, 3.372] | 253.289 [156.706, 296.531] | 0.083 [0.007, 0.229] | 0.999 [0.992, 1.000] | pass | 1 | 121 |
| sim_S0_r3 | copying simulator S0, leak 1.0, T1, primary operating point | G1 | leak 1.00 | 0.0012 | 1.379 [0.166, 4.165] | 1.736 [1.051, 2.296] | 256.465 [123.641, 296.662] | 0.037 [0.003, 0.112] | 1.000 [0.995, 1.000] | pass | 2 | 558 |
| sim_S1_r1 | copying simulator S1, leak 0.1, T1, primary operating point | Phillips | leak 0.10 | 0.0048 | 0.275 [0.034, 0.647] | 1.601 [1.013, 2.167] | 262.388 [165.766, 296.658] | 0.084 [0.008, 0.197] | 0.997 [0.979, 1.000] | pass | 1 | 92 |
| sim_S1_r2 | copying simulator S1, leak 0.1, T1, primary operating point | Phillips | leak 0.10 | 0.0191 | 1.618 [0.832, 3.772] | 0.887 [0.304, 1.383] | 220.831 [29.410, 293.949] | 0.034 [0.003, 0.104] | 0.998 [0.968, 1.000] | pass | 2 | 345 |
| sim_S1_r3 | copying simulator S1, leak 0.1, T1, primary operating point | Phillips | leak 0.10 | 0.0057 | 0.352 [0.036, 0.916] | 1.918 [1.384, 2.439] | 269.507 [187.263, 297.743] | 0.047 [0.004, 0.137] | 0.999 [0.994, 1.000] | pass | 1 | 94 |
| sim_S2_r1 | copying simulator S2, leak 0.03, T1, primary operating point | Phillips | leak 0.03 | 0.0016 | 2.600 [1.141, 5.106] | 1.061 [0.192, 1.825] | 224.550 [4.249, 295.310] | 0.085 [0.009, 0.199] | 0.993 [0.639, 1.000] | pass | 2 | 441 |
| sim_S2_r2 | copying simulator S2, leak 0.03, T1, primary operating point | Phillips | leak 0.03 | 0.0017 | 0.723 [0.160, 1.576] | 2.025 [1.422, 2.599] | 264.522 [170.990, 297.524] | 0.083 [0.008, 0.203] | 0.998 [0.986, 1.000] | pass | 1 | 81 |
| sim_S2_r3 | copying simulator S2, leak 0.03, T1, primary operating point | Phillips | leak 0.03 | 0.0129 | 1.241 [0.752, 2.861] | 1.585 [0.285, 2.215] | 260.700 [10.810, 297.021] | 0.075 [0.007, 0.212] | 0.998 [0.857, 1.000] | pass | 2 | 324 |
| slope_rho30_r1 | step-against-slope (no barrier, moderate range) | G1 | omega 0.00, rho 30 | 0.1034 | 2.376 [0.368, 4.674] | 3.951 [2.635, 5.662] | 36.790 [12.431, 85.941] | 0.530 [0.060, 1.163] | 0.983 [0.909, 1.000] | pass | 2 | 498 |
| slope_rho30_r2 | step-against-slope (no barrier, moderate range) | G1 | omega 0.00, rho 30 | 0.1340 | 2.075 [0.239, 4.240] | 3.896 [2.816, 5.617] | 25.023 [9.394, 60.272] | 0.432 [0.042, 1.235] | 0.988 [0.906, 1.000] | pass | 2 | 467 |
| slope_rho30_r3 | step-against-slope (no barrier, moderate range) | G1 | omega 0.00, rho 30 | 0.0164 | 0.996 [0.092, 2.820] | 5.030 [3.899, 6.517] | 43.075 [23.497, 76.407] | 0.275 [0.025, 0.798] | 0.997 [0.973, 1.000] | pass | 2 | 457 |
| slope_rho60_r1 | step-against-slope (no barrier, moderate range) | G1 | omega 0.00, rho 60 | 0.0794 | 0.878 [0.095, 2.260] | 4.057 [2.912, 5.391] | 130.899 [59.604, 251.701] | 1.016 [0.687, 1.313] | 0.942 [0.868, 0.978] | pass | 2 | 461 |
| slope_rho60_r2 | step-against-slope (no barrier, moderate range) | G1 | omega 0.00, rho 60 | 0.1290 | 0.757 [0.071, 2.253] | 4.869 [3.786, 6.308] | 73.044 [41.768, 132.812] | 0.340 [0.034, 0.826] | 0.995 [0.970, 1.000] | pass | 2 | 486 |
| slope_rho60_r3 | step-against-slope (no barrier, moderate range) | G1 | omega 0.00, rho 60 | 0.2860 | 0.826 [0.078, 2.623] | 4.290 [3.012, 5.748] | 90.649 [39.172, 193.274] | 0.711 [0.133, 1.150] | 0.974 [0.916, 0.999] | pass | 2 | 456 |
