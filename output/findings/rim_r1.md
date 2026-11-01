# Rim attributes R1: readings beside decoration on the 19 rim sites

Produced by `analyses/rims/rim_r1.py --readings` (git 31b41830fcb9) from the cell CSVs and the rim and decorated tables. Three kinds of number appear and are named: the observed count-based moment F (binomial sampling variance removed; `mls_emergence.signatures.variance.moment_f`); the generator band of the same statistic (10th to 90th percentile of 50 records, NumPy linear rule, endpoints inclusive, eligible at 45 finite values), with status below, inside, above or unavailable, of which only 'above' is read as more variation than the baseline; and posterior summaries of latent site frequencies under Beta(k + 1/2, N - k + 1/2) per site and Dirichlet(counts + 1/2) per composition (4,000 direct draws, two seeds, the Monte Carlo error of each median shown), with Beta(1, 1) and Dirichlet(1) sensitivities. The yardstick is the fixed-target copying simulator with all 28 populations transmitting and the 19 rim sites observed, at the frozen step 0 time scales, two windows (T1: ends over the last L_d; 8LR: the last 8 L_R). Everything is conditional on the 19-site sample (absent: Commerce, Cramor_Place, Cummins, Grant, Holden_Lake, Hollywood, Lake_Cormorant, Nickel, Turnbow), on that yardstick and on that population model; Woodlyn (27 rims) is below Mainfort's 45-rim minimum and Mound Place (48) just above it, so every statistic is also given on the 18 sites without Woodlyn. Nothing here is phrased as a choice among causes (spec section 1). Descriptive summaries under rule 18.

## A. Per-attribute F against the band, beside the common decorated classes on the same sites

### primary point, window T1

| trait | observed F (19) | band (19) | status | observed F (18) | band (18) | status (18) |
|---|---|---|---|---|---|---|
| interior bevel | 0.105 | [-0.001, 0.025] n 50 | **above** | 0.083 | [-0.001, 0.028] n 50 | above |
| exterior notching | 0.069 | [-0.002, 0.024] n 50 | **above** | 0.075 | [-0.001, 0.022] n 50 | above |
| lip flat | 0.031 | [-0.001, 0.027] n 50 | **above** | 0.031 | [-0.001, 0.026] n 50 | above |
| paste Bell (of N) (read separately) | 0.266 | [-0.001, 0.028] n 50 | **above** | 0.243 | [-0.001, 0.022] n 50 | above |
| Parkin_Punctated, current counts | 0.107 | [0.002, 0.027] n 50 | **above** | 0.107 | [0.002, 0.022] n 50 | above |
| Barton/Kent/MPI, current counts | 0.044 | [0.002, 0.029] n 50 | **above** | 0.038 | [0.002, 0.026] n 50 | above |
| Painted, current counts | 0.208 | [0.003, 0.026] n 50 | **above** | 0.211 | [0.002, 0.026] n 50 | above |
| Parkin_Punctated, 1940 PFG parts | 0.109 | [0.001, 0.031] n 50 | **above** | 0.110 | [0.001, 0.028] n 50 | above |
| Barton/Kent/MPI, 1940 PFG parts | 0.062 | [0.002, 0.030] n 50 | **above** | 0.056 | [0.001, 0.029] n 50 | above |
| Painted, 1940 PFG parts | 0.208 | [0.002, 0.018] n 50 | **above** | 0.211 | [0.002, 0.017] n 50 | above |

### primary point, window 8LR

| trait | observed F (19) | band (19) | status | observed F (18) | band (18) | status (18) |
|---|---|---|---|---|---|---|
| interior bevel | 0.105 | [0.004, 0.031] n 50 | **above** | 0.083 | [0.003, 0.029] n 50 | above |
| exterior notching | 0.069 | [0.007, 0.027] n 50 | **above** | 0.075 | [0.005, 0.023] n 50 | above |
| lip flat | 0.031 | [0.004, 0.027] n 50 | **above** | 0.031 | [0.005, 0.027] n 50 | above |
| paste Bell (of N) (read separately) | 0.266 | [0.005, 0.030] n 50 | **above** | 0.243 | [0.005, 0.026] n 50 | above |
| Parkin_Punctated, current counts | 0.107 | [0.005, 0.026] n 50 | **above** | 0.107 | [0.006, 0.025] n 50 | above |
| Barton/Kent/MPI, current counts | 0.044 | [0.006, 0.023] n 50 | **above** | 0.038 | [0.006, 0.020] n 50 | above |
| Painted, current counts | 0.208 | [0.005, 0.024] n 50 | **above** | 0.211 | [0.005, 0.024] n 50 | above |
| Parkin_Punctated, 1940 PFG parts | 0.109 | [0.005, 0.029] n 50 | **above** | 0.110 | [0.006, 0.029] n 50 | above |
| Barton/Kent/MPI, 1940 PFG parts | 0.062 | [0.004, 0.021] n 50 | **above** | 0.056 | [0.004, 0.021] n 50 | above |
| Painted, 1940 PFG parts | 0.208 | [0.006, 0.027] n 50 | **above** | 0.211 | [0.005, 0.025] n 50 | above |

### secondary point, window T1

| trait | observed F (19) | band (19) | status | observed F (18) | band (18) | status (18) |
|---|---|---|---|---|---|---|
| interior bevel | 0.105 | [0.004, 0.025] n 50 | **above** | 0.083 | [0.004, 0.021] n 50 | above |
| exterior notching | 0.069 | [0.002, 0.023] n 50 | **above** | 0.075 | [0.003, 0.023] n 50 | above |
| lip flat | 0.031 | [0.002, 0.031] n 50 | **above** | 0.031 | [0.003, 0.028] n 50 | above |
| paste Bell (of N) (read separately) | 0.266 | [0.003, 0.022] n 50 | **above** | 0.243 | [0.001, 0.023] n 50 | above |
| Parkin_Punctated, current counts | 0.107 | [0.004, 0.025] n 50 | **above** | 0.107 | [0.003, 0.026] n 50 | above |
| Barton/Kent/MPI, current counts | 0.044 | [0.005, 0.023] n 50 | **above** | 0.038 | [0.004, 0.023] n 50 | above |
| Painted, current counts | 0.208 | [0.004, 0.018] n 50 | **above** | 0.211 | [0.004, 0.019] n 50 | above |
| Parkin_Punctated, 1940 PFG parts | 0.109 | [0.004, 0.023] n 50 | **above** | 0.110 | [0.004, 0.022] n 50 | above |
| Barton/Kent/MPI, 1940 PFG parts | 0.062 | [0.004, 0.019] n 50 | **above** | 0.056 | [0.004, 0.020] n 50 | above |
| Painted, 1940 PFG parts | 0.208 | [0.004, 0.021] n 50 | **above** | 0.211 | [0.004, 0.022] n 50 | above |

### secondary point, window 8LR

| trait | observed F (19) | band (19) | status | observed F (18) | band (18) | status (18) |
|---|---|---|---|---|---|---|
| interior bevel | 0.105 | [0.007, 0.030] n 50 | **above** | 0.083 | [0.006, 0.032] n 50 | above |
| exterior notching | 0.069 | [0.005, 0.029] n 50 | **above** | 0.075 | [0.005, 0.029] n 50 | above |
| lip flat | 0.031 | [0.003, 0.029] n 50 | **above** | 0.031 | [0.004, 0.029] n 50 | above |
| paste Bell (of N) (read separately) | 0.266 | [0.005, 0.033] n 50 | **above** | 0.243 | [0.005, 0.035] n 50 | above |
| Parkin_Punctated, current counts | 0.107 | [0.007, 0.026] n 50 | **above** | 0.107 | [0.007, 0.026] n 50 | above |
| Barton/Kent/MPI, current counts | 0.044 | [0.006, 0.025] n 50 | **above** | 0.038 | [0.007, 0.025] n 50 | above |
| Painted, current counts | 0.208 | [0.008, 0.025] n 50 | **above** | 0.211 | [0.007, 0.027] n 50 | above |
| Parkin_Punctated, 1940 PFG parts | 0.109 | [0.009, 0.024] n 50 | **above** | 0.110 | [0.006, 0.025] n 50 | above |
| Barton/Kent/MPI, 1940 PFG parts | 0.062 | [0.009, 0.027] n 50 | **above** | 0.056 | [0.010, 0.028] n 50 | above |
| Painted, 1940 PFG parts | 0.208 | [0.010, 0.028] n 50 | **above** | 0.211 | [0.010, 0.029] n 50 | above |

Lip form is a reconstructed field: the file's flat and round percentages sum to 99, 101 and 98 at Beck, Irby, Starkley; over the 27 combinations of one rim more or fewer at those sites the lip-flat F ranges 0.030 to 0.032.

Smallest margin of the observed F to the band edge it was read against, over the points, windows and both site sets (observed minus the upper edge when above; four decimals, since a status at the third decimal can rest on the fourth):

- interior bevel: +0.0512 at the secondary point, window 8LR, 18 sites (observed 0.0828, band [0.0058, 0.0316]); 8 of 8 cells above.
- exterior notching: +0.0395 at the secondary point, window 8LR, 19 sites (observed 0.0686, band [0.0050, 0.0291]); 8 of 8 cells above.
- lip flat: +0.0001 at the secondary point, window T1, 19 sites (observed 0.0309, band [0.0020, 0.0308]); 8 of 8 cells above; the lip reconstruction range 0.0302 to 0.0317 contains that band edge (0.0308), so the status at this cell is not resolved by the record.
- paste Bell (of N): +0.2077 at the secondary point, window 8LR, 18 sites (observed 0.2427, band [0.0048, 0.0350]); 8 of 8 cells above.

## B. Association with geography

Spearman of |p_i - p_j| with straight-line distance over the 171 pairs (positive: pairs farther apart differ more), and of p with latitude (positive: higher to the north) and with longitude (positive: higher to the east). Observed value; posterior median [5%, 95%]; the Beta(1, 1) or Dirichlet(1) posterior median in the last column. Rows for the 18 sites give observed values only.

| trait | rho pair difference vs distance | rho vs latitude | rho vs longitude | alpha 1 medians (lat, lon, dist) |
|---|---|---|---|---|
| interior bevel | 0.22; posterior 0.21 [0.11, 0.31] (MC error of the median 0.001) | -0.61; posterior -0.58 [-0.68, -0.47] (MC error of the median 0.004) | 0.33; posterior 0.30 [0.18, 0.42] (MC error of the median 0.002) | -0.57, +0.30, +0.20 |
| interior bevel, 18 sites | 0.23 | -0.62 | 0.23 | |
| exterior notching | 0.16; posterior 0.12 [0.04, 0.22] (MC error of the median 0.001) | -0.48; posterior -0.39 [-0.53, -0.25] (MC error of the median 0.002) | 0.51; posterior 0.51 [0.37, 0.65] (MC error of the median 0.004) | -0.38, +0.52, +0.12 |
| exterior notching, 18 sites | 0.14 | -0.46 | 0.44 | |
| lip flat | -0.09; posterior -0.05 [-0.12, 0.05] (MC error of the median 0.002) | -0.23; posterior -0.20 [-0.38, -0.01] (MC error of the median 0.001) | 0.39; posterior 0.39 [0.21, 0.57] (MC error of the median 0.006) | -0.19, +0.39, -0.05 |
| lip flat, 18 sites | -0.04 | -0.23 | 0.52 | |
| paste Bell (of N) | 0.31; posterior 0.28 [0.20, 0.36] (MC error of the median 0.001) | -0.41; posterior -0.39 [-0.49, -0.27] (MC error of the median 0.000) | 0.61; posterior 0.61 [0.53, 0.71] (MC error of the median 0.002) | -0.39, +0.61, +0.28 |
| paste Bell (of N), 18 sites | 0.27 | -0.42 | 0.54 | |
| Parkin_Punctated, current counts | 0.29; posterior 0.29 [0.25, 0.33] (MC error of the median 0.001) | 0.66; posterior 0.65 [0.58, 0.72] (MC error of the median 0.002) | 0.18; posterior 0.15 [0.06, 0.25] (MC error of the median 0.002) | +0.67, +0.13, +0.29 |
| Barton/Kent/MPI, current counts | 0.05; posterior 0.07 [-0.01, 0.13] (MC error of the median 0.002) | 0.08; posterior 0.16 [0.02, 0.29] (MC error of the median 0.002) | -0.31; posterior -0.26 [-0.36, -0.16] (MC error of the median 0.002) | +0.16, -0.26, +0.07 |
| Painted, current counts | 0.35; posterior 0.33 [0.28, 0.38] (MC error of the median 0.001) | -0.75; posterior -0.74 [-0.78, -0.70] (MC error of the median 0.000) | -0.28; posterior -0.29 [-0.35, -0.24] (MC error of the median 0.004) | -0.74, -0.29, +0.33 |
| Parkin_Punctated, 1940 PFG parts | 0.32; posterior 0.32 [0.28, 0.37] (MC error of the median 0.001) | 0.64; posterior 0.66 [0.57, 0.73] (MC error of the median 0.001) | 0.15; posterior 0.12 [0.01, 0.24] (MC error of the median 0.002) | +0.67, +0.10, +0.33 |
| Barton/Kent/MPI, 1940 PFG parts | 0.04; posterior 0.05 [-0.03, 0.12] (MC error of the median 0.000) | 0.08; posterior 0.13 [0.01, 0.28] (MC error of the median 0.005) | -0.26; posterior -0.22 [-0.32, -0.11] (MC error of the median 0.002) | +0.14, -0.22, +0.05 |
| Painted, 1940 PFG parts | 0.33; posterior 0.31 [0.26, 0.38] (MC error of the median 0.002) | -0.72; posterior -0.72 [-0.77, -0.63] (MC error of the median 0.000) | -0.28; posterior -0.30 [-0.36, -0.24] (MC error of the median 0.000) | -0.72, -0.30, +0.32 |

## C. Rim and decoration pair differences after linear adjustment on ranks for distance

Rim difference |p_i - p_j|; decoration difference the total variation 0.5 sum_c |p_ic - p_jc| over the ten classes (zero for identical compositions). The statistic is the Spearman correlation of the two after least-squares removal of ranked distance from each: association after linear adjustment on ranks, not 'at a given distance' (a nonlinear dependence on distance survives it). Positive: pairs that differ more in the attribute also differ more in decoration. Each rim attribute separately; no combined rim measure, since the marginal table carries no co-occurrence and the two collections' physical overlap is unknown.

| trait | against current counts: plain, adjusted | against 1940 PFG parts: plain, adjusted | alpha 1 adjusted medians |
|---|---|---|---|
| interior bevel | 0.17; posterior 0.14 [0.05, 0.23] (MC error of the median 0.001); adjusted 0.09; posterior 0.06 [-0.04, 0.16] (MC error of the median 0.002) | 0.25; posterior 0.20 [0.10, 0.31] (MC error of the median 0.002); adjusted 0.18; posterior 0.14 [0.02, 0.26] (MC error of the median 0.002) | +0.05, +0.13 |
| exterior notching | 0.10; posterior 0.07 [-0.01, 0.15] (MC error of the median 0.000); adjusted 0.04; posterior 0.02 [-0.06, 0.10] (MC error of the median 0.001) | 0.15; posterior 0.12 [0.03, 0.21] (MC error of the median 0.001); adjusted 0.09; posterior 0.08 [-0.01, 0.17] (MC error of the median 0.000) | +0.01, +0.07 |
| lip flat | -0.09; posterior -0.08 [-0.15, 0.00] (MC error of the median 0.001); adjusted -0.06; posterior -0.07 [-0.15, 0.02] (MC error of the median 0.000) | -0.06; posterior -0.06 [-0.14, 0.03] (MC error of the median 0.001); adjusted -0.03; posterior -0.04 [-0.13, 0.05] (MC error of the median 0.002) | -0.07, -0.05 |
| paste Bell (of N) | 0.25; posterior 0.21 [0.12, 0.30] (MC error of the median 0.000); adjusted 0.14; posterior 0.11 [0.02, 0.21] (MC error of the median 0.002) | 0.26; posterior 0.22 [0.13, 0.31] (MC error of the median 0.000); adjusted 0.15; posterior 0.12 [0.02, 0.22] (MC error of the median 0.002) | +0.11, +0.12 |

## D. Paste, separately

Paste is a type-level distinction of the 1951 typology (Bell Plain against Neeley's Ferry Plain) and is not independent of how the phases were drawn; the file's paste percentages leave 550 of 2065 rims unresolved, which are not negative Bell observations. The tables above read Bell of N; here Bell is also given with the named pastes as denominator (two denominators). Bell of N: mean site frequency 0.283; Bell of named: 0.371.

## E. Published attributes reported but not read

| attribute | rims | sites with any (of 19) | pooled frequency |
|---|---|---|---|
| beaded | 28 | 10 | 0.0136 |
| noded | 13 | 7 | 0.0063 |
| wide_nodes | 10 | 8 | 0.0048 |
| notched_top | 23 | 13 | 0.0111 |

A recorded scope choice (at least 100 rims and 12 sites with any), not a demonstrated resolution limit. Bevel degrees as shares of the beveled rims: sharp 90, moderate 103, slight 55 of 248.

## F. The priors (rule 20)

Beta(1/2, 1/2) on each site frequency: one pseudo-observation in total, weight 1/(N + 1) on the posterior mean (0.036 at 27 rims, 0.0036 at 277); direction toward 1/2. Observable prior predictive, Beta-Binomial(N, 1/2, 1/2) at each site's N: 5th, 50th and 95th percentiles of the predicted count against the observed count for the three read attributes:

| site | N | prior predictive count (5%, 50%, 95%) | bevel | notching | lip flat |
|---|---|---|---|---|---|
| Barton_Ranch | 102 | 0, 51, 102 | 2 | 6 | 35 |
| Beck | 183 | 1, 92, 182 | 42 | 22 | 79 |
| Belle_Meade | 143 | 1, 72, 142 | 44 | 29 | 74 |
| Big_Eddy | 88 | 0, 44, 88 | 14 | 1 | 20 |
| Castile_Landing | 72 | 0, 36, 72 | 1 | 0 | 12 |
| Clay_Hill | 159 | 1, 79, 158 | 19 | 10 | 59 |
| Davis | 96 | 0, 48, 96 | 2 | 2 | 31 |
| Fortune | 84 | 0, 42, 84 | 0 | 0 | 16 |
| Irby | 277 | 1, 139, 276 | 48 | 57 | 101 |
| Kent_Place | 153 | 1, 77, 152 | 14 | 4 | 54 |
| Mound_Place | 48 | 0, 24, 48 | 4 | 0 | 17 |
| Neeleys_Ferry | 108 | 0, 54, 108 | 1 | 1 | 23 |
| Parkin | 67 | 0, 33, 67 | 5 | 4 | 30 |
| Rose_Mound | 76 | 0, 38, 76 | 2 | 0 | 24 |
| Starkley | 108 | 0, 54, 108 | 21 | 7 | 26 |
| Vernon_Paul | 99 | 0, 50, 99 | 3 | 5 | 37 |
| Walls | 107 | 0, 53, 107 | 14 | 8 | 45 |
| Williamson | 68 | 0, 34, 68 | 2 | 0 | 19 |
| Woodlyn | 27 | 0, 13, 27 | 10 | 4 | 6 |

The 5th to 95th percentile interval of the predicted count spans 99 to 100 percent of N across the sites and contains 57 of the 57 observed counts: the prior concentrates the predicted count nowhere, so it does not already resemble the data. The Beta(1, 1) sensitivity medians are in the last columns of sections B and C.

Recovery program: 20 replicates at three true configurations on these 19 sites' totals, 1,000 draws per fit, 90 percent intervals; a coverage proportion at 20 replicates has binomial standard deviation 0.067 at a true rate of 0.9, so the table resolves gross failure, not the nominal rate.

| configuration | true rho latitude | coverage of rho latitude | coverage of rho pair distance | sign recovered |
|---|---|---|---|---|
| positive | +0.75 | 0.85 | 0.90 | 1.00 |
| opposite | -0.69 | 0.85 | 0.95 | 1.00 |
| zero | +0.09 | 0.90 | 0.80 | 0.90 |

## G. Limits stated in advance

19 sites of 28; rims from collections Mainfort names only in aggregate, never per site and never as the 1940 survey collections, so whether a site's rims and its decorated counts share material is unknown; one observer, which helps consistency but does not make recording independent across sites; totals 27 to 277, so per-site frequencies carry wide intervals and the F estimator's own variance at 19 sites is large; the yardstick is the fixed-target baseline (no regional drift) under one population model; interior beveling is spatially as well as temporally patterned by Mainfort's own account (p. 37), so no clock reading is made.
