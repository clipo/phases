# Per-fit sampler diagnostics for the latent-group probe's Stan fits (rule 16)

Produced by `analyses/probe_latent_groups/export_fit_diagnostics.py` from the fit records in `output/gp_cache/` (gitignored; 297 records with diagnostics). Every value is the one the gate read at the final rung; the gate thresholds are R-hat at most 1.01 on gated variables, bulk and tail ESS at least 400, no divergent transitions, E-BFMI at least 0.3 and, for A1, a per-chain co-membership split of at most 0.5. Treedepth saturation is reported, not gated. Rung settings are iterations per chain / warmup, adapt_delta and max_treedepth; four chains throughout. A label with two rows was fitted under two versions of the Stan source, with the class-intercept prior hard-coded as Normal(0, 2) and, after D25 (2026-10-18), read from the data (default Normal(0, 10)); the findings tables named per family report the data-driven version under the plain name and the hard-coded version as `*_alphaN0-2.md`. Descriptive export; exploratory design validation, not evidence about the record. Columns in `fit_diagnostics.csv`.

## A0 spatial field: basin fit, prior predictive, sensitivity, recovery, intercept priors (refit_a0_basin.R, a0_intercept_prior.R)

23 records, 22 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a0_alpha_P10 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0059 (lp__) | 874 (tau[1]) | 1666 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.672 |  | pass | 146 |
| a0_alpha_P2 | Normal(0, 2) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0034 (lp__) | 608 (tau[1]) | 1249 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.676 |  | pass | 89 |
| a0_alpha_PC | Normal(pooled log-odds, 2) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0050 (lp__) | 1053 (tau[1]) | 1928 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.720 |  | pass | 53 |
| a0_basin | Normal(0, 10) | 2026-11-03 | 1 (3000/1500, 0.95, 12) | 1.0037 (tau[1]) | 628 (tau[1]) | 1535 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.675 |  | pass | 125 |
| a0_basin | Normal(0, 2), hard-coded | 2026-10-16 | 1 (3000/1500, 0.95, 12) | 1.0034 (lp__) | 608 (tau[1]) | 1249 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.676 |  | pass | 87 |
| a0_prior_predictive | Normal(0, 10) | 2026-11-03 | 1 (2000/1000, 0.9, 10) | 1.0056 (alpha[1]) | 1469 (lp__) | 1831 (tau[1]) | 0/4000 (0.00%) | 0/4000 (0.0%; deepest 5) | 0.963 |  | pass | 4 |
| a0_prior_predictive | Normal(0, 2), hard-coded | 2026-10-16 | 1 (2000/1000, 0.9, 10) | 1.0042 (alpha[4]) | 1740 (lp__) | 1783 (sigma) | 0/4000 (0.00%) | 0/4000 (0.0%; deepest 5) | 0.958 |  | pass | 4 |
| a0_recovery_share0.2_rho5 | Normal(0, 10) | 2026-11-03 | 2 (6000/3000, 0.99, 13) | 1.0037 (f[3,3]) | 1127 (alpha[2]) | 1467 (alpha[2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.667 |  | pass | 247 |
| a0_recovery_share0.2_rho5 | Normal(0, 10) | 2026-11-03 | 2 (6000/3000, 0.99, 13) | 1.0037 (f[3,3]) | 1127 (alpha[2]) | 1467 (alpha[2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.667 |  | pass | 246 |
| a0_recovery_share0.2_rho5 | Normal(0, 2), hard-coded | 2026-10-16 | 2 (6000/3000, 0.99, 13) | 1.0076 (alpha[2]) | 1142 (alpha[2]) | 1919 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.657 |  | pass | 179 |
| a0_recovery_share0.2_rho60 | Normal(0, 10) | 2026-11-03 | 2 (6000/3000, 0.99, 13) | 1.0062 (alpha[3]) | 905 (rho) | 1046 (alpha[2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.665 |  | pass | 222 |
| a0_recovery_share0.2_rho60 | Normal(0, 10) | 2026-11-03 | 2 (6000/3000, 0.99, 13) | 1.0062 (alpha[3]) | 905 (rho) | 1046 (alpha[2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.665 |  | pass | 221 |
| a0_recovery_share0.2_rho60 | Normal(0, 2), hard-coded | 2026-10-16 | 2 (6000/3000, 0.99, 13) | 1.0039 (rho) | 1125 (alpha[2]) | 1684 (alpha[6]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.656 |  | pass | 181 |
| a0_recovery_share0.9_rho5 | Normal(0, 10) | 2026-11-03 | 2 (6000/3000, 0.99, 13) | 1.0087 (rho) | 491 (rho) | 205 (rho) | 11/12000 (0.09%) | 0/12000 (0.0%; deepest 11) | 0.646 |  | fail: ess_tail, divergences | 115 |
| a0_recovery_share0.9_rho5 | Normal(0, 10) | 2026-11-03 | 3 (12000/6000, 0.999, 14) | 1.0028 (alpha[3]) | 1323 (rho) | 596 (rho) | 0/24000 (0.00%) | 0/24000 (0.0%; deepest 13) | 0.664 |  | pass | 657 |
| a0_recovery_share0.9_rho5 | Normal(0, 2), hard-coded | 2026-10-16 | 2 (6000/3000, 0.99, 13) | 1.0041 (alpha[1]) | 1233 (alpha[2]) | 2041 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.617 |  | pass | 94 |
| a0_recovery_share0.9_rho60 | Normal(0, 10) | 2026-11-03 | 1 (3000/1500, 0.95, 12) | 1.0024 (f[5,6]) | 1075 (rho) | 1410 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.705 |  | pass | 64 |
| a0_recovery_share0.9_rho60 | Normal(0, 10) | 2026-11-03 | 1 (3000/1500, 0.95, 12) | 1.0024 (f[5,6]) | 1075 (rho) | 1410 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.705 |  | pass | 66 |
| a0_recovery_share0.9_rho60 | Normal(0, 2), hard-coded | 2026-10-16 | 1 (3000/1500, 0.95, 12) | 1.0036 (f[8,3]) | 873 (rho) | 1430 (f[28,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.650 |  | pass | 44 |
| a0_sens_rho_uniform | Normal(0, 10) | 2026-11-03 | 2 (6000/3000, 0.99, 13) | 1.0028 (f[4,1]) | 1349 (tau[1]) | 2667 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.643 |  | pass | 702 |
| a0_sens_rho_uniform | Normal(0, 2), hard-coded | 2026-10-16 | 2 (6000/3000, 0.99, 13) | 1.0045 (tau[1]) | 1215 (tau[1]) | 2236 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.660 |  | pass | 263 |
| a0_sens_scale05 | Normal(0, 10) | 2026-11-03 | 1 (3000/1500, 0.95, 12) | 1.0043 (rho) | 905 (f[3,1]) | 1748 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.715 |  | pass | 52 |
| a0_sens_scale05 | Normal(0, 2), hard-coded | 2026-10-16 | 1 (3000/1500, 0.95, 12) | 1.0024 (lp__) | 919 (tau[1]) | 1708 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.696 |  | pass | 53 |

## A1 latent-unit recovery and equivalence fits (validate_a1.R; a1_validation.md)

17 records, 4 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a1val_a0_on_zero | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0041 (f[5,3]) | 960 (tau[1]) | 1729 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.746 |  | pass | 58 |
| a1val_a1_fixed0_on_zero | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0056 (lp__) | 1069 (tau[1]) | 1940 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.763 |  | pass | 165 |
| a1val_a1_on_zero | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0244 (comembership[27,25]) | 277 (comembership[25,21]) | 450 (comembership[14,11]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.632 | 0.144 | fail: rhat, ess_bulk | 633 |
| a1val_sg0.25_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0056 (comembership[13,2]) | 722 (comembership[3,2]) | 1073 (comembership[5,3]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.706 | 0.074 | pass | 544 |
| a1val_sg0.25_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0282 (comembership[14,2]) | 260 (comembership[13,8]) | 244 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.673 | 0.112 | fail: rhat, ess_bulk, ess_tail | 633 |
| a1val_sg0.25_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0687 (comembership[5,2]) | 84 (comembership[5,4]) | 241 (comembership[20,15]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.756 | 0.247 | fail: rhat, ess_bulk, ess_tail | 614 |
| a1val_sg0.25_r4 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0258 (comembership[21,3]) | 250 (comembership[21,3]) | 427 (comembership[10,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.658 | 0.163 | fail: rhat, ess_bulk | 694 |
| a1val_sg0.25_r5 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0204 (comembership[19,2]) | 200 (comembership[19,2]) | 429 (comembership[28,25]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.680 | 0.165 | fail: rhat, ess_bulk | 636 |
| a1val_sg0.5_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0097 (comembership[20,2]) | 737 (comembership[20,3]) | 888 (comembership[26,13]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 11) | 0.669 | 0.117 | fail: divergences | 684 |
| a1val_sg0.5_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0280 (comembership[22,13]) | 240 (comembership[13,4]) | 84 (comembership[13,4]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.675 | 0.141 | fail: rhat, ess_bulk, ess_tail | 649 |
| a1val_sg0.5_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0183 (comembership[28,2]) | 487 (comembership[5,2]) | 543 (comembership[22,15]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.696 | 0.148 | fail: rhat | 625 |
| a1val_sg0.5_r4 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0246 (comembership[8,2]) | 223 (comembership[25,21]) | 287 (comembership[23,11]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.701 | 0.137 | fail: rhat, ess_bulk, ess_tail | 656 |
| a1val_sg0.5_r5 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0076 (comembership[5,2]) | 700 (comembership[20,19]) | 632 (comembership[20,15]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.659 | 0.077 | pass | 670 |
| a1val_sg1_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0437 (comembership[22,8]) | 126 (comembership[25,22]) | 422 (comembership[19,12]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.653 | 0.222 | fail: rhat, ess_bulk | 661 |
| a1val_sg1_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0236 (comembership[21,11]) | 237 (comembership[24,21]) | 217 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.669 | 0.178 | fail: rhat, ess_bulk, ess_tail | 1106 |
| a1val_sg2_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.2900 (comembership[18,2]) | 11 (comembership[18,2]) | 11 (comembership[24,14]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.733 | 0.452 | fail: rhat, ess_bulk, ess_tail | 650 |
| a1val_sg2_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0357 (comembership[3,2]) | 280 (comembership[3,2]) | 424 (comembership[24,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.667 | 0.237 | fail: rhat, ess_bulk | 506 |

## A1 refits under the Normal(0, 10) intercept prior (a1_prior_check.R; a1_prior_check.md)

3 records, 0 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a1prior_sg0.5_r1 | Normal(0, 10) | 2026-10-19 | 2 (6000/3000, 0.99, 13) | 1.0128 (comembership[22,21]) | 656 (comembership[20,3]) | 922 (comembership[26,3]) | 0/12000 (0.00%) | 1/12000 (0.0%; deepest 13) | 0.677 | 0.101 | fail: rhat | 994 |
| a1prior_sg0.5_r5 | Normal(0, 10) | 2026-10-19 | 2 (6000/3000, 0.99, 13) | 1.0086 (comembership[22,20]) | 575 (comembership[27,25]) | 372 (comembership[20,10]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.701 | 0.093 | fail: ess_tail | 793 |
| a1prior_zero | Normal(0, 10) | 2026-10-19 | 2 (6000/3000, 0.99, 13) | 1.0297 (comembership[22,19]) | 231 (comembership[21,19]) | 423 (comembership[9,7]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.697 | 0.143 | fail: rhat, ess_bulk | 773 |

## A3 barrier recovery and equivalence fits, range free (validate_a3.R; a3_validation.md and a3_validation_alphaN0-2.md)

82 records, 82 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3val_a0_on_om0 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0040 (alpha[5]) | 1135 (rho) | 1708 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.742 |  | pass | 101 |
| a3val_a0_on_om0 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0054 (f[15,1]) | 1058 (tau[1]) | 1712 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.762 |  | pass | 76 |
| a3val_a3_fixed0_on_om0 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0053 (rho) | 973 (rho) | 1500 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.728 |  | pass | 117 |
| a3val_a3_fixed0_on_om0 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0038 (lp__) | 1179 (tau[1]) | 2007 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.725 |  | pass | 102 |
| a3val_own_G1_om0.25_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0066 (sigma) | 514 (rho) | 1036 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.735 |  | pass | 110 |
| a3val_own_G1_om0.25_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0076 (tau[1]) | 768 (tau[1]) | 1022 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.737 |  | pass | 85 |
| a3val_own_G1_om0.25_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0030 (tau[1]) | 1352 (lp__) | 2105 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.723 |  | pass | 136 |
| a3val_own_G1_om0.25_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0025 (f[1,2,2]) | 835 (tau[1]) | 1472 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 109 |
| a3val_own_G1_om0.25_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0030 (rho) | 1699 (rho) | 2816 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.742 |  | pass | 412 |
| a3val_own_G1_om0.25_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0031 (tau[1]) | 720 (tau[1]) | 955 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.755 |  | pass | 105 |
| a3val_own_G1_om0.5_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0062 (sigma) | 937 (sigma) | 622 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.748 |  | pass | 133 |
| a3val_own_G1_om0.5_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0090 (omega) | 478 (omega) | 739 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.704 |  | pass | 91 |
| a3val_own_G1_om0.5_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0033 (f[1,23,7]) | 1077 (tau[1]) | 807 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.712 |  | pass | 114 |
| a3val_own_G1_om0.5_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0051 (tau[1]) | 851 (omega) | 1182 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.716 |  | pass | 72 |
| a3val_own_G1_om0.5_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0038 (tau[1]) | 1202 (tau[1]) | 1781 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.732 |  | pass | 135 |
| a3val_own_G1_om0.5_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0059 (tau[1]) | 604 (tau[1]) | 598 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.738 |  | pass | 108 |
| a3val_own_G1_om0_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0045 (rho) | 739 (omega) | 725 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.760 |  | pass | 114 |
| a3val_own_G1_om0_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0060 (omega) | 1301 (omega) | 2547 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.751 |  | pass | 387 |
| a3val_own_G1_om0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0025 (lp__) | 1340 (lp__) | 1799 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.726 |  | pass | 138 |
| a3val_own_G1_om0_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0028 (tau[1]) | 1116 (tau[1]) | 1414 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.717 |  | pass | 96 |
| a3val_own_G1_om0_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0056 (lp__) | 628 (rho) | 962 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.715 |  | pass | 73 |
| a3val_own_G1_om0_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0018 (sigma) | 1549 (tau[1]) | 1878 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.697 |  | pass | 399 |
| a3val_own_G1_om1_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0076 (rho) | 634 (rho) | 1042 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.749 |  | pass | 82 |
| a3val_own_G1_om1_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0025 (f[1,8,8]) | 1168 (tau[1]) | 1189 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.761 |  | pass | 105 |
| a3val_own_G1_om1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0048 (f[1,28,1]) | 1291 (rho) | 1495 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.732 |  | pass | 129 |
| a3val_own_G1_om1_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0059 (tau[1]) | 658 (tau[1]) | 1311 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.705 |  | pass | 78 |
| a3val_own_G1_om1_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0033 (tau[1]) | 1323 (tau[1]) | 1486 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.749 |  | pass | 127 |
| a3val_own_G1_om1_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0082 (omega) | 706 (omega) | 1218 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.813 |  | pass | 86 |
| a3val_own_L_om0.25_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0045 (rho) | 1301 (tau[1]) | 1508 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.770 |  | pass | 122 |
| a3val_own_L_om0.25_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0027 (omega) | 973 (tau[1]) | 1091 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.702 |  | pass | 108 |
| a3val_own_L_om0.25_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0045 (tau[1]) | 1006 (rho) | 1619 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.705 |  | pass | 130 |
| a3val_own_L_om0.25_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0066 (tau[1]) | 637 (tau[1]) | 742 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.729 |  | pass | 104 |
| a3val_own_L_om0.25_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0044 (sigma) | 2033 (sigma) | 2767 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.724 |  | pass | 253 |
| a3val_own_L_om0.25_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0078 (lp__) | 463 (tau[1]) | 887 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.682 |  | pass | 80 |
| a3val_own_L_om0.5_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0063 (sigma) | 706 (rho) | 1557 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.799 |  | pass | 115 |
| a3val_own_L_om0.5_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0027 (tau[1]) | 1358 (sigma) | 1042 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.727 |  | pass | 449 |
| a3val_own_L_om0.5_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0016 (f[1,22,3]) | 2258 (rho) | 3756 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.753 |  | pass | 445 |
| a3val_own_L_om0.5_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0087 (lp__) | 1025 (tau[1]) | 1663 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.720 |  | pass | 101 |
| a3val_own_L_om0.5_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0036 (tau[1]) | 1239 (rho) | 1870 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.733 |  | pass | 131 |
| a3val_own_L_om0.5_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0038 (omega) | 919 (tau[1]) | 1337 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.732 |  | pass | 71 |
| a3val_own_L_om0_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0038 (sigma) | 1208 (omega) | 2385 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.683 |  | pass | 122 |
| a3val_own_L_om0_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0044 (lp__) | 809 (omega) | 1477 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.730 |  | pass | 72 |
| a3val_own_L_om0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0037 (lp__) | 1179 (lp__) | 2031 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.661 |  | pass | 125 |
| a3val_own_L_om0_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0056 (tau[1]) | 1053 (tau[1]) | 1788 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.730 |  | pass | 86 |
| a3val_own_L_om0_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0086 (rho) | 872 (sigma) | 1411 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.754 |  | pass | 106 |
| a3val_own_L_om0_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0086 (omega) | 998 (omega) | 1100 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.782 |  | pass | 87 |
| a3val_own_L_om1_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0076 (omega) | 1286 (omega) | 1345 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.752 |  | pass | 132 |
| a3val_own_L_om1_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0031 (omega) | 650 (omega) | 934 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.761 |  | pass | 72 |
| a3val_own_L_om1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0042 (rho) | 1352 (omega) | 2438 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.744 |  | pass | 130 |
| a3val_own_L_om1_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0059 (omega) | 1059 (omega) | 2269 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.733 |  | pass | 92 |
| a3val_own_L_om1_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0018 (rho) | 1569 (rho) | 3031 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.729 |  | pass | 402 |
| a3val_own_L_om1_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0059 (lp__) | 780 (tau[1]) | 1302 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.737 |  | pass | 113 |
| a3val_sim_S0_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0035 (lp__) | 1508 (lp__) | 1653 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.675 |  | pass | 62 |
| a3val_sim_S0_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0084 (omega) | 910 (omega) | 1202 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.693 |  | pass | 391 |
| a3val_sim_S0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0036 (lp__) | 1417 (lp__) | 2479 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.674 |  | pass | 131 |
| a3val_sim_S0_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0059 (lp__) | 1066 (omega) | 1431 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.686 |  | pass | 88 |
| a3val_sim_S0_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0024 (f[1,3,2]) | 1496 (omega) | 1875 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.753 |  | pass | 74 |
| a3val_sim_S0_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0025 (omega) | 1167 (omega) | 3513 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.683 |  | pass | 403 |
| a3val_sim_S1_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0023 (sigma) | 1075 (sigma) | 1342 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.594 |  | pass | 214 |
| a3val_sim_S1_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0040 (lp__) | 850 (lp__) | 1460 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.573 |  | pass | 67 |
| a3val_sim_S1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0029 (f[1,4,5]) | 1607 (lp__) | 2562 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.775 |  | pass | 63 |
| a3val_sim_S1_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0074 (omega) | 1065 (omega) | 2788 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.707 |  | pass | 248 |
| a3val_sim_S1_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0088 (omega) | 764 (lp__) | 1235 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.590 |  | pass | 68 |
| a3val_sim_S1_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0037 (omega) | 1251 (omega) | 1908 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.670 |  | pass | 71 |
| a3val_sim_S2_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0034 (sigma) | 1009 (sigma) | 1745 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.814 |  | pass | 54 |
| a3val_sim_S2_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0069 (sigma) | 405 (sigma) | 1202 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 100 |
| a3val_sim_S2_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0040 (sigma) | 1299 (sigma) | 2181 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.682 |  | pass | 240 |
| a3val_sim_S2_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0036 (omega) | 1145 (omega) | 1620 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.696 |  | pass | 71 |
| a3val_sim_S2_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0047 (f[1,19,1]) | 896 (sigma) | 1165 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.764 |  | pass | 41 |
| a3val_sim_S2_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0068 (sigma) | 820 (sigma) | 719 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.670 |  | pass | 257 |
| a3val_slope_rho30_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0019 (lp__) | 1433 (rho) | 2472 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.740 |  | pass | 485 |
| a3val_slope_rho30_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0034 (rho) | 964 (rho) | 2302 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.757 |  | pass | 392 |
| a3val_slope_rho30_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0018 (f[1,10,4]) | 1556 (rho) | 2732 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.724 |  | pass | 411 |
| a3val_slope_rho30_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0067 (sigma) | 810 (rho) | 1840 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.774 |  | pass | 397 |
| a3val_slope_rho30_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0030 (lp__) | 1226 (rho) | 2334 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.739 |  | pass | 129 |
| a3val_slope_rho30_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0014 (lp__) | 2222 (rho) | 3468 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.786 |  | pass | 372 |
| a3val_slope_rho60_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0014 (f[1,21,9]) | 1918 (rho) | 3001 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.758 |  | pass | 399 |
| a3val_slope_rho60_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0040 (tau[1]) | 1302 (tau[1]) | 2432 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.713 |  | pass | 267 |
| a3val_slope_rho60_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0025 (log_lik[24]) | 1319 (lp__) | 2123 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.767 |  | pass | 132 |
| a3val_slope_rho60_r2 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0080 (lp__) | 1281 (lp__) | 1871 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.795 |  | pass | 108 |
| a3val_slope_rho60_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0045 (lp__) | 1063 (tau[1]) | 1697 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.739 |  | pass | 111 |
| a3val_slope_rho60_r3 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0035 (tau[1]) | 1228 (tau[1]) | 2374 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.728 |  | pass | 380 |

## A3 with the range pinned at 6 km (a3_pinned_range.R 6; a3_pinned_rho6.md and _alphaN0-2.md)

78 records, 78 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3pin6_own_G1_om0.25_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0087 (omega) | 1070 (omega) | 1694 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.745 |  | pass | 60 |
| a3pin6_own_G1_om0.25_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0039 (f[1,11,3]) | 865 (sigma) | 1208 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.760 |  | pass | 72 |
| a3pin6_own_G1_om0.25_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0046 (f[1,23,1]) | 1226 (f[1,14,1]) | 2394 (f[1,22,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.722 |  | pass | 68 |
| a3pin6_own_G1_om0.25_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0050 (omega) | 866 (sigma) | 1022 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.703 |  | pass | 100 |
| a3pin6_own_G1_om0.25_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0019 (lp__) | 2228 (lp__) | 4114 (f[1,19,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.737 |  | pass | 256 |
| a3pin6_own_G1_om0.25_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0066 (f[1,13,2]) | 1087 (sigma) | 1356 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.749 |  | pass | 73 |
| a3pin6_own_G1_om0.5_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0049 (lp__) | 1422 (lp__) | 1804 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.713 |  | pass | 72 |
| a3pin6_own_G1_om0.5_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0014 (f[1,19,3]) | 1971 (sigma) | 2973 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.756 |  | pass | 366 |
| a3pin6_own_G1_om0.5_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0029 (sigma) | 1239 (sigma) | 1092 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.739 |  | pass | 68 |
| a3pin6_own_G1_om0.5_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0036 (lp__) | 974 (lp__) | 1459 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.692 |  | pass | 72 |
| a3pin6_own_G1_om0.5_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0049 (f[1,15,2]) | 1015 (omega) | 1791 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.743 |  | pass | 73 |
| a3pin6_own_G1_om0.5_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0021 (f[1,8,9]) | 1025 (lp__) | 1553 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.713 |  | pass | 75 |
| a3pin6_own_G1_om0_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0041 (f[1,5,2]) | 1250 (f[1,27,1]) | 2197 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.758 |  | pass | 69 |
| a3pin6_own_G1_om0_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0026 (f[1,11,1]) | 1130 (omega) | 1068 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.725 |  | pass | 73 |
| a3pin6_own_G1_om0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0040 (lp__) | 835 (f[1,21,1]) | 1555 (f[1,6,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.747 |  | pass | 63 |
| a3pin6_own_G1_om0_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0065 (f[1,26,1]) | 980 (sigma) | 1716 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.710 |  | pass | 72 |
| a3pin6_own_G1_om0_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0049 (omega) | 1262 (omega) | 2152 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.702 |  | pass | 185 |
| a3pin6_own_G1_om0_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0034 (omega) | 1124 (omega) | 984 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.667 |  | pass | 320 |
| a3pin6_own_G1_om1_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0063 (lp__) | 1216 (lp__) | 2266 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.778 |  | pass | 66 |
| a3pin6_own_G1_om1_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0039 (alpha[4]) | 906 (sigma) | 998 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.727 |  | pass | 74 |
| a3pin6_own_G1_om1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0053 (omega) | 1169 (lp__) | 1773 (f[1,18,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.749 |  | pass | 69 |
| a3pin6_own_G1_om1_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0043 (sigma) | 997 (sigma) | 853 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.744 |  | pass | 73 |
| a3pin6_own_G1_om1_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0036 (f[1,23,2]) | 913 (sigma) | 1298 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.756 |  | pass | 76 |
| a3pin6_own_G1_om1_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0043 (tau[1]) | 1026 (sigma) | 1495 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.699 |  | pass | 82 |
| a3pin6_own_L_om0.25_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0038 (sigma) | 1238 (alpha[1]) | 2269 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.713 |  | pass | 205 |
| a3pin6_own_L_om0.25_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0071 (omega) | 834 (omega) | 839 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.634 |  | pass | 393 |
| a3pin6_own_L_om0.25_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0030 (f[1,4,1]) | 1155 (omega) | 1709 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.742 |  | pass | 191 |
| a3pin6_own_L_om0.25_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0040 (lp__) | 700 (omega) | 672 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.668 |  | pass | 74 |
| a3pin6_own_L_om0.25_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0032 (f[1,27,2]) | 1682 (alpha[1]) | 2784 (f[1,18,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.709 |  | pass | 142 |
| a3pin6_own_L_om0.25_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0030 (alpha[5]) | 851 (omega) | 970 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.682 |  | pass | 68 |
| a3pin6_own_L_om0.5_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0040 (omega) | 925 (sigma) | 1537 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.725 |  | pass | 76 |
| a3pin6_own_L_om0.5_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0042 (f[1,25,2]) | 879 (sigma) | 2092 (f[1,19,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.749 |  | pass | 101 |
| a3pin6_own_L_om0.5_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0083 (lp__) | 1158 (lp__) | 1791 (f[1,11,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.688 |  | pass | 67 |
| a3pin6_own_L_om0.5_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0047 (omega) | 864 (sigma) | 1075 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.692 |  | pass | 84 |
| a3pin6_own_L_om0.5_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0060 (omega) | 1010 (alpha[1]) | 1940 (f[1,1,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.734 |  | pass | 58 |
| a3pin6_own_L_om0.5_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0095 (lp__) | 1185 (lp__) | 1217 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.709 |  | pass | 72 |
| a3pin6_own_L_om0_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0036 (lp__) | 1984 (omega) | 3202 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.756 |  | pass | 339 |
| a3pin6_own_L_om0_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0069 (sigma) | 907 (sigma) | 1542 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.666 |  | pass | 73 |
| a3pin6_own_L_om0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0053 (f[1,14,1]) | 548 (f[1,26,1]) | 1055 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.729 |  | pass | 53 |
| a3pin6_own_L_om0_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0069 (lp__) | 1022 (lp__) | 1128 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.699 |  | pass | 72 |
| a3pin6_own_L_om0_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0028 (f[1,7,2]) | 1081 (omega) | 1761 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.696 |  | pass | 65 |
| a3pin6_own_L_om0_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0041 (lp__) | 934 (sigma) | 1526 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.709 |  | pass | 73 |
| a3pin6_own_L_om1_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0037 (f[1,11,5]) | 987 (f[1,8,1]) | 1694 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.749 |  | pass | 96 |
| a3pin6_own_L_om1_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0034 (lp__) | 1086 (sigma) | 1534 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.794 |  | pass | 72 |
| a3pin6_own_L_om1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0047 (lp__) | 1167 (lp__) | 2122 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.739 |  | pass | 95 |
| a3pin6_own_L_om1_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0044 (sigma) | 674 (sigma) | 1165 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.758 |  | pass | 77 |
| a3pin6_own_L_om1_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0040 (f[1,18,2]) | 1217 (omega) | 1973 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.735 |  | pass | 229 |
| a3pin6_own_L_om1_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0086 (lp__) | 735 (sigma) | 1107 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.725 |  | pass | 76 |
| a3pin6_sim_S0_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0034 (f[1,8,2]) | 1583 (omega) | 1798 (f[1,28,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.736 |  | pass | 45 |
| a3pin6_sim_S0_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0028 (omega) | 1306 (omega) | 1263 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.713 |  | pass | 108 |
| a3pin6_sim_S0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0034 (f[1,1,1]) | 1171 (f[1,21,1]) | 2181 (f[1,12,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.727 |  | pass | 73 |
| a3pin6_sim_S0_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0021 (f[1,28,6]) | 1212 (lp__) | 2061 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.650 |  | pass | 75 |
| a3pin6_sim_S0_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0040 (f[1,4,2]) | 969 (sigma) | 1852 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.606 |  | pass | 57 |
| a3pin6_sim_S0_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0039 (f[1,10,2]) | 945 (lp__) | 1662 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.566 |  | pass | 106 |
| a3pin6_sim_S1_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0044 (lp__) | 809 (lp__) | 1454 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.570 |  | pass | 24 |
| a3pin6_sim_S1_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0060 (lp__) | 774 (lp__) | 1534 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.576 |  | pass | 26 |
| a3pin6_sim_S1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0040 (f[1,28,4]) | 1660 (lp__) | 1799 (f[1,23,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.805 |  | pass | 50 |
| a3pin6_sim_S1_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0048 (omega) | 585 (omega) | 1210 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.654 |  | pass | 70 |
| a3pin6_sim_S1_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0043 (sigma) | 986 (lp__) | 1501 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.621 |  | pass | 46 |
| a3pin6_sim_S1_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0013 (alpha[6]) | 1341 (omega) | 3226 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.630 |  | pass | 251 |
| a3pin6_sim_S2_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0044 (f[1,8,2]) | 1357 (sigma) | 795 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.842 |  | pass | 47 |
| a3pin6_sim_S2_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0021 (omega) | 2403 (omega) | 2515 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.796 |  | pass | 334 |
| a3pin6_sim_S2_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0030 (f[1,18,7]) | 962 (lp__) | 912 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.650 |  | pass | 41 |
| a3pin6_sim_S2_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0079 (omega) | 410 (omega) | 541 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.560 |  | pass | 72 |
| a3pin6_sim_S2_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0032 (f[1,2,4]) | 1056 (sigma) | 753 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.773 |  | pass | 49 |
| a3pin6_sim_S2_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0052 (f[1,22,8]) | 507 (omega) | 938 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.695 |  | pass | 69 |
| a3pin6_slope_rho30_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0063 (f[1,5,4]) | 844 (f[1,27,1]) | 1454 (f[1,10,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.769 |  | pass | 104 |
| a3pin6_slope_rho30_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0057 (sigma) | 757 (sigma) | 1111 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.797 |  | pass | 73 |
| a3pin6_slope_rho30_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0035 (lp__) | 803 (omega) | 792 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.767 |  | pass | 81 |
| a3pin6_slope_rho30_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0067 (lp__) | 1022 (omega) | 1095 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.759 |  | pass | 69 |
| a3pin6_slope_rho30_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0043 (f[1,8,2]) | 1084 (f[1,22,1]) | 2339 (f[1,22,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.801 |  | pass | 264 |
| a3pin6_slope_rho30_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0077 (f[1,19,2]) | 787 (omega) | 1024 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.722 |  | pass | 71 |
| a3pin6_slope_rho60_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0070 (f[1,9,2]) | 825 (alpha[1]) | 1080 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.712 |  | pass | 68 |
| a3pin6_slope_rho60_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0038 (omega) | 932 (lp__) | 1606 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.685 |  | pass | 71 |
| a3pin6_slope_rho60_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0079 (lp__) | 690 (f[1,28,1]) | 1180 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.748 |  | pass | 77 |
| a3pin6_slope_rho60_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0048 (f[1,3,7]) | 890 (lp__) | 1522 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.731 |  | pass | 72 |
| a3pin6_slope_rho60_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0049 (lp__) | 934 (omega) | 1439 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.728 |  | pass | 72 |
| a3pin6_slope_rho60_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0052 (lp__) | 1239 (lp__) | 1911 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.730 |  | pass | 70 |

## A3 with the range pinned at 12 km (a3_pinned_range.R 12; a3_pinned_rho12.md and _alphaN0-2.md)

78 records, 76 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3pin12_own_G1_om0.25_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0018 (sigma) | 2382 (omega) | 3132 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.733 |  | pass | 223 |
| a3pin12_own_G1_om0.25_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0036 (sigma) | 1402 (sigma) | 1800 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.746 |  | pass | 271 |
| a3pin12_own_G1_om0.25_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0029 (f[1,25,2]) | 1952 (omega) | 2390 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.727 |  | pass | 276 |
| a3pin12_own_G1_om0.25_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0099 (omega) | 805 (omega) | 1028 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.690 |  | pass | 75 |
| a3pin12_own_G1_om0.25_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0047 (f[1,19,2]) | 1163 (lp__) | 1982 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.741 |  | pass | 68 |
| a3pin12_own_G1_om0.25_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0056 (lp__) | 867 (lp__) | 1284 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.751 |  | pass | 74 |
| a3pin12_own_G1_om0.5_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0067 (tau[1]) | 782 (tau[1]) | 1196 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.729 |  | pass | 71 |
| a3pin12_own_G1_om0.5_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0030 (log_lik[11]) | 867 (tau[1]) | 1432 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.733 |  | pass | 77 |
| a3pin12_own_G1_om0.5_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0051 (f[1,24,1]) | 1281 (lp__) | 2378 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.718 |  | pass | 69 |
| a3pin12_own_G1_om0.5_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0015 (lp__) | 2464 (sigma) | 3843 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.730 |  | pass | 264 |
| a3pin12_own_G1_om0.5_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0062 (f[1,27,2]) | 1056 (f[1,22,1]) | 2019 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.732 |  | pass | 70 |
| a3pin12_own_G1_om0.5_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0063 (tau[1]) | 812 (tau[1]) | 1281 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.710 |  | pass | 76 |
| a3pin12_own_G1_om0_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0089 (lp__) | 1309 (omega) | 1426 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.745 |  | pass | 73 |
| a3pin12_own_G1_om0_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0035 (sigma) | 702 (sigma) | 817 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.685 |  | pass | 74 |
| a3pin12_own_G1_om0_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0030 (lp__) | 1896 (f[1,22,1]) | 3143 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.757 |  | pass | 211 |
| a3pin12_own_G1_om0_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0049 (omega) | 981 (omega) | 892 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.725 |  | pass | 72 |
| a3pin12_own_G1_om0_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0068 (f[1,27,1]) | 1353 (f[1,27,1]) | 2606 (f[1,8,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.732 |  | fail: divergences | 204 |
| a3pin12_own_G1_om0_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0074 (omega) | 1109 (omega) | 893 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.666 |  | pass | 292 |
| a3pin12_own_G1_om1_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0033 (f[1,19,1]) | 1792 (f[1,22,1]) | 2664 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.777 |  | pass | 244 |
| a3pin12_own_G1_om1_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0061 (tau[1]) | 1126 (sigma) | 1502 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.702 |  | pass | 102 |
| a3pin12_own_G1_om1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0027 (f[1,18,5]) | 1202 (omega) | 1768 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.733 |  | pass | 70 |
| a3pin12_own_G1_om1_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0025 (f[1,12,4]) | 969 (lp__) | 1621 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.729 |  | pass | 75 |
| a3pin12_own_G1_om1_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0062 (f[1,4,1]) | 753 (tau[1]) | 1013 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.778 |  | pass | 95 |
| a3pin12_own_G1_om1_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0063 (tau[1]) | 709 (tau[1]) | 1246 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.755 |  | pass | 95 |
| a3pin12_own_L_om0.25_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0057 (f[1,27,1]) | 1110 (f[1,21,1]) | 2139 (f[1,1,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.743 |  | pass | 178 |
| a3pin12_own_L_om0.25_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0091 (omega) | 773 (omega) | 739 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.660 |  | pass | 321 |
| a3pin12_own_L_om0.25_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0030 (f[1,15,2]) | 1672 (alpha[1]) | 3110 (f[1,22,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.725 |  | pass | 211 |
| a3pin12_own_L_om0.25_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0059 (omega) | 676 (omega) | 879 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.665 |  | pass | 267 |
| a3pin12_own_L_om0.25_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0057 (f[1,18,1]) | 1239 (f[1,20,1]) | 2527 (f[1,6,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.707 |  | fail: divergences | 204 |
| a3pin12_own_L_om0.25_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0027 (f[1,21,6]) | 1026 (omega) | 1211 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.701 |  | pass | 70 |
| a3pin12_own_L_om0.5_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0037 (f[1,12,1]) | 992 (lp__) | 1834 (f[1,2,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.722 |  | pass | 77 |
| a3pin12_own_L_om0.5_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0034 (f[1,4,2]) | 849 (sigma) | 1698 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.758 |  | pass | 75 |
| a3pin12_own_L_om0.5_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0065 (omega) | 961 (omega) | 1140 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.728 |  | pass | 68 |
| a3pin12_own_L_om0.5_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0084 (f[1,11,9]) | 827 (sigma) | 1654 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.744 |  | pass | 77 |
| a3pin12_own_L_om0.5_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0033 (f[1,6,1]) | 1926 (f[1,22,1]) | 2601 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.758 |  | pass | 229 |
| a3pin12_own_L_om0.5_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0065 (omega) | 746 (omega) | 1324 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.705 |  | pass | 71 |
| a3pin12_own_L_om0_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0040 (omega) | 1508 (omega) | 2413 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.746 |  | pass | 235 |
| a3pin12_own_L_om0_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0084 (omega) | 979 (omega) | 719 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.758 |  | pass | 72 |
| a3pin12_own_L_om0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0049 (omega) | 636 (f[1,22,1]) | 989 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.763 |  | pass | 60 |
| a3pin12_own_L_om0_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0054 (omega) | 1122 (omega) | 1014 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.744 |  | pass | 71 |
| a3pin12_own_L_om0_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0053 (f[1,24,2]) | 761 (omega) | 1089 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.736 |  | pass | 66 |
| a3pin12_own_L_om0_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0042 (sigma) | 842 (sigma) | 1872 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.671 |  | pass | 72 |
| a3pin12_own_L_om1_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0040 (f[1,18,1]) | 807 (f[1,8,1]) | 1293 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.712 |  | pass | 77 |
| a3pin12_own_L_om1_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0035 (sigma) | 1199 (sigma) | 1898 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.748 |  | pass | 72 |
| a3pin12_own_L_om1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0060 (f[1,6,1]) | 873 (omega) | 941 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.748 |  | pass | 70 |
| a3pin12_own_L_om1_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0080 (omega) | 1024 (sigma) | 851 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.731 |  | pass | 95 |
| a3pin12_own_L_om1_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0057 (f[1,20,2]) | 1473 (alpha[1]) | 2751 (f[1,2,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.746 |  | pass | 243 |
| a3pin12_own_L_om1_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0044 (f[1,7,2]) | 713 (sigma) | 1003 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.719 |  | pass | 80 |
| a3pin12_sim_S0_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0038 (f[1,8,5]) | 1489 (omega) | 2366 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.771 |  | pass | 57 |
| a3pin12_sim_S0_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0032 (omega) | 1983 (omega) | 1604 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.682 |  | pass | 290 |
| a3pin12_sim_S0_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0043 (omega) | 1248 (f[1,22,1]) | 2407 (f[1,25,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.695 |  | pass | 73 |
| a3pin12_sim_S0_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0031 (lp__) | 1270 (lp__) | 2084 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.639 |  | pass | 96 |
| a3pin12_sim_S0_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0030 (omega) | 1277 (lp__) | 2125 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.638 |  | pass | 58 |
| a3pin12_sim_S0_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0037 (lp__) | 1191 (lp__) | 1296 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.656 |  | pass | 108 |
| a3pin12_sim_S1_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0056 (sigma) | 671 (sigma) | 1257 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.594 |  | pass | 32 |
| a3pin12_sim_S1_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0035 (sigma) | 598 (sigma) | 1284 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.596 |  | pass | 28 |
| a3pin12_sim_S1_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0028 (omega) | 1365 (omega) | 2039 (f[1,6,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.770 |  | pass | 43 |
| a3pin12_sim_S1_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0025 (lp__) | 1271 (omega) | 2811 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.680 |  | pass | 268 |
| a3pin12_sim_S1_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0041 (lp__) | 831 (sigma) | 944 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.632 |  | pass | 39 |
| a3pin12_sim_S1_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0051 (sigma) | 1081 (omega) | 2584 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.590 |  | pass | 244 |
| a3pin12_sim_S2_r1 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0039 (f[1,24,6]) | 712 (sigma) | 416 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.813 |  | pass | 56 |
| a3pin12_sim_S2_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0028 (sigma) | 1024 (omega) | 1165 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.728 |  | pass | 72 |
| a3pin12_sim_S2_r2 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0091 (lp__) | 910 (sigma) | 842 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.637 |  | pass | 51 |
| a3pin12_sim_S2_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0063 (omega) | 778 (omega) | 1716 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.599 |  | pass | 254 |
| a3pin12_sim_S2_r3 | Normal(0, 10) | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0050 (lp__) | 1225 (sigma) | 1131 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.693 |  | pass | 39 |
| a3pin12_sim_S2_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0056 (omega) | 626 (omega) | 1558 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.683 |  | pass | 67 |
| a3pin12_slope_rho30_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0014 (sigma) | 2154 (alpha[1]) | 3014 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.763 |  | pass | 430 |
| a3pin12_slope_rho30_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0048 (f[1,15,1]) | 739 (sigma) | 1097 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.744 |  | pass | 70 |
| a3pin12_slope_rho30_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0024 (f[1,2,1]) | 2130 (omega) | 4280 (f[1,20,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.761 |  | pass | 384 |
| a3pin12_slope_rho30_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0071 (omega) | 827 (omega) | 898 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.723 |  | pass | 67 |
| a3pin12_slope_rho30_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0030 (f[1,18,1]) | 1562 (f[1,26,1]) | 2733 (f[1,3,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.761 |  | pass | 290 |
| a3pin12_slope_rho30_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0082 (sigma) | 591 (omega) | 848 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.743 |  | pass | 74 |
| a3pin12_slope_rho60_r1 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0029 (f[1,16,4]) | 1404 (f[1,27,1]) | 2292 (f[1,1,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.719 |  | pass | 256 |
| a3pin12_slope_rho60_r1 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0050 (f[1,6,1]) | 970 (omega) | 963 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.669 |  | pass | 69 |
| a3pin12_slope_rho60_r2 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0018 (f[1,11,9]) | 1531 (omega) | 2617 (f[1,17,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.752 |  | pass | 277 |
| a3pin12_slope_rho60_r2 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0033 (lp__) | 894 (omega) | 825 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.742 |  | pass | 70 |
| a3pin12_slope_rho60_r3 | Normal(0, 10) | 2026-10-18 | 2 (6000/3000, 0.99, 13) | 1.0026 (f[1,8,2]) | 2102 (omega) | 3007 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.760 |  | pass | 388 |
| a3pin12_slope_rho60_r3 | Normal(0, 2), hard-coded | 2026-10-18 | 1 (3000/1500, 0.95, 12) | 1.0039 (lp__) | 982 (sigma) | 1389 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.734 |  | pass | 72 |

## Other probe fits (mechanism, mode and nugget checks; exploratory, not reported)

16 records, 8 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mech_a1_rho12_sg2_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | nan (comembership[22,20]) | 4 (comembership[22,20]) | 11 (coherence) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.740 | 1.000 | fail: rhat, ess_bulk, ess_tail, stuck, comembership_split | 704 |
| mech_a1_scattered_sg0.5 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0621 (comembership[22,21]) | 47 (comembership[22,21]) | 192 (comembership[24,4]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.627 | 0.254 | fail: rhat, ess_bulk, ess_tail | 660 |
| mech_a1_scattered_sg1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0828 (comembership[25,20]) | 39 (comembership[25,20]) | 89 (comembership[28,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.639 | 0.254 | fail: rhat, ess_bulk, ess_tail | 447 |
| modecheck_a0_scattered_sg1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0038 (sigma) | 965 (tau[1]) | 1639 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.737 |  | pass | 76 |
| modecheck_a0_scattered_sg1_nonug | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0065 (rho) | 1001 (rho) | 1282 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.730 |  | pass | 77 |
| modecheck_a0_sg2_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0074 (rho) | 872 (tau[1]) | 1740 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.765 |  | pass | 65 |
| nonug_a0_scattered_sg0.5 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0036 (sigma) | 996 (tau[1]) | 1602 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.722 |  | pass | 78 |
| nonug_a0_scattered_sg1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0038 (sigma) | 965 (tau[1]) | 1639 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.737 |  | pass | 76 |
| nonug_a0_sg0.5_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0046 (tau[1]) | 992 (tau[1]) | 1806 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.764 |  | pass | 74 |
| nonug_a0_sg2_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0074 (rho) | 872 (tau[1]) | 1740 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.765 |  | pass | 64 |
| nonug_a0_zero | Normal(0, 2), hard-coded | 2026-10-17 | 1 (3000/1500, 0.95, 12) | 1.0041 (f[5,3]) | 960 (tau[1]) | 1729 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.746 |  | pass | 59 |
| nonug_a1_scattered_sg0.5 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.1699 (comembership[21,19]) | 16 (comembership[22,11]) | 58 (comembership[27,24]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.727 | 0.552 | fail: rhat, ess_bulk, ess_tail, comembership_split | 645 |
| nonug_a1_scattered_sg1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.1574 (comembership[22,21]) | 20 (comembership[22,21]) | 12 (comembership[22,3]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.738 | 0.415 | fail: rhat, ess_bulk, ess_tail | 619 |
| nonug_a1_sg0.5_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0128 (comembership[25,21]) | 537 (comembership[22,21]) | 1029 (comembership[27,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.691 | 0.115 | fail: rhat | 614 |
| nonug_a1_sg2_r1 | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0558 (comembership[25,3]) | 58 (comembership[25,3]) | 150 (comembership[11,7]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.794 | 0.219 | fail: rhat, ess_bulk, ess_tail | 621 |
| nonug_a1_zero | Normal(0, 2), hard-coded | 2026-10-17 | 2 (6000/3000, 0.99, 13) | 1.0941 (comembership[25,2]) | 34 (comembership[25,2]) | 126 (comembership[25,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.720 | 0.379 | fail: rhat, ess_bulk, ess_tail | 618 |
