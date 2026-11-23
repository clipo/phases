# Per-fit sampler diagnostics for the latent-group probe's Stan fits (rule 16)

Produced by `analyses/probe_latent_groups/export_fit_diagnostics.py` from the fit records in `output/gp_cache/` (gitignored; 147 records with diagnostics). Every value is the one the gate read at the final rung; the gate thresholds are R-hat at most 1.01 on gated variables, bulk and tail ESS at least 400, no divergent transitions, E-BFMI at least 0.3 and, for A1, a per-chain co-membership split of at most 0.5. Treedepth saturation is reported, not gated. Rung settings are iterations per chain / warmup, adapt_delta and max_treedepth; four chains throughout. A label with two rows was fitted under two versions of the Stan source, with the class-intercept prior hard-coded as Normal(0, 2) and, after D25 (2026-10-18), read from the data (default Normal(0, 10)); the findings tables named per family report the data-driven version under the plain name and the hard-coded version as `*_alphaN0-2.md`. Descriptive export; exploratory design validation, not evidence about the record. Columns in `fit_diagnostics.csv`.

## A0 spatial field: basin fit, prior predictive, sensitivity, recovery, intercept priors (refit_a0_basin.R, a0_intercept_prior.R)

11 records, 11 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a0_alpha_P10 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0025 (lp__) | 898 (tau[1]) | 1778 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.673 |  | pass | 110 |
| a0_alpha_P2 | Normal(0, 2) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0027 (tau[1]) | 1260 (tau[1]) | 2241 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.711 |  | pass | 359 |
| a0_alpha_PC | Normal(pooled log-odds, 2) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0040 (lp__) | 878 (lp__) | 1894 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.696 |  | pass | 60 |
| a0_basin | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0025 (lp__) | 898 (tau[1]) | 1778 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.673 |  | pass | 116 |
| a0_prior_predictive | Normal(0, 10) | 2026-11-22 | 1 (2000/1000, 0.9, 10) | 1.0048 (f[4,1]) | 1315 (lp__) | 1525 (tau[1]) | 0/4000 (0.00%) | 0/4000 (0.0%; deepest 4) | 0.965 |  | pass | 3 |
| a0_recovery_share0.2_rho5 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0044 (alpha[2]) | 955 (alpha[2]) | 1236 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.706 |  | pass | 186 |
| a0_recovery_share0.2_rho60 | Normal(0, 10) | 2026-11-22 | 3 (12000/6000, 0.999, 14) | 1.0023 (alpha[9]) | 2054 (rho) | 2411 (alpha[6]) | 0/24000 (0.00%) | 0/24000 (0.0%; deepest 13) | 0.674 |  | pass | 733 |
| a0_recovery_share0.9_rho5 | Normal(0, 10) | 2026-11-22 | 4 (16000/8000, 0.9999, 15) | 1.0023 (rho) | 2416 (rho) | 1444 (rho) | 0/32000 (0.00%) | 0/32000 (0.0%; deepest 12) | 0.672 |  | pass | 486 |
| a0_recovery_share0.9_rho60 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0061 (tau[1]) | 1157 (rho) | 2405 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.727 |  | pass | 59 |
| a0_sens_rho_uniform | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0068 (tau[1]) | 930 (tau[1]) | 1817 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.654 |  | pass | 127 |
| a0_sens_scale05 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0063 (tau[1]) | 961 (tau[1]) | 1813 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.676 |  | pass | 56 |

## A1 latent-unit recovery and equivalence fits (validate_a1.R; a1_validation.md)

17 records, 4 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a1val_a0_on_zero | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0029 (tau[1]) | 1150 (tau[1]) | 1324 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.769 |  | pass | 132 |
| a1val_a1_fixed0_on_zero | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0030 (sigma) | 1446 (tau[1]) | 2020 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.735 |  | pass | 257 |
| a1val_a1_on_zero | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0179 (comembership[5,2]) | 386 (comembership[27,19]) | 821 (comembership[7,5]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.674 | 0.096 | fail: rhat, ess_bulk | 975 |
| a1val_sg0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0060 (Delta_w) | 533 (comembership[3,2]) | 1411 (comembership[14,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.721 | 0.070 | pass | 917 |
| a1val_sg0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0415 (comembership[20,2]) | 123 (comembership[20,2]) | 292 (Delta_w) | 4/12000 (0.03%) | 0/12000 (0.0%; deepest 10) | 0.667 | 0.119 | fail: rhat, ess_bulk, ess_tail, divergences | 444 |
| a1val_sg0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0146 (comembership[26,16]) | 419 (comembership[3,2]) | 763 (comembership[27,9]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.693 | 0.105 | fail: rhat | 966 |
| a1val_sg0.25_r4 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0200 (comembership[22,3]) | 380 (comembership[22,3]) | 572 (comembership[25,5]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.634 | 0.117 | fail: rhat, ess_bulk | 891 |
| a1val_sg0.25_r5 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0116 (comembership[22,11]) | 472 (comembership[25,19]) | 1114 (comembership[13,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.664 | 0.116 | fail: rhat | 877 |
| a1val_sg0.5_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0106 (comembership[19,8]) | 873 (comembership[19,3]) | 1108 (comembership[8,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.715 | 0.090 | fail: rhat | 976 |
| a1val_sg0.5_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0154 (comembership[20,2]) | 403 (comembership[13,5]) | 372 (comembership[13,7]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.670 | 0.101 | fail: rhat, ess_tail | 847 |
| a1val_sg0.5_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0389 (comembership[19,5]) | 150 (comembership[5,3]) | 291 (comembership[28,3]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.722 | 0.184 | fail: rhat, ess_bulk, ess_tail | 699 |
| a1val_sg0.5_r4 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0088 (comembership[22,20]) | 632 (comembership[21,3]) | 1045 (comembership[21,18]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.656 | 0.071 | pass | 1183 |
| a1val_sg0.5_r5 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0110 (comembership[25,22]) | 593 (comembership[25,21]) | 916 (comembership[20,16]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.677 | 0.117 | fail: rhat | 733 |
| a1val_sg1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0900 (comembership[22,16]) | 38 (comembership[22,19]) | 47 (comembership[17,2]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 11) | 0.626 | 0.374 | fail: rhat, ess_bulk, ess_tail, divergences | 610 |
| a1val_sg1_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0127 (comembership[25,21]) | 484 (comembership[21,19]) | 992 (comembership[9,3]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.656 | 0.137 | fail: rhat | 415 |
| a1val_sg2_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.6111 (comembership[27,20]) | 7 (comembership[27,20]) | 6 (comembership[20,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.745 | 0.726 | fail: rhat, ess_bulk, ess_tail, comembership_split | 1014 |
| a1val_sg2_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0206 (comembership[21,8]) | 205 (comembership[25,22]) | 548 (comembership[23,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.710 | 0.164 | fail: rhat, ess_bulk | 690 |

## A3 barrier recovery and equivalence fits, range free (validate_a3.R; a3_validation.md and a3_validation_alphaN0-2.md)

41 records, 41 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3val_a0_on_om0 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0022 (f[4,6]) | 1129 (lp__) | 2042 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.715 |  | pass | 134 |
| a3val_a3_fixed0_on_om0 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0029 (f[1,15,1]) | 1182 (tau[1]) | 1495 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.760 |  | pass | 184 |
| a3val_own_G1_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0024 (rho) | 1270 (rho) | 1691 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.708 |  | pass | 554 |
| a3val_own_G1_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (sigma) | 1925 (rho) | 3460 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.713 |  | pass | 923 |
| a3val_own_G1_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0096 (sigma) | 460 (rho) | 825 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.735 |  | pass | 122 |
| a3val_own_G1_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0073 (tau[1]) | 636 (rho) | 544 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.729 |  | pass | 192 |
| a3val_own_G1_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0073 (rho) | 877 (tau[1]) | 1293 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.714 |  | pass | 151 |
| a3val_own_G1_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0066 (lp__) | 968 (tau[1]) | 1448 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 190 |
| a3val_own_G1_om0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0035 (sigma) | 1291 (lp__) | 2090 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.738 |  | pass | 179 |
| a3val_own_G1_om0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0055 (sigma) | 1210 (tau[1]) | 1494 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.764 |  | pass | 190 |
| a3val_own_G1_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0039 (rho) | 1331 (rho) | 1880 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.679 |  | pass | 551 |
| a3val_own_G1_om1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0037 (rho) | 1694 (rho) | 2862 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.793 |  | pass | 384 |
| a3val_own_G1_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0032 (omega) | 1034 (rho) | 1773 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.735 |  | pass | 166 |
| a3val_own_G1_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0074 (lp__) | 1150 (tau[1]) | 1304 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.799 |  | pass | 167 |
| a3val_own_L_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0084 (f[1,17,1]) | 577 (tau[1]) | 1168 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 134 |
| a3val_own_L_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0048 (rho) | 1329 (rho) | 1558 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.755 |  | pass | 151 |
| a3val_own_L_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0025 (sigma) | 1617 (rho) | 1934 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.748 |  | pass | 356 |
| a3val_own_L_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0049 (f[1,18,1]) | 1141 (rho) | 1400 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.780 |  | pass | 161 |
| a3val_own_L_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0052 (f[1,15,1]) | 1029 (rho) | 1719 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.726 |  | pass | 152 |
| a3val_own_L_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0054 (tau[1]) | 1261 (lp__) | 2090 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.717 |  | pass | 148 |
| a3val_own_L_om0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0071 (lp__) | 1346 (tau[1]) | 1680 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.734 |  | pass | 155 |
| a3val_own_L_om0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0034 (lp__) | 1357 (lp__) | 2304 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.726 |  | pass | 175 |
| a3val_own_L_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0020 (f[1,5,2]) | 1873 (rho) | 3117 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.757 |  | pass | 503 |
| a3val_own_L_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0047 (sigma) | 1042 (rho) | 1001 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.759 |  | pass | 150 |
| a3val_own_L_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0035 (f[1,13,3]) | 1566 (omega) | 2771 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.778 |  | pass | 155 |
| a3val_own_L_om1_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0030 (f[1,14,2]) | 1186 (rho) | 1825 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.788 |  | pass | 413 |
| a3val_sim_S0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0011 (lp__) | 2665 (omega) | 3133 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.753 |  | pass | 268 |
| a3val_sim_S0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0029 (lp__) | 1526 (lp__) | 2325 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.737 |  | pass | 165 |
| a3val_sim_S0_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0026 (f[1,16,2]) | 1404 (lp__) | 1767 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.699 |  | pass | 94 |
| a3val_sim_S1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0015 (omega) | 1244 (sigma) | 1541 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.571 |  | pass | 274 |
| a3val_sim_S1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0025 (f[1,18,2]) | 1414 (sigma) | 1843 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.798 |  | pass | 83 |
| a3val_sim_S1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0085 (lp__) | 615 (sigma) | 1115 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.532 |  | pass | 89 |
| a3val_sim_S2_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0047 (sigma) | 1266 (sigma) | 2639 (f[1,10,5]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.784 |  | pass | 83 |
| a3val_sim_S2_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0035 (sigma) | 1280 (sigma) | 2085 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.662 |  | pass | 256 |
| a3val_sim_S2_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0046 (tau[1]) | 791 (sigma) | 1153 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.757 |  | pass | 47 |
| a3val_slope_rho30_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0019 (rho) | 1411 (rho) | 2753 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.772 |  | pass | 485 |
| a3val_slope_rho30_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0026 (f[1,10,4]) | 1516 (sigma) | 2432 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.745 |  | pass | 432 |
| a3val_slope_rho30_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0024 (f[1,24,2]) | 1609 (rho) | 3241 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.742 |  | pass | 553 |
| a3val_slope_rho60_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0044 (f[1,13,4]) | 720 (tau[1]) | 1232 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.727 |  | pass | 118 |
| a3val_slope_rho60_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0023 (f[1,4,7]) | 1415 (lp__) | 2228 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.752 |  | pass | 142 |
| a3val_slope_rho60_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0054 (tau[1]) | 1237 (tau[1]) | 2235 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.720 |  | pass | 555 |

## A3 with the range pinned at 6 km (a3_pinned_range.R 6; a3_pinned_rho6.md and _alphaN0-2.md)

39 records, 39 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3pin6_own_G1_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0092 (omega) | 837 (omega) | 1261 (f[1,20,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.692 |  | pass | 80 |
| a3pin6_own_G1_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0038 (lp__) | 1250 (lp__) | 1249 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.689 |  | pass | 86 |
| a3pin6_own_G1_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0050 (f[1,11,1]) | 822 (f[1,20,1]) | 1562 (f[1,17,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.760 |  | pass | 75 |
| a3pin6_own_G1_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0053 (f[1,8,1]) | 1166 (sigma) | 1260 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.742 |  | pass | 95 |
| a3pin6_own_G1_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0067 (omega) | 1140 (omega) | 1826 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.721 |  | pass | 84 |
| a3pin6_own_G1_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0075 (lp__) | 955 (lp__) | 1439 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.697 |  | pass | 101 |
| a3pin6_own_G1_om0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0065 (omega) | 556 (f[1,22,1]) | 1211 (f[1,20,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.755 |  | pass | 82 |
| a3pin6_own_G1_om0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0078 (f[1,12,1]) | 603 (f[1,3,1]) | 952 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.730 |  | pass | 73 |
| a3pin6_own_G1_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0056 (omega) | 1094 (omega) | 1719 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.728 |  | pass | 202 |
| a3pin6_own_G1_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0039 (f[1,26,1]) | 1204 (f[1,20,1]) | 1884 (f[1,27,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.776 |  | pass | 86 |
| a3pin6_own_G1_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0054 (lp__) | 1182 (lp__) | 2352 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.736 |  | pass | 96 |
| a3pin6_own_G1_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0042 (sigma) | 858 (omega) | 1120 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.723 |  | pass | 102 |
| a3pin6_own_L_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (f[1,20,1]) | 1076 (alpha[1]) | 1978 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.702 |  | pass | 203 |
| a3pin6_own_L_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0037 (f[1,18,2]) | 1499 (omega) | 2324 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.752 |  | pass | 293 |
| a3pin6_own_L_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0040 (f[1,11,1]) | 1517 (f[1,3,1]) | 3007 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.737 |  | pass | 295 |
| a3pin6_own_L_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0066 (f[1,21,2]) | 726 (f[1,27,1]) | 1566 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.761 |  | pass | 104 |
| a3pin6_own_L_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0079 (f[1,8,1]) | 649 (alpha[1]) | 999 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.668 |  | pass | 74 |
| a3pin6_own_L_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0081 (f[1,17,4]) | 894 (f[1,27,2]) | 1233 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.717 |  | pass | 85 |
| a3pin6_own_L_om0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0031 (f[1,23,1]) | 1488 (alpha[1]) | 2217 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.747 |  | pass | 268 |
| a3pin6_own_L_om0_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0031 (f[1,26,2]) | 1225 (omega) | 1609 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 9) | 0.691 |  | pass | 212 |
| a3pin6_own_L_om0_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0051 (omega) | 635 (omega) | 1019 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.702 |  | pass | 52 |
| a3pin6_own_L_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0086 (f[1,18,1]) | 745 (omega) | 1354 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.725 |  | pass | 116 |
| a3pin6_own_L_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0040 (f[1,25,1]) | 1150 (omega) | 1934 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.752 |  | pass | 102 |
| a3pin6_own_L_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0094 (sigma) | 596 (omega) | 845 (f[1,12,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.715 |  | pass | 99 |
| a3pin6_sim_S0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0032 (f[1,19,7]) | 1309 (lp__) | 2239 (f[1,17,5]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.744 |  | pass | 60 |
| a3pin6_sim_S0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0037 (f[1,12,2]) | 1266 (f[1,22,1]) | 2100 (f[1,8,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.679 |  | pass | 100 |
| a3pin6_sim_S0_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0049 (lp__) | 999 (lp__) | 1898 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.609 |  | pass | 91 |
| a3pin6_sim_S1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0057 (omega) | 742 (lp__) | 1524 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.537 |  | pass | 40 |
| a3pin6_sim_S1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0053 (f[1,23,1]) | 1538 (alpha[1]) | 1987 (f[1,8,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.771 |  | pass | 66 |
| a3pin6_sim_S1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0072 (lp__) | 1022 (lp__) | 1856 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.634 |  | pass | 54 |
| a3pin6_sim_S2_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0023 (f[1,23,5]) | 1782 (f[1,20,2]) | 1561 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.770 |  | pass | 73 |
| a3pin6_sim_S2_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0044 (sigma) | 1088 (sigma) | 1638 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.664 |  | pass | 53 |
| a3pin6_sim_S2_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0035 (f[1,26,4]) | 1115 (sigma) | 584 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.711 |  | pass | 56 |
| a3pin6_slope_rho30_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0044 (f[1,8,4]) | 905 (alpha[1]) | 1386 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.756 |  | pass | 101 |
| a3pin6_slope_rho30_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0050 (f[1,21,4]) | 928 (f[1,22,1]) | 1745 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 87 |
| a3pin6_slope_rho30_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0020 (lp__) | 1080 (f[1,22,1]) | 2300 (f[1,13,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.753 |  | pass | 309 |
| a3pin6_slope_rho60_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0036 (f[1,9,1]) | 1308 (lp__) | 2301 (f[1,12,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.713 |  | pass | 89 |
| a3pin6_slope_rho60_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0042 (f[1,8,1]) | 981 (omega) | 1171 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.767 |  | pass | 86 |
| a3pin6_slope_rho60_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0096 (lp__) | 1205 (lp__) | 1711 (f[1,6,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.657 |  | pass | 94 |

## A3 with the range pinned at 12 km (a3_pinned_range.R 12; a3_pinned_rho12.md and _alphaN0-2.md)

39 records, 31 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3pin12_own_G1_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0026 (f[1,15,2]) | 1446 (omega) | 2143 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.743 |  | pass | 273 |
| a3pin12_own_G1_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0018 (lp__) | 2199 (omega) | 4042 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.702 |  | pass | 346 |
| a3pin12_own_G1_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0031 (omega) | 1407 (lp__) | 1885 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.758 |  | pass | 81 |
| a3pin12_own_G1_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0027 (f[1,10,2]) | 520 (tau[1]) | 797 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.739 |  | pass | 93 |
| a3pin12_own_G1_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0059 (omega) | 958 (omega) | 1445 (f[1,19,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.719 |  | pass | 80 |
| a3pin12_own_G1_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0046 (tau[1]) | 964 (tau[1]) | 1366 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.738 |  | pass | 100 |
| a3pin12_own_G1_om0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0018 (lp__) | 1916 (f[1,26,1]) | 3007 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.743 |  | pass | 339 |
| a3pin12_own_G1_om0_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0061 (f[1,23,1]) | 1434 (f[1,20,1]) | 2798 (f[1,21,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.730 |  | pass | 258 |
| a3pin12_own_G1_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (tau[1]) | 1591 (omega) | 3270 (f[1,27,1]) | 2/12000 (0.02%) | 0/12000 (0.0%; deepest 10) | 0.706 |  | fail: divergences | 224 |
| a3pin12_own_G1_om1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0025 (f[1,1,1]) | 1742 (f[1,20,1]) | 2313 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.704 |  | pass | 321 |
| a3pin12_own_G1_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0046 (lp__) | 1308 (lp__) | 1817 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.753 |  | pass | 104 |
| a3pin12_own_G1_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0089 (tau[1]) | 891 (tau[1]) | 1788 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.742 |  | pass | 161 |
| a3pin12_own_L_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0016 (f[1,11,1]) | 1078 (alpha[1]) | 2159 (f[1,19,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 9) | 0.732 |  | fail: divergences | 205 |
| a3pin12_own_L_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0015 (tau[1]) | 2072 (alpha[1]) | 4038 (f[1,20,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.721 |  | fail: divergences | 270 |
| a3pin12_own_L_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0023 (f[1,3,1]) | 1902 (alpha[1]) | 2664 (f[1,20,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.769 |  | pass | 244 |
| a3pin12_own_L_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0034 (omega) | 929 (f[1,20,1]) | 2269 (f[1,10,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.810 |  | pass | 98 |
| a3pin12_own_L_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0026 (omega) | 1772 (omega) | 3084 (f[1,25,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.729 |  | pass | 281 |
| a3pin12_own_L_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0054 (omega) | 945 (f[1,22,1]) | 871 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.743 |  | pass | 75 |
| a3pin12_own_L_om0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0031 (f[1,3,2]) | 1839 (f[1,20,1]) | 3691 (alpha[2]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.722 |  | fail: divergences | 294 |
| a3pin12_own_L_om0_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0025 (omega) | 1833 (omega) | 3941 (f[1,21,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 11) | 0.729 |  | fail: divergences | 346 |
| a3pin12_own_L_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0036 (omega) | 1778 (omega) | 3691 (f[1,22,2]) | 2/12000 (0.02%) | 0/12000 (0.0%; deepest 10) | 0.711 |  | fail: divergences | 272 |
| a3pin12_own_L_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0068 (f[1,1,7]) | 711 (omega) | 1184 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.755 |  | pass | 99 |
| a3pin12_own_L_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0030 (f[1,25,1]) | 974 (omega) | 1445 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.784 |  | pass | 96 |
| a3pin12_own_L_om1_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (omega) | 1532 (alpha[2]) | 2803 (f[1,22,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.765 |  | fail: divergences | 333 |
| a3pin12_sim_S0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0022 (f[1,1,3]) | 1150 (omega) | 1469 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.703 |  | pass | 53 |
| a3pin12_sim_S0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (f[1,10,1]) | 1344 (lp__) | 2333 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.710 |  | pass | 95 |
| a3pin12_sim_S0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0019 (omega) | 2493 (lp__) | 4608 (lp__) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.685 |  | pass | 267 |
| a3pin12_sim_S1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0096 (sigma) | 594 (lp__) | 1137 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.557 |  | pass | 32 |
| a3pin12_sim_S1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0028 (f[1,24,1]) | 1486 (omega) | 2135 (f[1,2,4]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.783 |  | pass | 51 |
| a3pin12_sim_S1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0063 (lp__) | 1018 (sigma) | 1375 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.666 |  | pass | 46 |
| a3pin12_sim_S2_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (f[1,7,4]) | 915 (sigma) | 514 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.805 |  | pass | 55 |
| a3pin12_sim_S2_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (lp__) | 1050 (sigma) | 1222 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.707 |  | pass | 47 |
| a3pin12_sim_S2_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0037 (f[1,5,1]) | 1337 (sigma) | 1126 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.741 |  | pass | 49 |
| a3pin12_slope_rho30_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (sigma) | 808 (omega) | 831 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.760 |  | pass | 113 |
| a3pin12_slope_rho30_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0035 (sigma) | 1434 (omega) | 2407 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.724 |  | pass | 344 |
| a3pin12_slope_rho30_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0035 (f[1,20,1]) | 1946 (f[1,6,1]) | 3762 (f[1,25,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.693 |  | pass | 323 |
| a3pin12_slope_rho60_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0040 (f[1,24,1]) | 1588 (omega) | 2624 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.720 |  | pass | 277 |
| a3pin12_slope_rho60_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0030 (lp__) | 1623 (f[1,22,1]) | 3220 (f[1,26,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 11) | 0.774 |  | fail: divergences | 321 |
| a3pin12_slope_rho60_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0052 (f[1,20,4]) | 867 (omega) | 595 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.752 |  | pass | 93 |
