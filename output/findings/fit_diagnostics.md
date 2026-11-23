# Per-fit sampler diagnostics for the latent-group probe's Stan fits (rule 16)

Produced by `analyses/probe_latent_groups/export_fit_diagnostics.py` from the fit records in `output/gp_cache/` (gitignored; 283 records with diagnostics). Every value is the one the gate read at the final rung; the gate thresholds are R-hat at most 1.01 on gated variables, bulk and tail ESS at least 400, no divergent transitions, E-BFMI at least 0.3 and, for A1, a per-chain co-membership split of at most 0.5. Treedepth saturation is reported, not gated. Rung settings are iterations per chain / warmup, adapt_delta and max_treedepth; four chains throughout. A validation label has two rows, one per class-intercept prior, each fitted to records simulated at that prior's basin-fit values; the findings tables named per family report the Normal(0, 10) fits under the plain name and the Normal(0, 2) fits (the validation scripts run with `--alpha-sd=2`) as `*_alphaN0-2.md`. Descriptive export; exploratory design validation, not evidence about the record. Columns in `fit_diagnostics.csv`.

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

## A1 latent-unit recovery and equivalence fits (validate_a1.R; a1_validation.md and a1_validation_alphaN0-2.md)

34 records, 7 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a1val_a0_on_zero | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0029 (tau[1]) | 1150 (tau[1]) | 1324 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.769 |  | pass | 132 |
| a1val_a0_on_zero | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0051 (tau[1]) | 947 (tau[1]) | 1440 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.779 |  | pass | 99 |
| a1val_a1_fixed0_on_zero | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0030 (sigma) | 1446 (tau[1]) | 2020 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.735 |  | pass | 257 |
| a1val_a1_fixed0_on_zero | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0024 (f[5,4]) | 1165 (tau[1]) | 1806 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.759 |  | pass | 189 |
| a1val_a1_on_zero | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0179 (comembership[5,2]) | 386 (comembership[27,19]) | 821 (comembership[7,5]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.674 | 0.096 | fail: rhat, ess_bulk | 975 |
| a1val_a1_on_zero | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0457 (comembership[25,21]) | 120 (comembership[25,21]) | 262 (comembership[9,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.663 | 0.244 | fail: rhat, ess_bulk, ess_tail | 780 |
| a1val_sg0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0060 (Delta_w) | 533 (comembership[3,2]) | 1411 (comembership[14,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.721 | 0.070 | pass | 917 |
| a1val_sg0.25_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0087 (comembership[22,2]) | 859 (comembership[3,2]) | 1245 (comembership[11,7]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.743 | 0.098 | pass | 731 |
| a1val_sg0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0415 (comembership[20,2]) | 123 (comembership[20,2]) | 292 (Delta_w) | 4/12000 (0.03%) | 0/12000 (0.0%; deepest 10) | 0.667 | 0.119 | fail: rhat, ess_bulk, ess_tail, divergences | 444 |
| a1val_sg0.25_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0283 (comembership[14,3]) | 175 (comembership[25,2]) | 305 (sigma_g) | 6/12000 (0.05%) | 0/12000 (0.0%; deepest 12) | 0.675 | 0.168 | fail: rhat, ess_bulk, ess_tail, divergences | 1290 |
| a1val_sg0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0146 (comembership[26,16]) | 419 (comembership[3,2]) | 763 (comembership[27,9]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.693 | 0.105 | fail: rhat | 966 |
| a1val_sg0.25_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0191 (comembership[20,2]) | 279 (comembership[20,2]) | 670 (comembership[12,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.685 | 0.119 | fail: rhat, ess_bulk | 872 |
| a1val_sg0.25_r4 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0200 (comembership[22,3]) | 380 (comembership[22,3]) | 572 (comembership[25,5]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.634 | 0.117 | fail: rhat, ess_bulk | 891 |
| a1val_sg0.25_r4 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0102 (comembership[19,13]) | 416 (comembership[22,3]) | 701 (comembership[27,24]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.665 | 0.105 | fail: rhat | 842 |
| a1val_sg0.25_r5 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0116 (comembership[22,11]) | 472 (comembership[25,19]) | 1114 (comembership[13,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.664 | 0.116 | fail: rhat | 877 |
| a1val_sg0.25_r5 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0238 (Delta_w) | 266 (Delta_w) | 437 (comembership[17,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.703 | 0.186 | fail: rhat, ess_bulk | 741 |
| a1val_sg0.5_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0106 (comembership[19,8]) | 873 (comembership[19,3]) | 1108 (comembership[8,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.715 | 0.090 | fail: rhat | 976 |
| a1val_sg0.5_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0105 (comembership[21,19]) | 909 (comembership[25,2]) | 1147 (comembership[8,4]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.673 | 0.100 | fail: rhat | 696 |
| a1val_sg0.5_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0154 (comembership[20,2]) | 403 (comembership[13,5]) | 372 (comembership[13,7]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.670 | 0.101 | fail: rhat, ess_tail | 847 |
| a1val_sg0.5_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0165 (comembership[16,3]) | 440 (comembership[13,1]) | 414 (comembership[17,6]) | 10/12000 (0.08%) | 0/12000 (0.0%; deepest 11) | 0.664 | 0.127 | fail: rhat, divergences | 674 |
| a1val_sg0.5_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0389 (comembership[19,5]) | 150 (comembership[5,3]) | 291 (comembership[28,3]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.722 | 0.184 | fail: rhat, ess_bulk, ess_tail | 699 |
| a1val_sg0.5_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0118 (comembership[21,5]) | 362 (comembership[21,5]) | 454 (comembership[15,4]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.671 | 0.092 | fail: rhat, ess_bulk | 438 |
| a1val_sg0.5_r4 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0088 (comembership[22,20]) | 632 (comembership[21,3]) | 1045 (comembership[21,18]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.656 | 0.071 | pass | 1183 |
| a1val_sg0.5_r4 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0246 (comembership[21,8]) | 440 (comembership[21,8]) | 660 (comembership[22,9]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.661 | 0.099 | fail: rhat | 700 |
| a1val_sg0.5_r5 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0110 (comembership[25,22]) | 593 (comembership[25,21]) | 916 (comembership[20,16]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.677 | 0.117 | fail: rhat | 733 |
| a1val_sg0.5_r5 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0167 (comembership[27,25]) | 600 (comembership[27,21]) | 492 (comembership[18,6]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.645 | 0.131 | fail: rhat | 674 |
| a1val_sg1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0900 (comembership[22,16]) | 38 (comembership[22,19]) | 47 (comembership[17,2]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 11) | 0.626 | 0.374 | fail: rhat, ess_bulk, ess_tail, divergences | 610 |
| a1val_sg1_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0415 (comembership[25,19]) | 106 (comembership[25,22]) | 407 (comembership[7,3]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.660 | 0.250 | fail: rhat, ess_bulk | 582 |
| a1val_sg1_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0127 (comembership[25,21]) | 484 (comembership[21,19]) | 992 (comembership[9,3]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.656 | 0.137 | fail: rhat | 415 |
| a1val_sg1_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0411 (comembership[20,11]) | 114 (comembership[25,20]) | 215 (sigma_g) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.666 | 0.204 | fail: rhat, ess_bulk, ess_tail | 662 |
| a1val_sg2_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.6111 (comembership[27,20]) | 7 (comembership[27,20]) | 6 (comembership[20,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.745 | 0.726 | fail: rhat, ess_bulk, ess_tail, comembership_split | 1014 |
| a1val_sg2_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0695 (comembership[21,5]) | 48 (comembership[21,5]) | 187 (comembership[5,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.754 | 0.295 | fail: rhat, ess_bulk, ess_tail | 535 |
| a1val_sg2_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0206 (comembership[21,8]) | 205 (comembership[25,22]) | 548 (comembership[23,8]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.710 | 0.164 | fail: rhat, ess_bulk | 690 |
| a1val_sg2_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0152 (comembership[8,2]) | 319 (comembership[19,2]) | 265 (comembership[8,6]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.737 | 0.131 | fail: rhat, ess_bulk, ess_tail | 556 |

## A3 barrier recovery and equivalence fits, range free (validate_a3.R; a3_validation.md and a3_validation_alphaN0-2.md)

82 records, 82 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3val_a0_on_om0 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0022 (f[4,6]) | 1129 (lp__) | 2042 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.715 |  | pass | 134 |
| a3val_a0_on_om0 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0023 (tau[1]) | 1451 (lp__) | 2647 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.721 |  | pass | 119 |
| a3val_a3_fixed0_on_om0 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0029 (f[1,15,1]) | 1182 (tau[1]) | 1495 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.760 |  | pass | 184 |
| a3val_a3_fixed0_on_om0 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0057 (lp__) | 1047 (lp__) | 1791 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.725 |  | pass | 146 |
| a3val_own_G1_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0024 (rho) | 1270 (rho) | 1691 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.708 |  | pass | 554 |
| a3val_own_G1_om0.25_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0054 (rho) | 1153 (rho) | 1079 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.776 |  | pass | 119 |
| a3val_own_G1_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (sigma) | 1925 (rho) | 3460 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.713 |  | pass | 923 |
| a3val_own_G1_om0.25_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0026 (lp__) | 1079 (tau[1]) | 1838 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.727 |  | pass | 161 |
| a3val_own_G1_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0096 (sigma) | 460 (rho) | 825 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.735 |  | pass | 122 |
| a3val_own_G1_om0.25_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0020 (tau[1]) | 1657 (tau[1]) | 2662 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.742 |  | pass | 556 |
| a3val_own_G1_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0073 (tau[1]) | 636 (rho) | 544 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.729 |  | pass | 192 |
| a3val_own_G1_om0.5_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0048 (omega) | 605 (tau[1]) | 915 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.734 |  | pass | 134 |
| a3val_own_G1_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0073 (rho) | 877 (tau[1]) | 1293 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.714 |  | pass | 151 |
| a3val_own_G1_om0.5_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0052 (omega) | 1059 (tau[1]) | 1382 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.734 |  | pass | 108 |
| a3val_own_G1_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0066 (lp__) | 968 (tau[1]) | 1448 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 190 |
| a3val_own_G1_om0.5_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0035 (tau[1]) | 802 (tau[1]) | 1403 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.744 |  | pass | 153 |
| a3val_own_G1_om0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0035 (sigma) | 1291 (lp__) | 2090 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.738 |  | pass | 179 |
| a3val_own_G1_om0_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0014 (rho) | 2327 (omega) | 3782 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.731 |  | pass | 558 |
| a3val_own_G1_om0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0055 (sigma) | 1210 (tau[1]) | 1494 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.764 |  | pass | 190 |
| a3val_own_G1_om0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0023 (log_lik[10]) | 1088 (tau[1]) | 1443 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.784 |  | pass | 146 |
| a3val_own_G1_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0039 (rho) | 1331 (rho) | 1880 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.679 |  | pass | 551 |
| a3val_own_G1_om0_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0052 (lp__) | 875 (tau[1]) | 1190 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.734 |  | pass | 147 |
| a3val_own_G1_om1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0037 (rho) | 1694 (rho) | 2862 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.793 |  | pass | 384 |
| a3val_own_G1_om1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0049 (tau[1]) | 532 (tau[1]) | 720 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.782 |  | pass | 155 |
| a3val_own_G1_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0032 (omega) | 1034 (rho) | 1773 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.735 |  | pass | 166 |
| a3val_own_G1_om1_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0024 (omega) | 1814 (tau[1]) | 2322 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.716 |  | pass | 557 |
| a3val_own_G1_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0074 (lp__) | 1150 (tau[1]) | 1304 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.799 |  | pass | 167 |
| a3val_own_G1_om1_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0075 (tau[1]) | 600 (tau[1]) | 870 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.821 |  | pass | 113 |
| a3val_own_L_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0084 (f[1,17,1]) | 577 (tau[1]) | 1168 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 134 |
| a3val_own_L_om0.25_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0061 (tau[1]) | 927 (tau[1]) | 1508 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.778 |  | pass | 163 |
| a3val_own_L_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0048 (rho) | 1329 (rho) | 1558 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.755 |  | pass | 151 |
| a3val_own_L_om0.25_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0043 (lp__) | 810 (tau[1]) | 1148 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.723 |  | pass | 146 |
| a3val_own_L_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0025 (sigma) | 1617 (rho) | 1934 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.748 |  | pass | 356 |
| a3val_own_L_om0.25_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0056 (rho) | 610 (tau[1]) | 887 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.687 |  | pass | 115 |
| a3val_own_L_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0049 (f[1,18,1]) | 1141 (rho) | 1400 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.780 |  | pass | 161 |
| a3val_own_L_om0.5_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0039 (rho) | 670 (rho) | 829 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.762 |  | pass | 477 |
| a3val_own_L_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0052 (f[1,15,1]) | 1029 (rho) | 1719 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.726 |  | pass | 152 |
| a3val_own_L_om0.5_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0039 (tau[1]) | 897 (tau[1]) | 1413 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.777 |  | pass | 184 |
| a3val_own_L_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0054 (tau[1]) | 1261 (lp__) | 2090 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.717 |  | pass | 148 |
| a3val_own_L_om0.5_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0047 (omega) | 1052 (omega) | 1244 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.696 |  | pass | 117 |
| a3val_own_L_om0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0071 (lp__) | 1346 (tau[1]) | 1680 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.734 |  | pass | 155 |
| a3val_own_L_om0_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0055 (omega) | 879 (omega) | 916 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.743 |  | pass | 105 |
| a3val_own_L_om0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0034 (lp__) | 1357 (lp__) | 2304 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.726 |  | pass | 175 |
| a3val_own_L_om0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0032 (lp__) | 1285 (lp__) | 1712 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.754 |  | pass | 156 |
| a3val_own_L_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0020 (f[1,5,2]) | 1873 (rho) | 3117 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.757 |  | pass | 503 |
| a3val_own_L_om0_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0035 (tau[1]) | 934 (tau[1]) | 1721 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.750 |  | pass | 106 |
| a3val_own_L_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0047 (sigma) | 1042 (rho) | 1001 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.759 |  | pass | 150 |
| a3val_own_L_om1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0056 (tau[1]) | 775 (omega) | 703 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.717 |  | pass | 149 |
| a3val_own_L_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0035 (f[1,13,3]) | 1566 (omega) | 2771 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.778 |  | pass | 155 |
| a3val_own_L_om1_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0078 (omega) | 732 (omega) | 1849 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.710 |  | pass | 137 |
| a3val_own_L_om1_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0030 (f[1,14,2]) | 1186 (rho) | 1825 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.788 |  | pass | 413 |
| a3val_own_L_om1_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0033 (f[1,8,6]) | 624 (tau[1]) | 1102 (rho) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.742 |  | pass | 154 |
| a3val_sim_S0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0011 (lp__) | 2665 (omega) | 3133 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.753 |  | pass | 268 |
| a3val_sim_S0_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0049 (omega) | 892 (omega) | 1481 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.682 |  | pass | 491 |
| a3val_sim_S0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0029 (lp__) | 1526 (lp__) | 2325 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.737 |  | pass | 165 |
| a3val_sim_S0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0047 (sigma) | 837 (omega) | 1005 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.730 |  | pass | 121 |
| a3val_sim_S0_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0026 (f[1,16,2]) | 1404 (lp__) | 1767 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.699 |  | pass | 94 |
| a3val_sim_S0_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0017 (omega) | 990 (omega) | 3073 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.710 |  | pass | 558 |
| a3val_sim_S1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0015 (omega) | 1244 (sigma) | 1541 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.571 |  | pass | 274 |
| a3val_sim_S1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0021 (alpha[9]) | 864 (sigma) | 991 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.602 |  | pass | 92 |
| a3val_sim_S1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0025 (f[1,18,2]) | 1414 (sigma) | 1843 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.798 |  | pass | 83 |
| a3val_sim_S1_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0021 (f[1,8,7]) | 987 (omega) | 2492 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.718 |  | pass | 345 |
| a3val_sim_S1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0085 (lp__) | 615 (sigma) | 1115 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.532 |  | pass | 89 |
| a3val_sim_S1_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0044 (omega) | 1002 (sigma) | 1394 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.686 |  | pass | 94 |
| a3val_sim_S2_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0047 (sigma) | 1266 (sigma) | 2639 (f[1,10,5]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.784 |  | pass | 83 |
| a3val_sim_S2_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0060 (sigma) | 662 (sigma) | 2091 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.744 |  | pass | 441 |
| a3val_sim_S2_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0035 (sigma) | 1280 (sigma) | 2085 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.662 |  | pass | 256 |
| a3val_sim_S2_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0042 (sigma) | 1010 (omega) | 1363 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.650 |  | pass | 81 |
| a3val_sim_S2_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0046 (tau[1]) | 791 (sigma) | 1153 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.757 |  | pass | 47 |
| a3val_sim_S2_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0051 (lp__) | 714 (sigma) | 667 (sigma) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.651 |  | pass | 324 |
| a3val_slope_rho30_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0019 (rho) | 1411 (rho) | 2753 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.772 |  | pass | 485 |
| a3val_slope_rho30_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0032 (tau[1]) | 1341 (rho) | 2650 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.748 |  | pass | 498 |
| a3val_slope_rho30_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0026 (f[1,10,4]) | 1516 (sigma) | 2432 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.745 |  | pass | 432 |
| a3val_slope_rho30_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0013 (rho) | 1043 (rho) | 2233 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.742 |  | pass | 467 |
| a3val_slope_rho30_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0024 (f[1,24,2]) | 1609 (rho) | 3241 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.742 |  | pass | 553 |
| a3val_slope_rho30_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0014 (tau[1]) | 2226 (rho) | 3014 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.748 |  | pass | 457 |
| a3val_slope_rho60_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0044 (f[1,13,4]) | 720 (tau[1]) | 1232 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.727 |  | pass | 118 |
| a3val_slope_rho60_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0023 (lp__) | 1898 (tau[1]) | 3121 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.720 |  | pass | 461 |
| a3val_slope_rho60_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0023 (f[1,4,7]) | 1415 (lp__) | 2228 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.752 |  | pass | 142 |
| a3val_slope_rho60_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0018 (f[1,13,3]) | 2740 (lp__) | 4200 (tau[1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.751 |  | pass | 486 |
| a3val_slope_rho60_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0054 (tau[1]) | 1237 (tau[1]) | 2235 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.720 |  | pass | 555 |
| a3val_slope_rho60_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0032 (tau[1]) | 1117 (tau[1]) | 1579 (rho) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.757 |  | pass | 456 |

## A3 with the range pinned at 6 km (a3_pinned_range.R 6; a3_pinned_rho6.md and _alphaN0-2.md)

78 records, 78 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3pin6_own_G1_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0092 (omega) | 837 (omega) | 1261 (f[1,20,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.692 |  | pass | 80 |
| a3pin6_own_G1_om0.25_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0064 (sigma) | 1138 (sigma) | 1514 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.747 |  | pass | 102 |
| a3pin6_own_G1_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0038 (lp__) | 1250 (lp__) | 1249 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.689 |  | pass | 86 |
| a3pin6_own_G1_om0.25_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0071 (omega) | 1043 (omega) | 1177 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.721 |  | pass | 107 |
| a3pin6_own_G1_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0050 (f[1,11,1]) | 822 (f[1,20,1]) | 1562 (f[1,17,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.760 |  | pass | 75 |
| a3pin6_own_G1_om0.25_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0064 (lp__) | 912 (sigma) | 1641 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.678 |  | pass | 100 |
| a3pin6_own_G1_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0053 (f[1,8,1]) | 1166 (sigma) | 1260 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.742 |  | pass | 95 |
| a3pin6_own_G1_om0.5_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0054 (sigma) | 886 (tau[1]) | 926 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.722 |  | pass | 95 |
| a3pin6_own_G1_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0067 (omega) | 1140 (omega) | 1826 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.721 |  | pass | 84 |
| a3pin6_own_G1_om0.5_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0027 (f[1,25,2]) | 1260 (lp__) | 1636 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.756 |  | pass | 96 |
| a3pin6_own_G1_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0075 (lp__) | 955 (lp__) | 1439 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.697 |  | pass | 101 |
| a3pin6_own_G1_om0.5_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0070 (tau[1]) | 901 (sigma) | 1145 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.763 |  | pass | 102 |
| a3pin6_own_G1_om0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0065 (omega) | 556 (f[1,22,1]) | 1211 (f[1,20,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.755 |  | pass | 82 |
| a3pin6_own_G1_om0_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0034 (sigma) | 1088 (sigma) | 1729 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.748 |  | pass | 95 |
| a3pin6_own_G1_om0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0078 (f[1,12,1]) | 603 (f[1,3,1]) | 952 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.730 |  | pass | 73 |
| a3pin6_own_G1_om0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0050 (f[1,14,1]) | 810 (sigma) | 1573 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.732 |  | pass | 99 |
| a3pin6_own_G1_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0056 (omega) | 1094 (omega) | 1719 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.728 |  | pass | 202 |
| a3pin6_own_G1_om0_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0040 (omega) | 871 (omega) | 819 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.663 |  | pass | 526 |
| a3pin6_own_G1_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0039 (f[1,26,1]) | 1204 (f[1,20,1]) | 1884 (f[1,27,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.776 |  | pass | 86 |
| a3pin6_own_G1_om1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0034 (f[1,7,1]) | 1057 (sigma) | 1337 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.739 |  | pass | 104 |
| a3pin6_own_G1_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0054 (lp__) | 1182 (lp__) | 2352 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.736 |  | pass | 96 |
| a3pin6_own_G1_om1_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0063 (sigma) | 926 (lp__) | 1686 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.680 |  | pass | 100 |
| a3pin6_own_G1_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0042 (sigma) | 858 (omega) | 1120 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.723 |  | pass | 102 |
| a3pin6_own_G1_om1_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0044 (sigma) | 1007 (lp__) | 1884 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 142 |
| a3pin6_own_L_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (f[1,20,1]) | 1076 (alpha[1]) | 1978 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.702 |  | pass | 203 |
| a3pin6_own_L_om0.25_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0062 (omega) | 811 (omega) | 724 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.683 |  | pass | 380 |
| a3pin6_own_L_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0037 (f[1,18,2]) | 1499 (omega) | 2324 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.752 |  | pass | 293 |
| a3pin6_own_L_om0.25_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0048 (lp__) | 875 (omega) | 967 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.682 |  | pass | 572 |
| a3pin6_own_L_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0040 (f[1,11,1]) | 1517 (f[1,3,1]) | 3007 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.737 |  | pass | 295 |
| a3pin6_own_L_om0.25_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0029 (tau[1]) | 1182 (omega) | 1253 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.678 |  | pass | 90 |
| a3pin6_own_L_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0066 (f[1,21,2]) | 726 (f[1,27,1]) | 1566 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.761 |  | pass | 104 |
| a3pin6_own_L_om0.5_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0041 (omega) | 697 (sigma) | 1680 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.725 |  | pass | 104 |
| a3pin6_own_L_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0079 (f[1,8,1]) | 649 (alpha[1]) | 999 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.668 |  | pass | 74 |
| a3pin6_own_L_om0.5_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0034 (f[1,22,2]) | 1017 (sigma) | 1349 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.701 |  | pass | 101 |
| a3pin6_own_L_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0081 (f[1,17,4]) | 894 (f[1,27,2]) | 1233 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.717 |  | pass | 85 |
| a3pin6_own_L_om0.5_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0067 (f[1,4,1]) | 1085 (omega) | 1556 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.734 |  | pass | 100 |
| a3pin6_own_L_om0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0031 (f[1,23,1]) | 1488 (alpha[1]) | 2217 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.747 |  | pass | 268 |
| a3pin6_own_L_om0_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0029 (f[1,28,7]) | 981 (sigma) | 1639 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.701 |  | pass | 102 |
| a3pin6_own_L_om0_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0031 (f[1,26,2]) | 1225 (omega) | 1609 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 9) | 0.691 |  | pass | 212 |
| a3pin6_own_L_om0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0025 (omega) | 1151 (omega) | 887 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.722 |  | pass | 102 |
| a3pin6_own_L_om0_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0051 (omega) | 635 (omega) | 1019 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.702 |  | pass | 52 |
| a3pin6_own_L_om0_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0046 (omega) | 952 (omega) | 779 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.643 |  | pass | 95 |
| a3pin6_own_L_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0086 (f[1,18,1]) | 745 (omega) | 1354 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.725 |  | pass | 116 |
| a3pin6_own_L_om1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0044 (sigma) | 1055 (sigma) | 1617 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.772 |  | pass | 100 |
| a3pin6_own_L_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0040 (f[1,25,1]) | 1150 (omega) | 1934 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.752 |  | pass | 102 |
| a3pin6_own_L_om1_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0047 (lp__) | 968 (sigma) | 1937 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.771 |  | pass | 101 |
| a3pin6_own_L_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0094 (sigma) | 596 (omega) | 845 (f[1,12,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.715 |  | pass | 99 |
| a3pin6_own_L_om1_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0051 (sigma) | 702 (sigma) | 611 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.738 |  | pass | 113 |
| a3pin6_sim_S0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0032 (f[1,19,7]) | 1309 (lp__) | 2239 (f[1,17,5]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.744 |  | pass | 60 |
| a3pin6_sim_S0_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0086 (lp__) | 1253 (omega) | 1501 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.734 |  | pass | 107 |
| a3pin6_sim_S0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0037 (f[1,12,2]) | 1266 (f[1,22,1]) | 2100 (f[1,8,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.679 |  | pass | 100 |
| a3pin6_sim_S0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0054 (lp__) | 1341 (f[1,5,1]) | 1833 (f[1,27,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.716 |  | pass | 112 |
| a3pin6_sim_S0_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0049 (lp__) | 999 (lp__) | 1898 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.609 |  | pass | 91 |
| a3pin6_sim_S0_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0041 (lp__) | 1086 (sigma) | 1258 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.593 |  | pass | 155 |
| a3pin6_sim_S1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0057 (omega) | 742 (lp__) | 1524 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.537 |  | pass | 40 |
| a3pin6_sim_S1_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0041 (sigma) | 1532 (lp__) | 2612 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 9) | 0.555 |  | pass | 813 |
| a3pin6_sim_S1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0053 (f[1,23,1]) | 1538 (alpha[1]) | 1987 (f[1,8,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.771 |  | pass | 66 |
| a3pin6_sim_S1_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0064 (omega) | 675 (omega) | 1655 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.664 |  | pass | 195 |
| a3pin6_sim_S1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0072 (lp__) | 1022 (lp__) | 1856 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.634 |  | pass | 54 |
| a3pin6_sim_S1_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0036 (omega) | 1172 (omega) | 2866 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.577 |  | pass | 580 |
| a3pin6_sim_S2_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0023 (f[1,23,5]) | 1782 (f[1,20,2]) | 1561 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.770 |  | pass | 73 |
| a3pin6_sim_S2_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0077 (f[1,28,4]) | 863 (omega) | 634 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.770 |  | pass | 123 |
| a3pin6_sim_S2_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0044 (sigma) | 1088 (sigma) | 1638 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.664 |  | pass | 53 |
| a3pin6_sim_S2_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0072 (lp__) | 487 (omega) | 1294 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.559 |  | pass | 107 |
| a3pin6_sim_S2_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0035 (f[1,26,4]) | 1115 (sigma) | 584 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.711 |  | pass | 56 |
| a3pin6_sim_S2_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0046 (omega) | 619 (omega) | 1250 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.687 |  | pass | 103 |
| a3pin6_slope_rho30_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0044 (f[1,8,4]) | 905 (alpha[1]) | 1386 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.756 |  | pass | 101 |
| a3pin6_slope_rho30_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0072 (sigma) | 930 (sigma) | 1785 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.719 |  | pass | 101 |
| a3pin6_slope_rho30_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0050 (f[1,21,4]) | 928 (f[1,22,1]) | 1745 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.740 |  | pass | 87 |
| a3pin6_slope_rho30_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0021 (omega) | 1199 (lp__) | 1562 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.727 |  | pass | 97 |
| a3pin6_slope_rho30_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0020 (lp__) | 1080 (f[1,22,1]) | 2300 (f[1,13,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.753 |  | pass | 309 |
| a3pin6_slope_rho30_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0069 (sigma) | 834 (omega) | 727 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.719 |  | pass | 104 |
| a3pin6_slope_rho60_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0036 (f[1,9,1]) | 1308 (lp__) | 2301 (f[1,12,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.713 |  | pass | 89 |
| a3pin6_slope_rho60_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0057 (sigma) | 1019 (lp__) | 1960 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.687 |  | pass | 109 |
| a3pin6_slope_rho60_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0042 (f[1,8,1]) | 981 (omega) | 1171 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.767 |  | pass | 86 |
| a3pin6_slope_rho60_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0039 (omega) | 878 (omega) | 1025 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.766 |  | pass | 107 |
| a3pin6_slope_rho60_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0096 (lp__) | 1205 (lp__) | 1711 (f[1,6,2]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.657 |  | pass | 94 |
| a3pin6_slope_rho60_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0038 (f[1,6,2]) | 1023 (lp__) | 1962 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.745 |  | pass | 110 |

## A3 with the range pinned at 12 km (a3_pinned_range.R 12; a3_pinned_rho12.md and _alphaN0-2.md)

78 records, 70 pass the gate.

| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a3pin12_own_G1_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0026 (f[1,15,2]) | 1446 (omega) | 2143 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.743 |  | pass | 273 |
| a3pin12_own_G1_om0.25_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0032 (f[1,14,1]) | 902 (sigma) | 1130 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.645 |  | pass | 112 |
| a3pin12_own_G1_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0018 (lp__) | 2199 (omega) | 4042 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.702 |  | pass | 346 |
| a3pin12_own_G1_om0.25_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0022 (f[1,16,8]) | 1126 (sigma) | 1418 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.669 |  | pass | 102 |
| a3pin12_own_G1_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0031 (omega) | 1407 (lp__) | 1885 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.758 |  | pass | 81 |
| a3pin12_own_G1_om0.25_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0030 (f[1,17,3]) | 990 (sigma) | 1668 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.739 |  | pass | 107 |
| a3pin12_own_G1_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0027 (f[1,10,2]) | 520 (tau[1]) | 797 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.739 |  | pass | 93 |
| a3pin12_own_G1_om0.5_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0055 (sigma) | 495 (tau[1]) | 791 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.708 |  | pass | 101 |
| a3pin12_own_G1_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0059 (omega) | 958 (omega) | 1445 (f[1,19,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.719 |  | pass | 80 |
| a3pin12_own_G1_om0.5_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0075 (lp__) | 920 (lp__) | 1499 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.727 |  | pass | 94 |
| a3pin12_own_G1_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0046 (tau[1]) | 964 (tau[1]) | 1366 (tau[1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.738 |  | pass | 100 |
| a3pin12_own_G1_om0.5_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0049 (lp__) | 779 (sigma) | 1445 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.718 |  | pass | 104 |
| a3pin12_own_G1_om0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0018 (lp__) | 1916 (f[1,26,1]) | 3007 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.743 |  | pass | 339 |
| a3pin12_own_G1_om0_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0059 (sigma) | 1025 (sigma) | 1192 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.725 |  | pass | 96 |
| a3pin12_own_G1_om0_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0061 (f[1,23,1]) | 1434 (f[1,20,1]) | 2798 (f[1,21,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.730 |  | pass | 258 |
| a3pin12_own_G1_om0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0065 (omega) | 974 (omega) | 973 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.716 |  | pass | 111 |
| a3pin12_own_G1_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (tau[1]) | 1591 (omega) | 3270 (f[1,27,1]) | 2/12000 (0.02%) | 0/12000 (0.0%; deepest 10) | 0.706 |  | fail: divergences | 224 |
| a3pin12_own_G1_om0_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0029 (sigma) | 1031 (omega) | 965 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.690 |  | pass | 449 |
| a3pin12_own_G1_om1_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0025 (f[1,1,1]) | 1742 (f[1,20,1]) | 2313 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.704 |  | pass | 321 |
| a3pin12_own_G1_om1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0093 (tau[1]) | 662 (tau[1]) | 1051 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.694 |  | pass | 103 |
| a3pin12_own_G1_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0046 (lp__) | 1308 (lp__) | 1817 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.753 |  | pass | 104 |
| a3pin12_own_G1_om1_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0063 (sigma) | 905 (lp__) | 1417 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.654 |  | pass | 110 |
| a3pin12_own_G1_om1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0089 (tau[1]) | 891 (tau[1]) | 1788 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.742 |  | pass | 161 |
| a3pin12_own_G1_om1_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0061 (tau[1]) | 439 (tau[1]) | 774 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.798 |  | pass | 132 |
| a3pin12_own_L_om0.25_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0016 (f[1,11,1]) | 1078 (alpha[1]) | 2159 (f[1,19,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 9) | 0.732 |  | fail: divergences | 205 |
| a3pin12_own_L_om0.25_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0042 (omega) | 1106 (omega) | 1264 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.661 |  | pass | 548 |
| a3pin12_own_L_om0.25_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0015 (tau[1]) | 2072 (alpha[1]) | 4038 (f[1,20,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.721 |  | fail: divergences | 270 |
| a3pin12_own_L_om0.25_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0040 (omega) | 498 (omega) | 598 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.706 |  | pass | 141 |
| a3pin12_own_L_om0.25_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0023 (f[1,3,1]) | 1902 (alpha[1]) | 2664 (f[1,20,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.769 |  | pass | 244 |
| a3pin12_own_L_om0.25_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0051 (sigma) | 764 (omega) | 1088 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.722 |  | pass | 93 |
| a3pin12_own_L_om0.5_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0034 (omega) | 929 (f[1,20,1]) | 2269 (f[1,10,1]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.810 |  | pass | 98 |
| a3pin12_own_L_om0.5_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0055 (lp__) | 899 (sigma) | 1887 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.745 |  | pass | 131 |
| a3pin12_own_L_om0.5_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0026 (omega) | 1772 (omega) | 3084 (f[1,25,2]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.729 |  | pass | 281 |
| a3pin12_own_L_om0.5_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0080 (omega) | 897 (omega) | 831 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.745 |  | pass | 128 |
| a3pin12_own_L_om0.5_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0054 (omega) | 945 (f[1,22,1]) | 871 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.743 |  | pass | 75 |
| a3pin12_own_L_om0.5_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0047 (omega) | 1042 (omega) | 1203 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.733 |  | pass | 100 |
| a3pin12_own_L_om0_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0031 (f[1,3,2]) | 1839 (f[1,20,1]) | 3691 (alpha[2]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.722 |  | fail: divergences | 294 |
| a3pin12_own_L_om0_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0030 (sigma) | 1389 (omega) | 1236 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.750 |  | pass | 142 |
| a3pin12_own_L_om0_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0025 (omega) | 1833 (omega) | 3941 (f[1,21,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 11) | 0.729 |  | fail: divergences | 346 |
| a3pin12_own_L_om0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0045 (sigma) | 1182 (lp__) | 1231 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.753 |  | pass | 101 |
| a3pin12_own_L_om0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0036 (omega) | 1778 (omega) | 3691 (f[1,22,2]) | 2/12000 (0.02%) | 0/12000 (0.0%; deepest 10) | 0.711 |  | fail: divergences | 272 |
| a3pin12_own_L_om0_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0016 (f[1,13,3]) | 2053 (lp__) | 2501 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.727 |  | pass | 538 |
| a3pin12_own_L_om1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0068 (f[1,1,7]) | 711 (omega) | 1184 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.755 |  | pass | 99 |
| a3pin12_own_L_om1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0029 (f[1,12,5]) | 1048 (sigma) | 1944 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.753 |  | pass | 107 |
| a3pin12_own_L_om1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0030 (f[1,25,1]) | 974 (omega) | 1445 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.784 |  | pass | 96 |
| a3pin12_own_L_om1_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0032 (omega) | 948 (lp__) | 1185 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.703 |  | pass | 104 |
| a3pin12_own_L_om1_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0034 (omega) | 1532 (alpha[2]) | 2803 (f[1,22,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 10) | 0.765 |  | fail: divergences | 333 |
| a3pin12_own_L_om1_r3 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0024 (omega) | 1729 (sigma) | 1348 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 12) | 0.691 |  | pass | 507 |
| a3pin12_sim_S0_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0022 (f[1,1,3]) | 1150 (omega) | 1469 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.703 |  | pass | 53 |
| a3pin12_sim_S0_r1 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0038 (f[1,10,3]) | 2618 (lp__) | 2005 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.722 |  | pass | 570 |
| a3pin12_sim_S0_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (f[1,10,1]) | 1344 (lp__) | 2333 (lp__) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.710 |  | pass | 95 |
| a3pin12_sim_S0_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0064 (lp__) | 1288 (lp__) | 1878 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.658 |  | pass | 159 |
| a3pin12_sim_S0_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0019 (omega) | 2493 (lp__) | 4608 (lp__) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.685 |  | pass | 267 |
| a3pin12_sim_S0_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0032 (f[1,1,4]) | 1062 (lp__) | 1286 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.665 |  | pass | 163 |
| a3pin12_sim_S1_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0096 (sigma) | 594 (lp__) | 1137 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.557 |  | pass | 32 |
| a3pin12_sim_S1_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0067 (lp__) | 528 (lp__) | 1143 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.589 |  | pass | 35 |
| a3pin12_sim_S1_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0028 (f[1,24,1]) | 1486 (omega) | 2135 (f[1,2,4]) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.783 |  | pass | 51 |
| a3pin12_sim_S1_r2 | Normal(0, 2) | 2026-11-23 | 2 (6000/3000, 0.99, 13) | 1.0024 (omega) | 1163 (omega) | 2804 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 11) | 0.691 |  | pass | 365 |
| a3pin12_sim_S1_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0063 (lp__) | 1018 (sigma) | 1375 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.666 |  | pass | 46 |
| a3pin12_sim_S1_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0083 (sigma) | 574 (omega) | 1190 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.569 |  | pass | 78 |
| a3pin12_sim_S2_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (f[1,7,4]) | 915 (sigma) | 514 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.805 |  | pass | 55 |
| a3pin12_sim_S2_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0034 (alpha[9]) | 1225 (omega) | 1281 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.769 |  | pass | 93 |
| a3pin12_sim_S2_r2 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (lp__) | 1050 (sigma) | 1222 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.707 |  | pass | 47 |
| a3pin12_sim_S2_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0066 (f[1,12,5]) | 501 (omega) | 1275 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.567 |  | pass | 84 |
| a3pin12_sim_S2_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0037 (f[1,5,1]) | 1337 (sigma) | 1126 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 8) | 0.741 |  | pass | 49 |
| a3pin12_sim_S2_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0069 (lp__) | 607 (omega) | 1398 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.680 |  | pass | 91 |
| a3pin12_slope_rho30_r1 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0043 (sigma) | 808 (omega) | 831 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.760 |  | pass | 113 |
| a3pin12_slope_rho30_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0023 (tau[1]) | 856 (sigma) | 1646 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 10) | 0.770 |  | pass | 99 |
| a3pin12_slope_rho30_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0035 (sigma) | 1434 (omega) | 2407 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.724 |  | pass | 344 |
| a3pin12_slope_rho30_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0024 (omega) | 1315 (omega) | 1596 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.783 |  | pass | 95 |
| a3pin12_slope_rho30_r3 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0035 (f[1,20,1]) | 1946 (f[1,6,1]) | 3762 (f[1,25,1]) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.693 |  | pass | 323 |
| a3pin12_slope_rho30_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0072 (omega) | 578 (omega) | 919 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 11) | 0.698 |  | pass | 129 |
| a3pin12_slope_rho60_r1 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0040 (f[1,24,1]) | 1588 (omega) | 2624 (omega) | 0/12000 (0.00%) | 0/12000 (0.0%; deepest 10) | 0.720 |  | pass | 277 |
| a3pin12_slope_rho60_r1 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0037 (f[1,23,2]) | 830 (tau[1]) | 1343 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.735 |  | pass | 94 |
| a3pin12_slope_rho60_r2 | Normal(0, 10) | 2026-11-22 | 2 (6000/3000, 0.99, 13) | 1.0030 (lp__) | 1623 (f[1,22,1]) | 3220 (f[1,26,1]) | 1/12000 (0.01%) | 0/12000 (0.0%; deepest 11) | 0.774 |  | fail: divergences | 321 |
| a3pin12_slope_rho60_r2 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0055 (f[1,3,1]) | 1056 (omega) | 994 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.726 |  | pass | 87 |
| a3pin12_slope_rho60_r3 | Normal(0, 10) | 2026-11-22 | 1 (3000/1500, 0.95, 12) | 1.0052 (f[1,20,4]) | 867 (omega) | 595 (omega) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.752 |  | pass | 93 |
| a3pin12_slope_rho60_r3 | Normal(0, 2) | 2026-11-23 | 1 (3000/1500, 0.95, 12) | 1.0042 (lp__) | 1099 (sigma) | 1930 (sigma) | 0/6000 (0.00%) | 0/6000 (0.0%; deepest 9) | 0.753 |  | pass | 86 |
