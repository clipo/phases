# A1 validated on its own generator (pilot step 6; spec 4.1c)

Produced by `analyses/probe_latent_groups/validate_a1.R` on 2026-10-17. Records from `a1_contracts.py` (nuisances pinned to the A0 refit's identified quantities; field centered per class; membership = Phillips's partition 12/11/5; Gabriel edges for coherence). Exploratory design validation, not evidence about the record.

## (i) Code equivalence at zero group scale

A0 and A1 with the group scale fixed at zero, same record, same seed. Largest |z| over 41 shared posterior means (difference over combined MCSE): **2.27** at `alpha[8]`. Values above 3 would mean the two programs disagree.

| quantity | A0 mean | A1 (fixed 0) mean | z |
|---|---|---|---|
| sigma | 4.1438 | 4.1109 | -1.53 |
| tau[1] | 0.7105 | 0.7132 | 0.43 |
| rho | 190.9210 | 187.7594 | -1.77 |
| alpha[1] | -0.1912 | -0.1819 | 0.17 |
| alpha[2] | -0.3567 | -0.3681 | -0.25 |
| alpha[3] | -1.3004 | -1.2975 | 0.07 |
| alpha[4] | -0.5395 | -0.5740 | -0.79 |
| alpha[5] | -1.4412 | -1.4667 | -0.65 |
| alpha[6] | -1.4279 | -1.5143 | -2.24 |
| alpha[7] | -1.3603 | -1.3262 | 0.87 |
| alpha[8] | -1.4783 | -1.5756 | -2.27 |
| alpha[9] | -1.7000 | -1.7381 | -1.08 |
| spatial_share | 0.9680 | 0.9671 | -1.22 |
| log_lik[1..28] | | | max abs z 2.20 |

- A0 on zero: R-hat max 1.0041 (f[5,3]) | bulk ESS min 960 (tau[1]) | tail ESS min 1729 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.746 | stuck 0
- A1 fixed 0 on zero: R-hat max 1.0056 (lp__) | bulk ESS min 1069 (tau[1]) | tail ESS min 1940 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.763 | stuck 0

## (ii) Null calibration: unconstrained A1 on the zero-shift record

Truth: sigma_g = 0, no groups; the generating labels (Phillips's partition) are used only to read the margins. Realized between-'group' F_ST of the record from the field alone: 0.3100 (the record's is 0.0123; see the A0 adequacy note below).

| quantity | median [5%, 95%] |
|---|---|
| sigma_g | 0.517 [0.048, 1.386] |
| Delta | 0.975 [0.093, 2.190] |
| Delta_w | 0.910 [0.057, 1.971] |
| coherence (soft) | 0.003 [-0.062, 0.095] |
| coherence (hard) | 0.006 [-0.090, 0.116] |
| separation margins vs generating labels (12, 11, 5) | 0.007, -0.019, -0.005 |
| per-chain co-membership split (max abs diff) | 0.144 |

- A1 on zero: R-hat max 1.0244 (comembership[27,25]) | bulk ESS min 277 (comembership[25,21]) | tail ESS min 450 (comembership[14,11]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.632 | stuck 0
Visual set: `a1_null_rank.png`, `a1_null_trace.png`, `a1_null_pairs.png`, `a1_null_energy.png`.

## (iii) Recovery against the prior (sigma_g = 0.25 and 0.5, five records each)

| record | true sigma_g | realized between-group F_ST | Delta true | sigma_g post | Delta post | Delta_w post | coherence | margins (12, 11, 5) | gate | chain split |
|---|---|---|---|---|---|---|---|---|---|---|
| sg0.25_r1 | 0.25 | 0.2194 | 0.714 | 0.394 [0.040, 1.002] | 0.770 [0.076, 1.692] | 0.735 [0.060, 1.579] | -0.000 [-0.065, 0.068] | -0.00, 0.00, 0.00 | pass | 0.074 |
| sg0.25_r2 | 0.25 | 0.1366 | 0.540 | 0.323 [0.030, 3.846] | 0.613 [0.057, 3.856] | 0.447 [0.000, 1.418] | 0.000 [-0.067, 0.073] | 0.00, -0.00, 0.00 | fail: rhat,ess_bulk,ess_tail | 0.112 |
| sg0.25_r3 | 0.25 | 0.1305 | 0.742 | 0.855 [0.390, 1.643] | 1.442 [0.724, 2.334] | 1.334 [0.707, 1.992] | 0.037 [-0.047, 0.101] | 0.10, -0.09, -0.10 | fail: rhat,ess_bulk,ess_tail | 0.247 |
| sg0.25_r4 | 0.25 | 0.0192 | 0.561 | 0.324 [0.034, 1.217] | 0.632 [0.065, 1.605] | 0.555 [0.015, 1.273] | 0.000 [-0.071, 0.083] | -0.00, 0.01, -0.01 | fail: rhat,ess_bulk | 0.163 |
| sg0.25_r5 | 0.25 | 0.0760 | 0.544 | 0.410 [0.049, 1.153] | 0.778 [0.094, 1.814] | 0.717 [0.045, 1.586] | -0.000 [-0.077, 0.082] | 0.01, -0.01, -0.02 | fail: rhat,ess_bulk | 0.165 |
| sg0.5_r1 | 0.50 | 0.0435 | 1.298 | 0.272 [0.029, 0.924] | 0.522 [0.055, 1.374] | 0.479 [0.025, 1.146] | 0.000 [-0.069, 0.079] | 0.00, -0.00, -0.00 | fail: divergences | 0.117 |
| sg0.5_r2 | 0.50 | 0.0967 | 1.092 | 0.329 [0.035, 1.262] | 0.630 [0.067, 1.836] | 0.576 [0.037, 1.523] | -0.000 [-0.074, 0.084] | -0.01, 0.00, 0.00 | fail: rhat,ess_bulk,ess_tail | 0.141 |
| sg0.5_r3 | 0.50 | 0.5753 | 1.583 | 0.431 [0.042, 1.585] | 0.819 [0.078, 2.249] | 0.766 [0.064, 1.926] | -0.000 [-0.062, 0.070] | 0.01, -0.01, -0.00 | fail: rhat | 0.148 |
| sg0.5_r4 | 0.50 | 0.0301 | 1.024 | 0.484 [0.053, 1.406] | 0.931 [0.102, 2.370] | 0.850 [0.048, 2.152] | -0.004 [-0.112, 0.078] | -0.02, 0.01, -0.02 | fail: rhat,ess_bulk,ess_tail | 0.137 |
| sg0.5_r5 | 0.50 | 0.0488 | 1.153 | 0.268 [0.028, 1.098] | 0.506 [0.054, 1.451] | 0.448 [0.015, 1.135] | 0.000 [-0.066, 0.084] | -0.01, 0.00, 0.01 | pass | 0.077 |
| sg1_r1 | 1.00 | 0.0791 | 2.338 | 0.569 [0.136, 1.096] | 1.154 [0.259, 1.800] | 1.124 [0.197, 1.716] | -0.005 [-0.091, 0.082] | -0.03, 0.01, 0.01 | fail: rhat,ess_bulk | 0.222 |
| sg1_r2 | 1.00 | 0.0496 | 1.636 | 0.400 [0.038, 3.778] | 0.756 [0.076, 3.727] | 0.566 [0.000, 1.563] | 0.000 [-0.084, 0.074] | -0.01, -0.00, 0.00 | fail: rhat,ess_bulk,ess_tail | 0.178 |
| sg2_r1 | 2.00 | 0.3568 | 5.401 | 1.366 [0.155, 2.886] | 2.052 [0.302, 4.615] | 1.716 [0.284, 4.361] | 0.005 [-0.045, 0.315] | 0.12, -0.01, 0.11 | fail: rhat,ess_bulk,ess_tail | 0.452 |
| sg2_r2 | 2.00 | 0.2343 | 3.017 | 0.462 [0.051, 1.402] | 0.885 [0.097, 2.217] | 0.832 [0.083, 2.006] | -0.003 [-0.088, 0.072] | -0.00, 0.00, -0.00 | fail: rhat,ess_bulk | 0.237 |

## A0 adequacy note (found while pinning the generator, 2026-10-17)

The A0 posterior predictive with the fitted field reproduces the record's between-area F_ST (0.0122 against 0.0123), but fresh field draws at the same hyperparameters give a median of 0.085 (10th to 90th percentile 0.02 to 0.26). Between Phillips's Parkin and Kent areas the pooled logit share differs by +0.47 for Parkin_Punctated, -0.12 for Barton/Kent/MPI, -0.90 for Painted, and from -2.2 to +0.83 [corrected 2026-10-19 after the Codex audit of the regional-drift spec: Fortune_Noded is +0.83 in area_logits.md; first written as -0.17 to -2.2, which omitted it] for the rare classes, two of which (Walls_Engraved, Hull_Engraved) are absent from every Parkin-area site; per-class SDs of the logit share across sites are 0.8 to 1.8 (`analyses/probe_latent_groups/area_logits.py`, `output/findings/area_logits.md`). The reference class is not the flattest on this scale. A0 has one field scale for all classes and holds the record's pattern only in the realized field. [Corrected 2026-10-19 (Beau approved): this sentence read "In the record the dominant class, the reference, differs little between areas (Parkin minus Kent, -0.26 logits) while the rare classes differ by 1 to 4 logits". The -0.26 was the Parkin-minus-Kent difference in the site mean of log(Barton/Kent/MPI over Parkin_Punctated), not a property of Parkin_Punctated, and "1 to 4 logits" had no committed procedure and does not reproduce.] [Edited by hand 2026-10-19: regenerating this file would refit the thirteen A1 fits under the D25 contracts, which is the deferred A1 check.] A per-class field scale is the extension this points to. The own-generator records above therefore carry far more between-area contrast than the record (realized F_ST column), which is acceptable for a test of the software and of identification under the model (rule 20c) but is not a resemblance claim.

## Timing and gates

| fit | rung | seconds | gate | reasons | cached |
|---|---|---|---|---|---|
| A0 on zero | 1 | 58 | pass |  | yes |
| A1 fixed 0 on zero | 1 | 165 | pass |  | yes |
| A1 on zero | 2 | 633 | fail | rhat, ess_bulk | yes |
| sg0.25_r1 | 2 | 544 | pass |  | yes |
| sg0.25_r2 | 2 | 633 | fail | rhat, ess_bulk, ess_tail | yes |
| sg0.25_r3 | 2 | 614 | fail | rhat, ess_bulk, ess_tail | yes |
| sg0.25_r4 | 2 | 694 | fail | rhat, ess_bulk | yes |
| sg0.25_r5 | 2 | 636 | fail | rhat, ess_bulk | yes |
| sg0.5_r1 | 2 | 684 | fail | divergences | yes |
| sg0.5_r2 | 2 | 649 | fail | rhat, ess_bulk, ess_tail | yes |
| sg0.5_r3 | 2 | 625 | fail | rhat | yes |
| sg0.5_r4 | 2 | 656 | fail | rhat, ess_bulk, ess_tail | yes |
| sg0.5_r5 | 2 | 670 | pass |  | yes |
| sg1_r1 | 2 | 661 | fail | rhat, ess_bulk | no |
| sg1_r2 | 2 | 1106 | fail | rhat, ess_bulk, ess_tail | no |
| sg2_r1 | 2 | 650 | fail | rhat, ess_bulk, ess_tail | no |
| sg2_r2 | 2 | 506 | fail | rhat, ess_bulk | no |
