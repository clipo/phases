# A1 validated on its own generator (pilot step 6; spec 4.1c)

Produced by `analyses/probe_latent_groups/validate_a1.R --alpha-sd=2` on 2026-11-23, fitted under the Normal(0, 2) class-intercept prior. Records from `a1_contracts.py --alpha-sd 2` (nuisances pinned to the A0 refit's identified quantities; field centered per class; membership = Phillips's partition 12/11/5; Gabriel edges for coherence). Exploratory design validation, not evidence about the record.

## (i) Code equivalence at zero group scale

A0 and A1 with the group scale fixed at zero, same record, same seed. Largest |z| over 41 shared posterior means (difference over combined MCSE): **2.69** at `log_lik[28]`. Values above 3 would mean the two programs disagree.

| quantity | A0 mean | A1 (fixed 0) mean | z |
|---|---|---|---|
| sigma | 4.1671 | 4.1760 | 0.43 |
| tau[1] | 0.6416 | 0.6409 | -0.11 |
| rho | 192.0709 | 190.1850 | -1.07 |
| alpha[1] | -0.2075 | -0.1964 | 0.22 |
| alpha[2] | -0.3382 | -0.3401 | -0.04 |
| alpha[3] | -1.2516 | -1.2446 | 0.18 |
| alpha[4] | -0.5810 | -0.5695 | 0.27 |
| alpha[5] | -1.4595 | -1.4084 | 1.24 |
| alpha[6] | -1.6662 | -1.6679 | -0.04 |
| alpha[7] | -1.4945 | -1.4712 | 0.58 |
| alpha[8] | -1.2514 | -1.2828 | -0.79 |
| alpha[9] | -1.6653 | -1.6566 | 0.23 |
| spatial_share | 0.9736 | 0.9739 | 0.45 |
| log_lik[1..28] | | | max abs z 2.69 |

- A0 on zero: R-hat max 1.0051 (tau[1]) | bulk ESS min 947 (tau[1]) | tail ESS min 1440 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.779 | stuck 0
- A1 fixed 0 on zero: R-hat max 1.0024 (f[5,4]) | bulk ESS min 1165 (tau[1]) | tail ESS min 1806 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.759 | stuck 0

## (ii) Null calibration: unconstrained A1 on the zero-shift record

Truth: sigma_g = 0, no groups; the generating labels (Phillips's partition) are used only to read the margins. Realized between-'group' F_ST of the record from the field alone: 0.3076 (the record's is 0.0123; see the A0 adequacy note below).

| quantity | median [5%, 95%] |
|---|---|
| sigma_g | 0.442 [0.044, 1.201] |
| Delta | 0.850 [0.082, 2.027] |
| Delta_w | 0.807 [0.062, 1.909] |
| coherence (soft) | 0.003 [-0.062, 0.095] |
| coherence (hard) | 0.006 [-0.092, 0.119] |
| separation margins vs generating labels (12, 11, 5) | 0.008, -0.010, -0.010 |
| per-chain co-membership split (max abs diff) | 0.244 |

- A1 on zero: R-hat max 1.0457 (comembership[25,21]) | bulk ESS min 120 (comembership[25,21]) | tail ESS min 262 (comembership[9,8]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.663 | stuck 0
No figures are written for this arm.

## (iii) Recovery against the prior (sigma_g = 0.25 and 0.5, five records each)

| record | true sigma_g | realized between-group F_ST | Delta true | sigma_g post | Delta post | Delta_w post | coherence | margins (12, 11, 5) | gate | chain split |
|---|---|---|---|---|---|---|---|---|---|---|
| sg0.25_r1 | 0.25 | 0.2204 | 0.714 | 0.398 [0.044, 1.004] | 0.778 [0.082, 1.729] | 0.737 [0.059, 1.614] | -0.000 [-0.070, 0.069] | -0.01, 0.00, 0.00 | pass | 0.098 |
| sg0.25_r2 | 0.25 | 0.1345 | 0.540 | 0.372 [0.032, 4.674] | 0.696 [0.062, 5.305] | 0.407 [0.000, 1.428] | 0.000 [-0.063, 0.072] | -0.00, -0.00, 0.01 | fail: rhat,ess_bulk,ess_tail,divergences | 0.168 |
| sg0.25_r3 | 0.25 | 0.1329 | 0.742 | 0.605 [0.069, 1.526] | 1.023 [0.137, 2.054] | 0.907 [0.099, 1.700] | -0.000 [-0.068, 0.085] | 0.04, -0.04, -0.03 | fail: rhat,ess_bulk | 0.119 |
| sg0.25_r4 | 0.25 | 0.0187 | 0.561 | 0.278 [0.029, 1.012] | 0.534 [0.055, 1.431] | 0.479 [0.023, 1.175] | 0.000 [-0.066, 0.082] | -0.00, 0.00, 0.00 | fail: rhat | 0.105 |
| sg0.25_r5 | 0.25 | 0.0805 | 0.544 | 0.364 [0.040, 0.987] | 0.707 [0.076, 1.559] | 0.646 [0.035, 1.382] | -0.001 [-0.081, 0.075] | 0.00, -0.00, -0.03 | fail: rhat,ess_bulk | 0.186 |
| sg0.5_r1 | 0.50 | 0.0429 | 1.298 | 0.269 [0.028, 0.937] | 0.518 [0.052, 1.424] | 0.464 [0.023, 1.199] | 0.000 [-0.068, 0.079] | 0.00, -0.00, -0.00 | fail: rhat | 0.100 |
| sg0.5_r2 | 0.50 | 0.0967 | 1.092 | 0.338 [0.034, 1.824] | 0.639 [0.064, 2.128] | 0.568 [0.023, 1.459] | -0.000 [-0.073, 0.076] | -0.01, 0.01, -0.00 | fail: rhat,divergences | 0.127 |
| sg0.5_r3 | 0.50 | 0.5750 | 1.583 | 0.479 [0.052, 1.460] | 0.937 [0.099, 2.275] | 0.880 [0.083, 2.034] | 0.000 [-0.065, 0.074] | 0.01, -0.01, 0.00 | fail: rhat,ess_bulk | 0.092 |
| sg0.5_r4 | 0.50 | 0.0300 | 1.024 | 0.351 [0.034, 1.227] | 0.669 [0.067, 1.838] | 0.593 [0.019, 1.554] | -0.000 [-0.085, 0.079] | 0.00, -0.00, -0.01 | fail: rhat | 0.099 |
| sg0.5_r5 | 0.50 | 0.0474 | 1.153 | 0.290 [0.027, 1.357] | 0.551 [0.052, 1.665] | 0.478 [0.006, 1.262] | 0.000 [-0.071, 0.079] | 0.00, 0.00, -0.00 | fail: rhat | 0.131 |
| sg1_r1 | 1.00 | 0.0790 | 2.338 | 0.494 [0.058, 1.206] | 0.956 [0.113, 1.813] | 0.899 [0.070, 1.616] | -0.001 [-0.090, 0.084] | -0.03, 0.03, -0.00 | fail: rhat,ess_bulk | 0.250 |
| sg1_r2 | 1.00 | 0.0459 | 1.636 | 0.423 [0.040, 3.115] | 0.801 [0.076, 3.027] | 0.636 [0.000, 1.741] | -0.000 [-0.089, 0.070] | 0.00, 0.00, -0.01 | fail: rhat,ess_bulk,ess_tail | 0.204 |
| sg2_r1 | 2.00 | 0.3560 | 5.401 | 1.344 [0.141, 2.807] | 1.998 [0.278, 3.916] | 1.679 [0.260, 3.260] | 0.001 [-0.048, 0.120] | 0.05, -0.03, 0.07 | fail: rhat,ess_bulk,ess_tail | 0.295 |
| sg2_r2 | 2.00 | 0.2339 | 3.017 | 0.485 [0.041, 1.570] | 0.925 [0.079, 2.467] | 0.854 [0.060, 2.187] | -0.003 [-0.093, 0.067] | 0.00, 0.00, -0.01 | fail: rhat,ess_bulk,ess_tail | 0.131 |

## A0 adequacy note (found while pinning the generator, 2026-10-17)

The A0 posterior predictive with the fitted field reproduces the record's between-area F_ST (0.0122 against 0.0123), but fresh field draws at the same hyperparameters give a median of 0.085 (10th to 90th percentile 0.02 to 0.26). Between Phillips's Parkin and Kent areas the pooled logit share differs by +0.47 for Parkin_Punctated, -0.12 for Barton/Kent/MPI, -0.90 for Painted, and from -2.2 to +0.83 [corrected 2026-10-19 after the Codex audit of the regional-drift spec: Fortune_Noded is +0.83 in area_logits.md; first written as -0.17 to -2.2, which omitted it] for the rare classes, two of which (Walls_Engraved, Hull_Engraved) are absent from every Parkin-area site; per-class SDs of the logit share across sites are 0.8 to 1.8 (`analyses/probe_latent_groups/area_logits.py`, `output/findings/area_logits.md`). The reference class is not the flattest on this scale. A0 has one field scale for all classes and holds the record's pattern only in the realized field. A per-class field scale is the extension this points to. The own-generator records above therefore carry far more between-area contrast than the record (realized F_ST column), which is acceptable for a test of the software and of identification under the model (rule 20c) but is not a resemblance claim.

## Timing and gates

| fit | rung | seconds | gate | reasons | cached |
|---|---|---|---|---|---|
| A0 on zero | 1 | 99 | pass |  | yes |
| A1 fixed 0 on zero | 1 | 189 | pass |  | yes |
| A1 on zero | 2 | 780 | fail | rhat, ess_bulk, ess_tail | yes |
| sg0.25_r1 | 2 | 731 | pass |  | yes |
| sg0.25_r2 | 2 | 1290 | fail | rhat, ess_bulk, ess_tail, divergences | yes |
| sg0.25_r3 | 2 | 872 | fail | rhat, ess_bulk | no |
| sg0.25_r4 | 2 | 842 | fail | rhat | no |
| sg0.25_r5 | 2 | 741 | fail | rhat, ess_bulk | no |
| sg0.5_r1 | 2 | 696 | fail | rhat | no |
| sg0.5_r2 | 2 | 674 | fail | rhat, divergences | no |
| sg0.5_r3 | 2 | 438 | fail | rhat, ess_bulk | no |
| sg0.5_r4 | 2 | 700 | fail | rhat | no |
| sg0.5_r5 | 2 | 674 | fail | rhat | no |
| sg1_r1 | 2 | 582 | fail | rhat, ess_bulk | no |
| sg1_r2 | 2 | 662 | fail | rhat, ess_bulk, ess_tail | no |
| sg2_r1 | 2 | 535 | fail | rhat, ess_bulk, ess_tail | no |
| sg2_r2 | 2 | 556 | fail | rhat, ess_bulk, ess_tail | no |
