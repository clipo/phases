# A1 validated on its own generator (pilot step 6; spec 4.1c)

Produced by `analyses/probe_latent_groups/validate_a1.R` on 2026-11-22. Records from `a1_contracts.py` (nuisances pinned to the A0 refit's identified quantities; field centered per class; membership = Phillips's partition 12/11/5; Gabriel edges for coherence). Exploratory design validation, not evidence about the record.

## (i) Code equivalence at zero group scale

A0 and A1 with the group scale fixed at zero, same record, same seed. Largest |z| over 41 shared posterior means (difference over combined MCSE): **2.18** at `log_lik[6]`. Values above 3 would mean the two programs disagree.

| quantity | A0 mean | A1 (fixed 0) mean | z |
|---|---|---|---|
| sigma | 3.2601 | 3.2490 | -0.49 |
| tau[1] | 0.6455 | 0.6421 | -0.48 |
| rho | 107.1028 | 106.1449 | -0.50 |
| alpha[1] | -0.4832 | -0.5582 | -0.85 |
| alpha[2] | -1.3072 | -1.3461 | -0.56 |
| alpha[3] | -4.8839 | -4.8611 | 0.36 |
| alpha[4] | -2.3183 | -2.2878 | 0.45 |
| alpha[5] | -5.3407 | -5.2664 | 1.27 |
| alpha[6] | -5.8443 | -5.8980 | -0.92 |
| alpha[7] | -6.1446 | -6.1302 | 0.24 |
| alpha[8] | -4.6169 | -4.4963 | 2.06 |
| alpha[9] | -7.3225 | -7.2857 | 0.60 |
| spatial_share | 0.9547 | 0.9551 | 0.39 |
| log_lik[1..28] | | | max abs z 2.18 |

- A0 on zero: R-hat max 1.0029 (tau[1]) | bulk ESS min 1150 (tau[1]) | tail ESS min 1324 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.769 | stuck 0
- A1 fixed 0 on zero: R-hat max 1.0030 (sigma) | bulk ESS min 1446 (tau[1]) | tail ESS min 2020 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.735 | stuck 0

## (ii) Null calibration: unconstrained A1 on the zero-shift record

Truth: sigma_g = 0, no groups; the generating labels (Phillips's partition) are used only to read the margins. Realized between-'group' F_ST of the record from the field alone: 0.3127 (the record's is 0.0123; see the A0 adequacy note below).

| quantity | median [5%, 95%] |
|---|---|
| sigma_g | 0.461 [0.046, 1.297] |
| Delta | 0.881 [0.086, 2.191] |
| Delta_w | 0.841 [0.070, 2.066] |
| coherence (soft) | 0.001 [-0.067, 0.083] |
| coherence (hard) | 0.002 [-0.095, 0.110] |
| separation margins vs generating labels (12, 11, 5) | 0.019, -0.007, -0.007 |
| per-chain co-membership split (max abs diff) | 0.096 |

- A1 on zero: R-hat max 1.0179 (comembership[5,2]) | bulk ESS min 386 (comembership[27,19]) | tail ESS min 821 (comembership[7,5]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.674 | stuck 0
Visual set: `a1_null_rank.png`, `a1_null_trace.png`, `a1_null_pairs.png`, `a1_null_energy.png`.

## (iii) Recovery against the prior (sigma_g = 0.25 and 0.5, five records each)

| record | true sigma_g | realized between-group F_ST | Delta true | sigma_g post | Delta post | Delta_w post | coherence | margins (12, 11, 5) | gate | chain split |
|---|---|---|---|---|---|---|---|---|---|---|
| sg0.25_r1 | 0.25 | 0.2195 | 0.714 | 0.396 [0.039, 1.066] | 0.777 [0.076, 1.824] | 0.741 [0.063, 1.734] | 0.001 [-0.061, 0.081] | -0.01, -0.01, 0.03 | pass | 0.070 |
| sg0.25_r2 | 0.25 | 0.1034 | 0.540 | 0.347 [0.030, 1.272] | 0.646 [0.060, 2.031] | 0.599 [0.046, 1.839] | -0.000 [-0.064, 0.073] | -0.02, -0.02, 0.04 | fail: rhat,ess_bulk,ess_tail,divergences | 0.119 |
| sg0.25_r3 | 0.25 | 0.1279 | 0.742 | 0.411 [0.042, 1.213] | 0.788 [0.080, 1.863] | 0.731 [0.052, 1.670] | 0.001 [-0.064, 0.085] | 0.02, -0.01, -0.01 | fail: rhat | 0.105 |
| sg0.25_r4 | 0.25 | 0.0262 | 0.561 | 0.312 [0.031, 0.892] | 0.610 [0.057, 1.437] | 0.571 [0.044, 1.317] | 0.001 [-0.068, 0.087] | 0.00, -0.00, -0.01 | fail: rhat,ess_bulk | 0.117 |
| sg0.25_r5 | 0.25 | 0.0844 | 0.544 | 0.340 [0.034, 0.981] | 0.653 [0.067, 1.643] | 0.617 [0.052, 1.534] | 0.000 [-0.070, 0.083] | 0.00, -0.00, -0.00 | fail: rhat | 0.116 |
| sg0.5_r1 | 0.50 | 0.0383 | 1.298 | 0.246 [0.025, 0.844] | 0.467 [0.048, 1.304] | 0.429 [0.028, 1.131] | 0.000 [-0.065, 0.077] | 0.00, -0.00, -0.01 | fail: rhat | 0.090 |
| sg0.5_r2 | 0.50 | 0.1035 | 1.092 | 0.373 [0.039, 1.707] | 0.708 [0.074, 2.103] | 0.660 [0.056, 1.698] | -0.002 [-0.078, 0.070] | -0.02, 0.00, 0.00 | fail: rhat,ess_tail | 0.101 |
| sg0.5_r3 | 0.50 | 0.5506 | 1.583 | 0.579 [0.062, 1.596] | 1.139 [0.116, 2.735] | 1.095 [0.104, 2.534] | 0.002 [-0.064, 0.088] | 0.05, -0.02, -0.03 | fail: rhat,ess_bulk,ess_tail | 0.184 |
| sg0.5_r4 | 0.50 | 0.0290 | 1.024 | 0.367 [0.034, 1.059] | 0.697 [0.066, 1.779] | 0.655 [0.049, 1.657] | -0.000 [-0.083, 0.085] | -0.01, 0.00, -0.00 | pass | 0.071 |
| sg0.5_r5 | 0.50 | 0.0395 | 1.153 | 0.282 [0.030, 0.819] | 0.541 [0.058, 1.359] | 0.504 [0.040, 1.245] | 0.000 [-0.069, 0.079] | -0.00, -0.00, 0.01 | fail: rhat | 0.117 |
| sg1_r1 | 1.00 | 0.0837 | 2.338 | 0.569 [0.069, 1.272] | 1.115 [0.132, 2.115] | 1.074 [0.112, 1.952] | -0.001 [-0.087, 0.096] | -0.03, 0.02, -0.01 | fail: rhat,ess_bulk,ess_tail,divergences | 0.374 |
| sg1_r2 | 1.00 | 0.0490 | 1.636 | 0.325 [0.031, 0.988] | 0.618 [0.058, 1.611] | 0.589 [0.048, 1.488] | -0.000 [-0.075, 0.080] | 0.00, -0.00, 0.00 | fail: rhat | 0.137 |
| sg2_r1 | 2.00 | 0.3525 | 5.401 | 1.062 [0.099, 2.917] | 1.948 [0.189, 5.830] | 1.828 [0.183, 5.690] | 0.031 [-0.057, 0.324] | 0.25, -0.01, 0.05 | fail: rhat,ess_bulk,ess_tail,comembership_split | 0.726 |
| sg2_r2 | 2.00 | 0.2178 | 3.017 | 0.568 [0.062, 1.603] | 1.075 [0.119, 2.815] | 1.038 [0.113, 2.651] | -0.005 [-0.090, 0.074] | -0.00, 0.00, -0.00 | fail: rhat,ess_bulk | 0.164 |

## A0 adequacy note (found while pinning the generator, 2026-10-17)

The A0 posterior predictive with the fitted field reproduces the record's between-area F_ST (0.0122 against 0.0123), but fresh field draws at the same hyperparameters give a median of 0.085 (10th to 90th percentile 0.02 to 0.26). Between Phillips's Parkin and Kent areas the pooled logit share differs by +0.47 for Parkin_Punctated, -0.12 for Barton/Kent/MPI, -0.90 for Painted, and from -2.2 to +0.83 [corrected 2026-10-19 after the Codex audit of the regional-drift spec: Fortune_Noded is +0.83 in area_logits.md; first written as -0.17 to -2.2, which omitted it] for the rare classes, two of which (Walls_Engraved, Hull_Engraved) are absent from every Parkin-area site; per-class SDs of the logit share across sites are 0.8 to 1.8 (`analyses/probe_latent_groups/area_logits.py`, `output/findings/area_logits.md`). The reference class is not the flattest on this scale. A0 has one field scale for all classes and holds the record's pattern only in the realized field. A per-class field scale is the extension this points to. The own-generator records above therefore carry far more between-area contrast than the record (realized F_ST column), which is acceptable for a test of the software and of identification under the model (rule 20c) but is not a resemblance claim.

## Timing and gates

| fit | rung | seconds | gate | reasons | cached |
|---|---|---|---|---|---|
| A0 on zero | 1 | 132 | pass |  | no |
| A1 fixed 0 on zero | 1 | 257 | pass |  | no |
| A1 on zero | 2 | 975 | fail | rhat, ess_bulk | no |
| sg0.25_r1 | 2 | 917 | pass |  | no |
| sg0.25_r2 | 2 | 444 | fail | rhat, ess_bulk, ess_tail, divergences | no |
| sg0.25_r3 | 2 | 966 | fail | rhat | no |
| sg0.25_r4 | 2 | 891 | fail | rhat, ess_bulk | no |
| sg0.25_r5 | 2 | 877 | fail | rhat | no |
| sg0.5_r1 | 2 | 976 | fail | rhat | no |
| sg0.5_r2 | 2 | 847 | fail | rhat, ess_tail | no |
| sg0.5_r3 | 2 | 699 | fail | rhat, ess_bulk, ess_tail | no |
| sg0.5_r4 | 2 | 1183 | pass |  | no |
| sg0.5_r5 | 2 | 733 | fail | rhat | no |
| sg1_r1 | 2 | 610 | fail | rhat, ess_bulk, ess_tail, divergences | no |
| sg1_r2 | 2 | 415 | fail | rhat | no |
| sg2_r1 | 2 | 1014 | fail | rhat, ess_bulk, ess_tail, comembership_split | no |
| sg2_r2 | 2 | 690 | fail | rhat, ess_bulk | no |
