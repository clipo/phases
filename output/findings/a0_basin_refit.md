# A0 refit on the real basin data, with corrected diagnostics and the rule-20 program

Produced by `analyses/probe_latent_groups/refit_a0_basin.R` on 2026-10-16 (pilot step 4 of the latent-group probe; exploratory design validation, not evidence about the record).

## Data contract and provenance

Prior in the contract: `prior_rho_type = 2` (1 inverse-gamma, 2 proportional to 1/rho, 3 uniform), rho on [0.5, 300] km, half-Normal(0, 2.0) on sigma and tau, kernel 1. This settles the lineage question of inventory row 14: the contract Carl's basin fit read carries the 1/rho prior; only the two repair scripts switched to the inverse-gamma.
Software: R version 4.5.3 (2026-03-11); rstan 2.32.7; StanHeaders 2.39.1; posterior 1.7.0; loo 2.10.1; bayesplot 1.15.0. Stan source sha256 `83410f2962dec4d4ed0fad6dc255465fbea9568b7c9bedab683419ebb3b4364b`; fit layer sha256 `d0ccf69224fc13266b7484dc8ba60d1fb5bec829839b9da271f0be85eb08f00a`.

## Changes to the model file since bb2fd21 (model block unchanged except the third rho-prior case)

```diff
diff --git a/stan/spatial/composition_gp_marginal.stan b/stan/spatial/composition_gp_marginal.stan
index c4107ff..0c63d7c 100644
--- a/stan/spatial/composition_gp_marginal.stan
+++ b/stan/spatial/composition_gp_marginal.stan
@@ -64,10 +64,11 @@ data {
   array[N, K] int<lower=0> y;           // counts
   matrix[N, N] D;                       // pairwise distance, km
 
-  // rho prior. 1 = inv_gamma(a, b); 2 = uniform on log(rho) in [rho_lo, rho_hi].
-  // Both are offered because the choice is contested and rule 20 requires the
-  // sensitivity refit to lean the other way; neither is hard-coded.
-  int<lower=1, upper=2> prior_rho_type;
+  // rho prior. 1 = inv_gamma(a, b); 2 = uniform on log(rho) in [rho_lo, rho_hi];
+  // 3 = uniform on rho in [rho_lo, rho_hi] (added 2026-10-16 as the rule-20
+  // sensitivity prior leaning toward long ranges; latent-group probe spec 4.1).
+  // Offered as data because the choice is contested; none is hard-coded.
+  int<lower=1, upper=3> prior_rho_type;
   real<lower=0> prior_rho_a;
   real<lower=0> prior_rho_b;
   real<lower=0> rho_lo;
@@ -110,7 +111,8 @@ model {
   if (use_tau) tau ~ normal(0, prior_scale_sd);
   sigma ~ normal(0, prior_scale_sd);
   if (prior_rho_type == 1) rho ~ inv_gamma(prior_rho_a, prior_rho_b);
-  else target += -log(rho);             // p(rho) proportional to 1/rho
+  else if (prior_rho_type == 2) target += -log(rho);   // p(rho) proportional to 1/rho
+  // prior_rho_type == 3: uniform on [rho_lo, rho_hi], nothing to add.
   alpha ~ normal(0, 2);
   to_vector(z_f) ~ std_normal();
 
@@ -128,12 +130,14 @@ generated quantities {
   real half_distance = hc * rho;        // correlation 0.5; the primary scale
   real spatial_share = use_tau ? square(sigma) / (square(sigma) + square(tau[1])) : 1;
   array[N, K] int y_rep;
+  vector[N] log_lik;                    // pointwise, for held-out scoring (added 2026-10-16, probe spec 4)
   {
     for (i in 1:N) {
       vector[K] eta;
       eta[1] = 0;
       for (k in 2:K) eta[k] = alpha[k - 1] + f[i, k - 1];
       y_rep[i] = multinomial_logit_rng(eta, n_i[i]);
+      log_lik[i] = multinomial_logit_lpmf(y[i] | eta);
     }
   }
 }
```

## 1. The refit

Gate: PASS. Rung 1 (iter 3000, warmup 1500, adapt_delta 0.95, max_treedepth 12); 87 s wall time (4 chains in parallel); read from cache.
- a0_basin: R-hat max 1.0034 (lp__) | bulk ESS min 608 (tau[1]) | tail ESS min 1249 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.676 | stuck 0
Treedepth, measured against the configured maximum for the first time: 0 of 6000 post-warmup iterations (0.0%) at depth 12; largest depth observed 10. Gated variables: 293 (sigma, tau, rho, 9 alpha, lp__, 252 field coordinates, 28 pointwise log-likelihoods); exempt 0; stuck 0.

| quantity | median | 5% | 95% |
|---|---|---|---|
| sigma | 4.134 | 3.165 | 5.286 |
| tau[1] | 0.466 | 0.089 | 0.734 |
| rho | 204.142 | 109.260 | 287.255 |
| spatial_share | 0.988 | 0.961 | 1.000 |
| half_distance | 141.501 | 75.733 | 199.110 |
| practical_range | 612.427 | 327.780 | 861.764 |

Visual diagnostics (Liberation Sans, 300 dpi, 7 in): `a0_basin_rank.png`, `a0_basin_trace.png`, `a0_basin_pairs.png`, `a0_basin_energy.png`; plus `a0_basin_pairs_alpha.png` (pairs on the nine intercepts).
Posterior-predictive counts for 200 draws: `output/probe_latent_groups/a0_posterior_predictive_yrep.csv` (summarized by the Python script).

## 2a. Prior predictive (rule 20a)

Sampled with the model's own `likelihood = 0` switch. - a0_prior_predictive: R-hat max 1.0042 (alpha[4]) | bulk ESS min 1740 (lp__) | tail ESS min 1783 (sigma) | divergences 0/4000 (0.00%) | treedepth at max 10: 0/4000 (0.0%) | E-BFMI min 0.958 | stuck 0
Counts for 200 draws: `output/probe_latent_groups/a0_prior_predictive_yrep.csv`. The resemblance summaries (mean within-assemblage Gini-Simpson, mean richness, pooled Gini-Simpson, two-cluster plug-in F_ST) are computed beside the observed values by `analyses/probe_latent_groups/a0_predictive_summaries.py`, which owns those definitions (spec 3.3); its output is `output/findings/a0_predictive_summaries.md`.

## 2b. Prior sensitivity (rule 20b)

Direction of each default and what the alternative leans toward: half-Normal(0, 2) on sigma and tau admits large variance of both kinds; half-Normal(0, 0.5) leans toward little of either. The 1/rho prior puts half its mass below 12.2 km (the geometric mean of 0.5 and 300); the uniform prior on [0.5, 300] puts 90 percent of its mass above 30 km and leans toward long ranges.

| quantity | default median [5%, 95%] | half-Normal(0, 0.5) on sigma, tau | uniform rho |
|---|---|---|---|
| sigma | 4.134 [3.165, 5.286] | 2.148 [1.701, 2.663] | 4.228 [3.256, 5.298] |
| tau[1] | 0.466 [0.089, 0.734] | 0.498 [0.118, 0.741] | 0.479 [0.112, 0.739] |
| rho | 204.142 [109.260, 287.255] | 59.944 [30.134, 114.050] | 219.655 [124.035, 292.191] |
| spatial_share | 0.988 [0.961, 1.000] | 0.951 [0.878, 0.997] | 0.987 [0.961, 0.999] |
| half_distance | 141.501 [75.733, 199.110] | 41.550 [20.887, 79.053] | 152.254 [85.975, 202.531] |

- a0_sens_scale05: R-hat max 1.0024 (lp__) | bulk ESS min 919 (tau[1]) | tail ESS min 1708 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.696 | stuck 0
- a0_sens_rho_uniform: R-hat max 1.0045 (tau[1]) | bulk ESS min 1215 (tau[1]) | tail ESS min 2236 (tau[1]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.660 | stuck 0

### What is identified (D20)

The data see between-site contrasts, with variance 2 sigma^2 (1 - exp(-d/rho)). On 28 sites within 79 km the contrast scale sigma^2/rho and the contrast SD at distances inside the basin are informed; sigma, rho and the intercepts individually are not (the intercepts are anti-correlated with the field's per-class mean, correlations -0.96 to -1.00 on the default fit). Distances below are the nearest-neighbor median, the median pair and the maximum pair. Length scales in km and sigma as an amplitude are not reportable as estimates; these are.

| quantity | default prior | half-Normal(0, 0.5) on sigma, tau | uniform rho |
|---|---|---|---|
| sigma^2/rho (logit^2 per km) | 0.09 [0.06, 0.13] | 0.08 [0.05, 0.13] | 0.08 [0.06, 0.13] |
| field contrast SD at 4.3 km (logits) | 0.87 [0.69, 1.05] | 0.80 [0.62, 1.01] | 0.85 [0.69, 1.04] |
| contrast SD with nugget at 4.3 km | 1.10 [0.93, 1.31] | 1.07 [0.91, 1.28] | 1.10 [0.93, 1.32] |
| field contrast SD at 36.6 km (logits) | 2.42 [1.95, 2.92] | 2.04 [1.65, 2.45] | 2.39 [1.94, 2.89] |
| contrast SD with nugget at 36.6 km | 2.52 [2.14, 2.95] | 2.17 [1.88, 2.50] | 2.49 [2.13, 2.93] |
| field contrast SD at 79.1 km (logits) | 3.37 [2.74, 4.04] | 2.56 [2.14, 3.03] | 3.34 [2.73, 4.03] |
| contrast SD with nugget at 79.1 km | 3.44 [2.87, 4.07] | 2.66 [2.30, 3.08] | 3.41 [2.87, 4.05] |

## 2c. Recovery against the prior (rule 20c)

Data simulated from the model by `analyses/49_write_stan_data.simulate_from_model` (numpy, no Stan) at the real geometry and sherd totals; spatial share 0.2 and 0.9 crossed with rho 5 and 60 km; total latent variance 1 logit^2. Covered = truth inside the 5 to 95 percent interval.

| record | quantity | truth | median | 5% | 95% | covered |
|---|---|---|---|---|---|---|
| share0.2_rho5 | spatial_share | 0.200 | 0.361 | 0.004 | 0.990 | yes |
| share0.2_rho5 | rho | 5.000 | 1.896 | 0.599 | 131.846 | yes |
| share0.2_rho5 | sigma | 0.447 | 0.614 | 0.060 | 1.023 | yes |
| share0.2_rho5 | tau[1] | 0.894 | 0.811 | 0.101 | 1.046 | yes |
| share0.2_rho60 | spatial_share | 0.200 | 0.155 | 0.001 | 0.961 | yes |
| share0.2_rho60 | rho | 60.000 | 2.515 | 0.574 | 189.618 | yes |
| share0.2_rho60 | sigma | 0.447 | 0.428 | 0.038 | 1.060 | yes |
| share0.2_rho60 | tau[1] | 0.894 | 0.983 | 0.210 | 1.133 | yes |
| share0.9_rho5 | spatial_share | 0.900 | 0.923 | 0.529 | 0.999 | yes |
| share0.9_rho5 | rho | 5.000 | 4.242 | 2.450 | 8.361 | yes |
| share0.9_rho5 | sigma | 0.949 | 0.974 | 0.741 | 1.106 | yes |
| share0.9_rho5 | tau[1] | 0.316 | 0.285 | 0.028 | 0.715 | yes |
| share0.9_rho60 | spatial_share | 0.900 | 0.900 | 0.787 | 0.963 | yes |
| share0.9_rho60 | rho | 60.000 | 67.208 | 23.048 | 224.785 | yes |
| share0.9_rho60 | sigma | 0.949 | 0.945 | 0.658 | 1.581 | yes |
| share0.9_rho60 | tau[1] | 0.316 | 0.326 | 0.230 | 0.412 | yes |

- share0.2_rho5: R-hat max 1.0076 (alpha[2]) | bulk ESS min 1142 (alpha[2]) | tail ESS min 1919 (rho) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.657 | stuck 0
- share0.2_rho60: R-hat max 1.0039 (rho) | bulk ESS min 1125 (alpha[2]) | tail ESS min 1684 (alpha[6]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.656 | stuck 0
- share0.9_rho5: R-hat max 1.0041 (alpha[1]) | bulk ESS min 1233 (alpha[2]) | tail ESS min 2041 (rho) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.617 | stuck 0
- share0.9_rho60: R-hat max 1.0036 (f[8,3]) | bulk ESS min 873 (rho) | tail ESS min 1430 (f[28,2]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.650 | stuck 0

## Timing and gates

| fit | rung | seconds | gate | reasons | cached |
|---|---|---|---|---|---|
| a0_basin | 1 | 87 | pass |  | yes |
| a0_prior_predictive | 1 | 4 | pass |  | yes |
| a0_sens_scale05 | 1 | 53 | pass |  | yes |
| a0_sens_rho_uniform | 2 | 263 | pass |  | yes |
| a0_recovery_share0.2_rho5 | 2 | 179 | pass |  | yes |
| a0_recovery_share0.2_rho60 | 2 | 181 | pass |  | yes |
| a0_recovery_share0.9_rho5 | 2 | 94 | pass |  | yes |
| a0_recovery_share0.9_rho60 | 1 | 44 | pass |  | yes |

Fit records (JSON, with every diagnostic, hash and version) are beside the cached fits in `output/gp_cache/`.
