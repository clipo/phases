# A0 refit on the real basin data, with corrected diagnostics and the rule-20 program

Produced by `analyses/probe_latent_groups/refit_a0_basin.R` on 2026-11-03 (pilot step 4 of the latent-group probe; exploratory design validation, not evidence about the record).

## Data contract and provenance

Prior in the contract: `prior_rho_type = 2` (1 inverse-gamma, 2 proportional to 1/rho, 3 uniform), rho on [0.5, 300] km, half-Normal(0, 2.0) on sigma and tau, Normal(0, 10) on the class intercepts (`prior_alpha_mean`, `prior_alpha_sd`; D25 made the intercept prior data on 2026-10-18, default Normal(0, 10); the 2026-10-16 run of this script used the hard-coded Normal(0, 2) and is kept as `a0_basin_refit_alphaN0-2.md`), kernel 1. This settles the lineage question of inventory row 14: the contract Carl's basin fit read carries the 1/rho prior; only the two repair scripts switched to the inverse-gamma.
Software: R version 4.6.1 (2026-06-24); rstan 2.32.7; StanHeaders 2.39.1; posterior 1.7.0; loo 2.10.1; bayesplot 1.16.0. Stan source sha256 `35fbf97122c9ab244be792a09756b4997094defe3fa1ea181890c703dae5b902`; fit layer sha256 `79d141d356a09fd4cef38c97343b0485791af1ffcf18c65d489543f84a6ab1a1`.

## Changes to the model file since bb2fd21 (model block unchanged except the third rho-prior case)

```diff
diff --git a/stan/spatial/composition_gp_marginal.stan b/stan/spatial/composition_gp_marginal.stan
index c4107ff..34f6b01 100644
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
@@ -76,6 +77,11 @@ data {
   int<lower=1, upper=3> kernel;         // 1 exp, 2 Matern 3/2, 3 Matern 5/2
   int<lower=0, upper=1> use_tau;        // 0 removes the non-spatial term
   real<lower=0> prior_scale_sd;         // half-normal sd, tau and sigma alike
+  // Intercept prior as data (2026-10-18, D25): the hard-coded Normal(0, 2) was too narrow for the
+  // log-odds of rare classes against the reference, so the field's constant mode carried 0.3 to 5.6
+  // logits of the class baselines (discussion log, 2026-10-18).
+  vector[K - 1] prior_alpha_mean;
+  real<lower=0> prior_alpha_sd;
   int<lower=0, upper=1> likelihood;     // 0 = prior predictive
 }
 transformed data {
@@ -110,8 +116,9 @@ model {
   if (use_tau) tau ~ normal(0, prior_scale_sd);
   sigma ~ normal(0, prior_scale_sd);
   if (prior_rho_type == 1) rho ~ inv_gamma(prior_rho_a, prior_rho_b);
-  else target += -log(rho);             // p(rho) proportional to 1/rho
-  alpha ~ normal(0, 2);
+  else if (prior_rho_type == 2) target += -log(rho);   // p(rho) proportional to 1/rho
+  // prior_rho_type == 3: uniform on [rho_lo, rho_hi], nothing to add.
+  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
   to_vector(z_f) ~ std_normal();
 
   if (likelihood) {
@@ -128,12 +135,14 @@ generated quantities {
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

Gate: PASS. Rung 1 (iter 3000, warmup 1500, adapt_delta 0.95, max_treedepth 12); 125 s wall time (4 chains in parallel); read from cache.
- a0_basin: R-hat max 1.0037 (tau[1]) | bulk ESS min 628 (tau[1]) | tail ESS min 1535 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.675 | stuck 0
Treedepth, measured against the configured maximum for the first time: 0 of 6000 post-warmup iterations (0.0%) at depth 12; largest depth observed 10. Gated variables: 293 (sigma, tau, rho, 9 alpha, lp__, 252 field coordinates, 28 pointwise log-likelihoods); exempt 0; stuck 0.

| quantity | median | 5% | 95% |
|---|---|---|---|
| sigma | 3.022 | 2.070 | 4.563 |
| tau[1] | 0.454 | 0.081 | 0.742 |
| rho | 97.833 | 37.681 | 247.010 |
| spatial_share | 0.980 | 0.924 | 0.999 |
| half_distance | 67.813 | 26.118 | 171.214 |
| practical_range | 293.499 | 113.043 | 741.030 |

Visual diagnostics (Arial, 300 dpi, 7 in): `a0_basin_rank.png`, `a0_basin_trace.png`, `a0_basin_pairs.png`, `a0_basin_energy.png`; plus `a0_basin_pairs_alpha.png` (pairs on the nine intercepts).
Posterior-predictive counts for 200 draws: `output/probe_latent_groups/a0_posterior_predictive_yrep.csv` (summarized by the Python script).

## 2a. Prior predictive (rule 20a)

Sampled with the model's own `likelihood = 0` switch. - a0_prior_predictive: R-hat max 1.0056 (alpha[1]) | bulk ESS min 1469 (lp__) | tail ESS min 1831 (tau[1]) | divergences 0/4000 (0.00%) | treedepth at max 10: 0/4000 (0.0%) | E-BFMI min 0.963 | stuck 0
Counts for 200 draws: `output/probe_latent_groups/a0_prior_predictive_yrep.csv`. The resemblance summaries (mean within-assemblage Gini-Simpson, mean richness, pooled Gini-Simpson, two-cluster plug-in F_ST) are computed beside the observed values by `analyses/probe_latent_groups/a0_predictive_summaries.py`, which owns those definitions (spec 3.3); its output is `output/findings/a0_predictive_summaries.md`.

## 2b. Prior sensitivity (rule 20b)

Direction of each default and what the alternative leans toward: half-Normal(0, 2) on sigma and tau admits large variance of both kinds; half-Normal(0, 0.5) leans toward little of either. The 1/rho prior puts half its mass below 12.2 km (the geometric mean of 0.5 and 300); the uniform prior on [0.5, 300] puts 90 percent of its mass above 30 km and leans toward long ranges.

| quantity | default median [5%, 95%] | half-Normal(0, 0.5) on sigma, tau | uniform rho |
|---|---|---|---|
| sigma | 3.022 [2.070, 4.563] | 1.807 [1.505, 2.203] | 3.425 [2.288, 4.859] |
| tau[1] | 0.454 [0.081, 0.742] | 0.398 [0.057, 0.694] | 0.477 [0.104, 0.750] |
| rho | 97.833 [37.681, 247.010] | 32.108 [18.962, 58.391] | 138.834 [50.835, 270.121] |
| spatial_share | 0.980 [0.924, 0.999] | 0.955 [0.855, 0.999] | 0.982 [0.937, 0.999] |
| half_distance | 67.813 [26.118, 171.214] | 22.256 [13.143, 40.474] | 96.233 [35.236, 187.234] |

- a0_sens_scale05: R-hat max 1.0043 (rho) | bulk ESS min 905 (f[3,1]) | tail ESS min 1748 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.715 | stuck 0
- a0_sens_rho_uniform: R-hat max 1.0028 (f[4,1]) | bulk ESS min 1349 (tau[1]) | tail ESS min 2667 (tau[1]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.643 | stuck 0

### What is identified (D20)

The data see between-site contrasts, with variance 2 sigma^2 (1 - exp(-d/rho)). On 28 sites within 79 km the contrast scale sigma^2/rho and the contrast SD at distances inside the basin are informed; sigma, rho and the intercepts individually are not (the intercepts are anti-correlated with the field's per-class mean, correlations -0.96 to -1.00 on the default fit). Distances below are the nearest-neighbor median, the median pair and the maximum pair. Length scales in km and sigma as an amplitude are not reportable as estimates; these are.

| quantity | default prior | half-Normal(0, 0.5) on sigma, tau | uniform rho |
|---|---|---|---|
| sigma^2/rho (logit^2 per km) | 0.10 [0.06, 0.15] | 0.10 [0.06, 0.16] | 0.09 [0.05, 0.14] |
| field contrast SD at 4.3 km (logits) | 0.89 [0.70, 1.10] | 0.91 [0.71, 1.12] | 0.87 [0.68, 1.07] |
| contrast SD with nugget at 4.3 km | 1.12 [0.95, 1.33] | 1.09 [0.93, 1.29] | 1.11 [0.94, 1.33] |
| field contrast SD at 36.6 km (logits) | 2.38 [1.92, 2.85] | 2.10 [1.74, 2.44] | 2.36 [1.88, 2.86] |
| contrast SD with nugget at 36.6 km | 2.47 [2.11, 2.89] | 2.18 [1.92, 2.49] | 2.46 [2.09, 2.89] |
| field contrast SD at 79.1 km (logits) | 3.13 [2.52, 3.81] | 2.43 [2.05, 2.85] | 3.18 [2.55, 3.88] |
| contrast SD with nugget at 79.1 km | 3.20 [2.66, 3.84] | 2.50 [2.18, 2.90] | 3.26 [2.71, 3.91] |

## 2c. Recovery against the prior (rule 20c)

Data simulated from the model by `analyses/49_write_stan_data.simulate_from_model` (numpy, no Stan) at the real geometry and sherd totals; spatial share 0.2 and 0.9 crossed with rho 5 and 60 km; total latent variance 1 logit^2. Covered = truth inside the 5 to 95 percent interval.

| record | quantity | truth | median | 5% | 95% | covered |
|---|---|---|---|---|---|---|
| share0.2_rho5 | spatial_share | 0.200 | 0.348 | 0.005 | 0.987 | yes |
| share0.2_rho5 | rho | 5.000 | 2.005 | 0.601 | 154.042 | yes |
| share0.2_rho5 | sigma | 0.447 | 0.606 | 0.068 | 1.024 | yes |
| share0.2_rho5 | tau[1] | 0.894 | 0.817 | 0.113 | 1.053 | yes |
| share0.2_rho60 | spatial_share | 0.200 | 0.170 | 0.001 | 0.966 | yes |
| share0.2_rho60 | rho | 60.000 | 2.802 | 0.586 | 214.036 | yes |
| share0.2_rho60 | sigma | 0.447 | 0.452 | 0.041 | 1.058 | yes |
| share0.2_rho60 | tau[1] | 0.894 | 0.988 | 0.197 | 1.138 | yes |
| share0.9_rho5 | spatial_share | 0.900 | 0.917 | 0.458 | 0.999 | yes |
| share0.9_rho5 | rho | 5.000 | 4.314 | 2.443 | 11.623 | yes |
| share0.9_rho5 | sigma | 0.949 | 0.973 | 0.702 | 1.111 | yes |
| share0.9_rho5 | tau[1] | 0.316 | 0.299 | 0.029 | 0.803 | yes |
| share0.9_rho60 | spatial_share | 0.900 | 0.917 | 0.792 | 0.974 | yes |
| share0.9_rho60 | rho | 60.000 | 88.777 | 25.126 | 260.872 | yes |
| share0.9_rho60 | sigma | 0.949 | 1.069 | 0.675 | 1.794 | yes |
| share0.9_rho60 | tau[1] | 0.316 | 0.330 | 0.234 | 0.415 | yes |

- share0.2_rho5: R-hat max 1.0037 (f[3,3]) | bulk ESS min 1127 (alpha[2]) | tail ESS min 1467 (alpha[2]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.667 | stuck 0
- share0.2_rho60: R-hat max 1.0062 (alpha[3]) | bulk ESS min 905 (rho) | tail ESS min 1046 (alpha[2]) | divergences 0/12000 (0.00%) | treedepth at max 13: 0/12000 (0.0%) | E-BFMI min 0.665 | stuck 0
- share0.9_rho5: R-hat max 1.0028 (alpha[3]) | bulk ESS min 1323 (rho) | tail ESS min 596 (rho) | divergences 0/24000 (0.00%) | treedepth at max 14: 0/24000 (0.0%) | E-BFMI min 0.664 | stuck 0
- share0.9_rho60: R-hat max 1.0024 (f[5,6]) | bulk ESS min 1075 (rho) | tail ESS min 1410 (tau[1]) | divergences 0/6000 (0.00%) | treedepth at max 12: 0/6000 (0.0%) | E-BFMI min 0.705 | stuck 0

## Timing and gates

| fit | rung | seconds | gate | reasons | cached |
|---|---|---|---|---|---|
| a0_basin | 1 | 125 | pass |  | yes |
| a0_prior_predictive | 1 | 4 | pass |  | yes |
| a0_sens_scale05 | 1 | 52 | pass |  | yes |
| a0_sens_rho_uniform | 2 | 702 | pass |  | yes |
| a0_recovery_share0.2_rho5 | 2 | 247 | pass |  | no |
| a0_recovery_share0.2_rho60 | 2 | 222 | pass |  | no |
| a0_recovery_share0.9_rho5 | 3 | 657 | pass |  | no |
| a0_recovery_share0.9_rho60 | 1 | 66 | pass |  | no |

Fit records (JSON, with every diagnostic, hash and version) are beside the cached fits in `output/gp_cache/`.
