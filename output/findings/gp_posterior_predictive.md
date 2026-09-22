# Posterior predictive check: can the model reproduce the basin?

Produced by `analyses/02_spatial/04_posterior_predictive.R` on 2026-09-22.

Statistics chosen to discriminate: quantities the model was not fitted to directly and that the paper's claims rest on.

| model | statistic | observed | predictive median | 95% predictive interval | Bayesian p | verdict |
|---|---|---|---|---|---|---|
| composition_gp | fst | 0.0156 | 0.0155 | [0.0128, 0.0192] | 0.470 | pass |
| composition_gp | decay | -0.3397 | -0.3189 | [-0.3685, -0.2689] | 0.825 | pass |
| composition_gp | div_sd | 0.1123 | 0.1090 | [0.0977, 0.1231] | 0.310 | pass |
| composition_gp | div_mean | 0.5633 | 0.5641 | [0.5510, 0.5755] | 0.552 | pass |
| composition_gp_marginal | fst | 0.0156 | 0.0157 | [0.0130, 0.0191] | 0.517 | pass |
| composition_gp_marginal | decay | -0.3397 | -0.3210 | [-0.3605, -0.2712] | 0.818 | pass |
| composition_gp_marginal | div_sd | 0.1123 | 0.1089 | [0.0957, 0.1209] | 0.315 | pass |
| composition_gp_marginal | div_mean | 0.5633 | 0.5649 | [0.5525, 0.5763] | 0.580 | pass |

## Verdict

**The model reproduces every discriminating statistic.** Observed cultural F_ST, distance decay and diversity spread all fall inside the posterior predictive interval. Misspecification is therefore NOT the explanation for the treedepth saturation on real data, which leaves the flat-likelihood reading: rho is not identified above the design's resolving ceiling, and the posterior piling up against the top of its range is the data saying there is no decay within the basin. `spatial_share` is reportable; a length-scale NUMBER is not.

