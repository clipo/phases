# Posterior predictive check: can the model reproduce the basin?

Produced by `analyses/02_spatial/04_posterior_predictive.R` on 2026-09-22.

Statistics chosen to discriminate: quantities the model was not fitted to directly and that the paper's claims rest on.

| model | statistic | observed | predictive median | 95% predictive interval | Bayesian p | verdict |
|---|---|---|---|---|---|---|
| composition_gp | fst | 0.0156 | 0.0157 | [0.0128, 0.0191] | 0.517 | pass |
| composition_gp | decay | -0.3407 | -0.3185 | [-0.3594, -0.2763] | 0.810 | pass |
| composition_gp | div_sd | 0.1123 | 0.1090 | [0.0968, 0.1222] | 0.263 | pass |
| composition_gp | div_mean | 0.5633 | 0.5636 | [0.5522, 0.5754] | 0.527 | pass |
| composition_gp_marginal | fst | 0.0156 | 0.0157 | [0.0127, 0.0190] | 0.530 | pass |
| composition_gp_marginal | decay | -0.3407 | -0.3223 | [-0.3634, -0.2716] | 0.777 | pass |
| composition_gp_marginal | div_sd | 0.1123 | 0.1095 | [0.0960, 0.1224] | 0.343 | pass |
| composition_gp_marginal | div_mean | 0.5633 | 0.5645 | [0.5514, 0.5773] | 0.573 | pass |

## Verdict

**The model reproduces every discriminating statistic.** Observed cultural F_ST, distance decay and diversity spread all fall inside the posterior predictive interval. Misspecification is therefore NOT the explanation for the treedepth saturation on real data, which leaves the flat-likelihood reading: rho is not identified above the design's resolving ceiling, and the posterior piling up against the top of its range is the data saying there is no decay within the basin. `spatial_share` is reportable; a length-scale NUMBER is not.

