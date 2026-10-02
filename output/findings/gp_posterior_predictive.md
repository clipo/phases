# Posterior predictive check: can the model reproduce the basin?

Produced by `analyses/02_spatial/04_posterior_predictive.R` on 2026-10-01.

Statistics chosen to discriminate: quantities the model was not fitted to directly and that the paper's claims rest on.

| model | statistic | observed | predictive median | 95% predictive interval | Bayesian p | verdict |
|---|---|---|---|---|---|---|
| composition_gp | fst | 0.0156 | 0.0156 | [0.0130, 0.0190] | 0.512 | pass |
| composition_gp | decay | -0.3600 | -0.3400 | [-0.3824, -0.2913] | 0.848 | pass |
| composition_gp | div_sd | 0.1123 | 0.1084 | [0.0952, 0.1206] | 0.258 | pass |
| composition_gp | div_mean | 0.5633 | 0.5648 | [0.5521, 0.5774] | 0.605 | pass |
| composition_gp_marginal | fst | 0.0156 | 0.0155 | [0.0123, 0.0192] | 0.468 | pass |
| composition_gp_marginal | decay | -0.3600 | -0.3383 | [-0.3834, -0.2957] | 0.835 | pass |
| composition_gp_marginal | div_sd | 0.1123 | 0.1084 | [0.0965, 0.1217] | 0.258 | pass |
| composition_gp_marginal | div_mean | 0.5633 | 0.5646 | [0.5532, 0.5770] | 0.588 | pass |

## Verdict

**The model reproduces every discriminating statistic.** Observed cultural F_ST, distance decay and diversity spread all fall inside the posterior predictive interval. Misspecification is therefore NOT the explanation for the treedepth saturation on real data, which leaves the flat-likelihood reading: rho is not identified above the design's resolving ceiling, and the posterior piling up against the top of its range is the data saying there is no decay within the basin. `spatial_share` is reportable; a length-scale NUMBER is not.

