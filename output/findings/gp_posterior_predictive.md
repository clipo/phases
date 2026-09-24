# Posterior predictive check: can the model reproduce the basin?

Produced by `analyses/02_spatial/04_posterior_predictive.R` on 2026-09-23.

Statistics chosen to discriminate: quantities the model was not fitted to directly and that the paper's claims rest on.

| model | statistic | observed | predictive median | 95% predictive interval | Bayesian p | verdict |
|---|---|---|---|---|---|---|
| composition_gp | fst | 0.0156 | 0.0157 | [0.0128, 0.0186] | 0.545 | pass |
| composition_gp | decay | -0.3410 | -0.3213 | [-0.3628, -0.2762] | 0.772 | pass |
| composition_gp | div_sd | 0.1123 | 0.1092 | [0.0970, 0.1216] | 0.310 | pass |
| composition_gp | div_mean | 0.5633 | 0.5647 | [0.5512, 0.5772] | 0.573 | pass |
| composition_gp_marginal | fst | 0.0156 | 0.0157 | [0.0128, 0.0188] | 0.525 | pass |
| composition_gp_marginal | decay | -0.3410 | -0.3205 | [-0.3668, -0.2735] | 0.782 | pass |
| composition_gp_marginal | div_sd | 0.1123 | 0.1087 | [0.0962, 0.1219] | 0.295 | pass |
| composition_gp_marginal | div_mean | 0.5633 | 0.5641 | [0.5502, 0.5745] | 0.570 | pass |

## Verdict

**The model reproduces every discriminating statistic.** Observed cultural F_ST, distance decay and diversity spread all fall inside the posterior predictive interval. Misspecification is therefore NOT the explanation for the treedepth saturation on real data, which leaves the flat-likelihood reading: rho is not identified above the design's resolving ceiling, and the posterior piling up against the top of its range is the data saying there is no decay within the basin. `spatial_share` is reportable; a length-scale NUMBER is not.

