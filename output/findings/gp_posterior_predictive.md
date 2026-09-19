# Posterior predictive check: can the model reproduce the basin?

Produced by `analyses/02_spatial/04_posterior_predictive.R` on 2026-08-31.

Statistics chosen to discriminate: quantities the model was not fitted to directly and that the paper's claims rest on.

| model | statistic | observed | predictive median | 95% predictive interval | Bayesian p | verdict |
|---|---|---|---|---|---|---|
| composition_gp | fst | 0.0179 | 0.0178 | [0.0152, 0.0210] | 0.470 | pass |
| composition_gp | decay | -0.3267 | -0.3210 | [-0.3778, -0.2695] | 0.575 | pass |
| composition_gp | div_sd | 0.1444 | 0.1391 | [0.1263, 0.1538] | 0.233 | pass |
| composition_gp | div_mean | 0.5097 | 0.5098 | [0.4972, 0.5218] | 0.517 | pass |
| composition_gp_marginal | fst | 0.0179 | 0.0180 | [0.0152, 0.0212] | 0.512 | pass |
| composition_gp_marginal | decay | -0.3267 | -0.3233 | [-0.3699, -0.2652] | 0.557 | pass |
| composition_gp_marginal | div_sd | 0.1444 | 0.1396 | [0.1264, 0.1536] | 0.242 | pass |
| composition_gp_marginal | div_mean | 0.5097 | 0.5094 | [0.4952, 0.5217] | 0.487 | pass |

## Verdict

**The model reproduces every discriminating statistic.** Observed cultural F_ST, distance decay and diversity spread all fall inside the posterior predictive interval. Misspecification is therefore NOT the explanation for the treedepth saturation on real data, which leaves the flat-likelihood reading: rho is not identified above the design's resolving ceiling, and the posterior piling up against the top of its range is the data saying there is no decay within the basin. `spatial_share` is reportable; a length-scale NUMBER is not.

