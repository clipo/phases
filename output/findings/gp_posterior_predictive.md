# Posterior predictive check: can the model reproduce the basin?

Produced by `analyses/02_spatial/04_posterior_predictive.R` on 2026-09-20.

Statistics chosen to discriminate: quantities the model was not fitted to directly and that the paper's claims rest on.

| model | statistic | observed | predictive median | 95% predictive interval | Bayesian p | verdict |
|---|---|---|---|---|---|---|
| composition_gp | fst | 0.0179 | 0.0492 | [0.0472, 0.0508] | 1.000 | **FAIL** |
| composition_gp | decay | -0.3899 | -0.0464 | [-0.0620, -0.0304] | 1.000 | **FAIL** |
| composition_gp | div_sd | 0.1378 | 0.1933 | [0.1886, 0.1974] | 1.000 | **FAIL** |
| composition_gp | div_mean | 0.5082 | 0.5443 | [0.5399, 0.5491] | 1.000 | **FAIL** |
| composition_gp_marginal | fst | 0.0179 | 0.0491 | [0.0474, 0.0509] | 1.000 | **FAIL** |
| composition_gp_marginal | decay | -0.3899 | -0.0453 | [-0.0600, -0.0308] | 1.000 | **FAIL** |
| composition_gp_marginal | div_sd | 0.1378 | 0.1933 | [0.1884, 0.1977] | 1.000 | **FAIL** |
| composition_gp_marginal | div_mean | 0.5082 | 0.5443 | [0.5395, 0.5487] | 1.000 | **FAIL** |

## Verdict

**At least one discriminating statistic falls outside the posterior predictive interval.** The model does not reproduce the basin on a quantity the paper's claims depend on. No length scale and no variance share from this model is reportable until that is understood, whatever the convergence diagnostics say.

