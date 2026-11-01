# Hierarchical convergence model: validation battery

Config: full — real-panel fit 2000 draws x 4 chains (800 bootstrap reps for panel/seriation SEs); calibration fits 1000 draws x 4 chains over 3 synthetic cases (3 cases); prior-sensitivity refits 1000 draws x 4 chains over 3 prior settings.

## 1. p_convergence calibration (synthetic panels)

| case | true all-positive | P(all four > 0) | divergences |
|---|---|---|---|
| all_up | True | 1.000 | 0 |
| mixed | False | 0.000 | 0 |
| all_down | False | 0.000 | 0 |

## Real-panel fit (baseline)

- P(all four slopes > 0) = 0.000.
- MCMC health: max R-hat = 1.0018, min ESS = 1979, divergences = 0.

## 2. Prior sensitivity (real panel)

| mu_sd | tau_sd | sigma_sd | P(all four > 0) |
|---|---|---|---|
| 0.5 | 0.5 | 0.25 | 0.000 |
| 1.0 | 1.0 | 0.5 | 0.000 |
| 2.0 | 2.0 | 1.0 | 0.000 |

- Verdict STABLE across priors (range 0.000-0.000).

## 3. Posterior-predictive check (real panel)

- 18/18 observed panel cells (100.0%) fall within the model's 95% posterior-predictive interval.

## 4. Cross-read vs analysis 07 (recomputed, not parsed from output/empirical_refined.md)

| signature | 07 recomputed OLS slope | model posterior mean slope | same sign |
|---|---|---|---|
| neutral_departure | -0.19992 | -0.618 | True |
| fst | -0.01047 | -0.448 | True |
| spatial_boundary | +2.14117 | -0.348 | False |
| seriation | -0.64164 | -0.633 | True |

- All panel-signature signs match between the recomputed 07 slopes and the model posterior: False.
