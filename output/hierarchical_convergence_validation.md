# Hierarchical convergence model: validation battery

Config: full — real-panel fit 2000 draws x 4 chains (800 bootstrap reps for panel/seriation SEs); calibration fits 1000 draws x 4 chains over 3 synthetic cases (3 cases); prior-sensitivity refits 1000 draws x 4 chains over 3 prior settings.

## 1. p_convergence calibration (synthetic panels)

| case | true all-positive | P(all four > 0) | divergences |
|---|---|---|---|
| all_up | True | 1.000 | 1 |
| mixed | False | 0.000 | 0 |
| all_down | False | 0.000 | 1 |

## Real-panel fit (baseline)

- P(all four slopes > 0) = 0.715.
- MCMC health: max R-hat = nan, min ESS = 1984, divergences = 1.

## 2. Prior sensitivity (real panel)

| mu_sd | tau_sd | sigma_sd | P(all four > 0) |
|---|---|---|---|
| 0.5 | 0.5 | 0.25 | 0.770 |
| 1.0 | 1.0 | 0.5 | 0.726 |
| 2.0 | 2.0 | 1.0 | 0.655 |

- Verdict STABLE across priors (range 0.655-0.770).

## 3. Posterior-predictive check (real panel)

- 16/18 observed panel cells (88.9%) fall within the model's 95% posterior-predictive interval.

## 4. Cross-read vs analysis 07 (recomputed, not parsed from output/empirical_refined.md)

| signature | 07 recomputed OLS slope | model posterior mean slope | same sign |
|---|---|---|---|
| neutral_departure | +0.30512 | +0.576 | True |
| fst | +0.00076 | +0.272 | True |
| spatial_boundary | +0.51494 | +0.326 | True |
| seriation | +0.49578 | +0.489 | True |

- All panel-signature signs match between the recomputed 07 slopes and the model posterior: True.
