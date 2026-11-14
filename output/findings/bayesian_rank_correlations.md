# Bayesian rank correlations, replacing the "not significant" claims

Produced by `analyses/51_bayesian_rank_correlations.py` (full). Data preparation copied unchanged from `analyses/13_neiman_distance_and_fit.py`, so any movement is attributable to the inference and not to the inputs.

Estimand: the correlation of van der Waerden normal scores, the rank-based Bayesian analogue of Spearman's rho, with a symmetric Uniform(-1, 1) prior (rule 20: the prior leans neither way, so the burden of defence is light).

| claim | n | frequentist rho | p | **posterior median** | **95% CI** | **P(r > 0)** |
|---|---|---|---|---|---|---|
| divergence trajectory (MAIN TEXT) | 6 | +0.029 | 0.957 | **+0.087** | [-0.531, +0.634] | 0.597 |
| diversity-distance (SUPPLEMENTAL) | 28 | -0.176 | 0.371 | **-0.140** | [-0.435, +0.204] | 0.218 |

## Diagnostics (rule 16, every fit)

- divergence trajectory (MAIN TEXT): r = +0.087 [-0.531, +0.634], P(r > 0) = 0.597, n = 6 | R-hat 1.0012, bulk ESS 2626, tail ESS 2467, divergences 0/8000 (0.00%), treedepth>=10 0, E-BFMI 1.112
- diversity-distance (SUPPLEMENTAL): r = -0.140 [-0.435, +0.204], P(r > 0) = 0.218, n = 28 | R-hat 1.0012, bulk ESS 3054, tail ESS 3090, divergences 0/8000 (0.00%), treedepth>=10 0, E-BFMI 1.138

## Reading

The frequentist and Bayesian point estimates agree closely, as they should: the same rank information is being summarised. What changes is what can be said about it.

For the **divergence trajectory**, the manuscript currently says "no divergence trend along the sequence (Spearman's rho = +0.37, not significant)". With six bins the posterior is very wide and straddles zero comfortably. The defensible statement is that the data do not resolve the direction of any trend, which is weaker than "no divergence trend" and stronger than nothing: it says the record cannot distinguish growing divergence from none, rather than asserting there is none.

This is the substantive difference rule 18 is after. "Not significant" over six points reads as evidence of absence. A posterior that spans most of the interval says plainly that six points cannot settle it.
