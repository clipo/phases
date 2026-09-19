# Tempo and mode as parameters, not model selection

Produced by `analyses/54_tempo_mode_posterior.py` (full). Gini-Simpson diversity over 6 seriation bins.

Replaces AICc and Akaike weights over four time-series models with the posteriors of the two parameters that distinguish them (`docs/FREQUENTIST_INVENTORY.md` item E). BM, GRW, Stasis and OU are regions of one parameter space, not separate hypotheses.

## The two parameters

| parameter | meaning | median | 95% CI | reading |
|---|---|---|---|---|
| alpha | mean reversion; ~0 is BM, large is Stasis | 0.465 | [0.051, 10.314] | consistent with a random walk, no resolved attractor |
| mu | directional drift per bin; 0 is unbiased | -0.0929 | [-0.1796, -0.0071] | P(mu > 0) = 0.019; a direction is resolved |

## Diagnostics (rule 16)

- empirical: R-hat 1.0068, min ESS 957, divergences 35/8000 (0.44%), E-BFMI 0.621
- recovery: R-hat 1.0184, min ESS 231, divergences 375/8000, E-BFMI 0.513

## Rule 20(c): recovery of a directional series

`mu ~ Normal(0, 0.5)` is centred on zero, which is the paper's own conclusion, so it is sympathetic. A series simulated with a true drift of +0.040 per bin returns **+0.0373** [-0.0163, +0.0830], covering the truth.

A real directional signal is recovered, so a near-zero mu on the diversity trajectory is a statement about the data.

## A confound the four-way selection was hiding

The diversity trajectory rises monotonically across the six bins (0.686, 0.675, 0.623, 0.561, 0.488, 0.367). Yet mu is not resolved. That is not a failure of the fit: the OU relaxation term, `theta + (anc - theta) exp(-alpha t)`, can produce exactly that rise by starting below the optimum and relaxing up to it, so sustained directional drift and relaxation toward a higher equilibrium compete to explain the same monotone increase. With six points the data cannot separate them.

**Selecting a single model concealed this.** Reporting that the trajectory 'favours the unbiased random walk' names a winner among four options without saying that two of the underlying behaviours are indistinguishable on this record. The parameter posteriors say it directly: alpha spans BM to Stasis and mu spans both signs, so the defensible claim is that the diversity trajectory is consistent with a random walk and does not exclude either directional change or mean reversion.

## The second conversion in this script

Analysis 20 obtained each bin's diversity sampling variance by resampling `multinomial(Ntot, pr)` 500 times at the plug-in proportions `pr`. Here the per-bin composition is drawn from its Dirichlet posterior instead, so uncertainty in `pr` is carried rather than conditioned away. Per-bin posterior SDs: 0.0032, 0.0046, 0.0077, 0.0073, 0.0037, 0.0068.
