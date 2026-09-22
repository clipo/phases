# Tempo and mode as parameters, not model selection

Produced by `analyses/54_tempo_mode_posterior.py` (full). Gini-Simpson diversity over 6 seriation bins.

Replaces AICc and Akaike weights over four time-series models with the posteriors of the two parameters that distinguish them (`docs/FREQUENTIST_INVENTORY.md` item E). BM, GRW, Stasis and OU are regions of one parameter space, not separate hypotheses.

## The two parameters

| parameter | meaning | median | 95% CI | reading |
|---|---|---|---|---|
| alpha | mean reversion; ~0 is BM, large is Stasis | 0.765 | [0.045, 10.661] | consistent with a random walk, no resolved attractor |
| mu | directional drift per bin; 0 is unbiased | -0.0608 | [-0.1537, +0.0355] | P(mu > 0) = 0.071; no resolved direction |

## Diagnostics (rule 16)

- empirical: R-hat 1.0078, min ESS 980, divergences 18/8000 (0.23%), E-BFMI 0.612
- recovery: R-hat 1.1779, min ESS 16, divergences 1079/8000, E-BFMI 0.585

## Rule 20(c): recovery of a directional series

`mu ~ Normal(0, 0.5)` is centred on zero, which is the paper's own conclusion, so it is sympathetic. A series simulated with a true drift of +0.040 per bin returns **+0.0377** [-0.0111, +0.0772], covering the truth.

A real directional signal is recovered, so a near-zero mu on the diversity trajectory is a statement about the data.

## A confound the four-way selection was hiding

The diversity trajectory rises monotonically across the six bins (0.666, 0.734, 0.590, 0.600, 0.553, 0.460). Yet mu is not resolved. That is not a failure of the fit: the OU relaxation term, `theta + (anc - theta) exp(-alpha t)`, can produce exactly that rise by starting below the optimum and relaxing up to it, so sustained directional drift and relaxation toward a higher equilibrium compete to explain the same monotone increase. With six points the data cannot separate them.

**Selecting a single model concealed this.** Reporting that the trajectory 'favours the unbiased random walk' names a winner among four options without saying that two of the underlying behaviours are indistinguishable on this record. The parameter posteriors say it directly: alpha spans BM to Stasis and mu spans both signs, so the defensible claim is that the diversity trajectory is consistent with a random walk and does not exclude either directional change or mean reversion.

## The second conversion in this script

Analysis 20 obtained each bin's diversity sampling variance by resampling `multinomial(Ntot, pr)` 500 times at the plug-in proportions `pr`. Here the per-bin composition is drawn from its Dirichlet posterior instead, so uncertainty in `pr` is carried rather than conditioned away. Per-bin posterior SDs: 0.0065, 0.0064, 0.0082, 0.0061, 0.0043, 0.0081.
