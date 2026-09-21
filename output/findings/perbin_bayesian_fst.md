# Per-bin Bayesian F_ST, and the trend as a posterior

Produced by `analyses/50_perbin_bayesian_fst.py` (full). St. Francis basin, 43 assemblages, 5 spatial clusters, 6 seriation bins.

Replaces two frequentist objects with one generative one (`docs/FREQUENTIST_INVENTORY.md` items B and H): the rarefaction resampling interval on the reported trajectory, and the 800-draw assemblage-bootstrap standard errors that `analyses/40_hierarchical_convergence.py` feeds into its likelihood as if known.

## Per-bin Gini-Simpson F_ST posteriors

| bin | prior | median | 95% CI | posterior SD | plug-in | sherds | R-hat | min ESS | div |
|---|---|---|---|---|---|---|---|---|---|
| 0 | Beta(1,10) | 0.0103 | [0.0056, 0.0166] | 0.0028 | 0.0112 | 2373 | 1.0033 | 2310 | 0 |
| 1 | Beta(1,10) | 0.0398 | [0.0314, 0.0486] | 0.0044 | 0.0413 | 2829 | 1.0015 | 4343 | 0 |
| 2 | Beta(1,10) | 0.0830 | [0.0716, 0.0946] | 0.0058 | 0.0860 | 4585 | 1.0014 | 4626 | 0 |
| 3 | Beta(1,10) | 0.0105 | [0.0079, 0.0135] | 0.0014 | 0.0108 | 8650 | 1.0015 | 2697 | 0 |
| 4 | Beta(1,10) | 0.0043 | [0.0030, 0.0060] | 0.0008 | 0.0043 | 10209 | 1.0009 | 3356 | 0 |
| 5 | Beta(1,10) | 0.0164 | [0.0130, 0.0203] | 0.0019 | 0.0165 | 7745 | 1.0007 | 5113 | 0 |
| 0 | uniform | 0.0104 | [0.0057, 0.0166] | 0.0028 | 0.0112 | 2373 | 1.0012 | 2598 | 0 |
| 1 | uniform | 0.0402 | [0.0318, 0.0489] | 0.0044 | 0.0413 | 2829 | 1.0015 | 3902 | 0 |
| 2 | uniform | 0.0834 | [0.0719, 0.0949] | 0.0059 | 0.0860 | 4585 | 1.0008 | 4934 | 0 |
| 3 | uniform | 0.0106 | [0.0080, 0.0136] | 0.0014 | 0.0108 | 8650 | 1.0018 | 2567 | 0 |
| 4 | uniform | 0.0043 | [0.0030, 0.0060] | 0.0008 | 0.0043 | 10209 | 1.0011 | 3085 | 0 |
| 5 | uniform | 0.0165 | [0.0130, 0.0203] | 0.0019 | 0.0165 | 7745 | 1.0012 | 6166 | 0 |

## The trend, as a posterior slope

Slope of per-bin F_ST against standardized bin position. This is the quantity the manuscript currently reports as a Spearman rank correlation with a resampling interval.

| prior | slope median | 95% CI | P(slope > 0) | R-hat | min ESS | div |
|---|---|---|---|---|---|---|
| Beta(1,10) | -0.0071 | [-0.0426, +0.0264] | 0.313 | 1.0012 | 3260 | 0 |
| uniform | -0.0069 | [-0.0421, +0.0282] | 0.324 | 1.0025 | 3128 | 0 |

## Prior sensitivity (rule 20b)

The slope moves from -0.0071 under the Beta(1,10) prior to -0.0069 under uniform, a shift of 0.0002. P(slope > 0) moves from 0.313 to 0.324. Reported whether or not it is small, as the rule requires.

## Limitation, stated rather than buried

This is a two-stage analysis. Each bin's posterior is summarised by its median and standard deviation, and those summaries are carried into the trend model, so trend uncertainty is propagated through a normal approximation to each bin's posterior rather than exactly. A single joint model of bins and trend would be better and is not what this is. The improvement over what it replaces is nonetheless categorical: a posterior standard deviation has generative status, and a bootstrap standard deviation does not.
