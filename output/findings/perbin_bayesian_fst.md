# Per-bin Bayesian F_ST, and the trend as a posterior

Produced by `analyses/50_perbin_bayesian_fst.py` (full). St. Francis basin, 28 assemblages, 2 spatial clusters, 6 seriation bins.

Replaces two frequentist objects with one generative one (`docs/FREQUENTIST_INVENTORY.md` items B and H): the rarefaction resampling interval on the reported trajectory, and the 800-draw assemblage-bootstrap standard errors that `analyses/40_hierarchical_convergence.py` feeds into its likelihood as if known.

## Per-bin Gini-Simpson F_ST posteriors

| bin | prior | median | 95% CI | posterior SD | plug-in | sherds | R-hat | min ESS | div |
|---|---|---|---|---|---|---|---|---|---|
| 0 | Beta(1,10) | (fewer than two clusters) | - | - | - | - | - | - | - |
| 1 | Beta(1,10) | 0.0015 | [0.0005, 0.0034] | 0.0007 | 0.0015 | 5871 | 1.0014 | 2590 | 0 |
| 2 | Beta(1,10) | 0.0058 | [0.0027, 0.0112] | 0.0022 | 0.0058 | 1874 | 1.0019 | 2752 | 0 |
| 3 | Beta(1,10) | 0.0307 | [0.0213, 0.0409] | 0.0050 | 0.0326 | 2338 | 1.0013 | 4749 | 0 |
| 4 | Beta(1,10) | (fewer than two clusters) | - | - | - | - | - | - | - |
| 5 | Beta(1,10) | 0.0030 | [0.0004, 0.0089] | 0.0023 | 0.0024 | 813 | 1.0011 | 3877 | 0 |
| 0 | uniform | (fewer than two clusters) | - | - | - | - | - | - | - |
| 1 | uniform | 0.0015 | [0.0005, 0.0034] | 0.0008 | 0.0015 | 5871 | 1.0016 | 2253 | 0 |
| 2 | uniform | 0.0058 | [0.0028, 0.0112] | 0.0022 | 0.0058 | 1874 | 1.0014 | 2501 | 0 |
| 3 | uniform | 0.0312 | [0.0219, 0.0411] | 0.0049 | 0.0326 | 2338 | 1.0008 | 5119 | 0 |
| 4 | uniform | (fewer than two clusters) | - | - | - | - | - | - | - |
| 5 | uniform | 0.0031 | [0.0004, 0.0094] | 0.0024 | 0.0024 | 813 | 1.0011 | 4601 | 0 |

## The trend, as a posterior slope

Slope of per-bin F_ST against standardized bin position. This is the quantity the manuscript currently reports as a Spearman rank correlation with a resampling interval.

| prior | slope median | 95% CI | P(slope > 0) | R-hat | min ESS | div |
|---|---|---|---|---|---|---|
| Beta(1,10) | +0.0009 | [-0.0341, +0.0330] | 0.535 | 1.0058 | 1445 | 2 |
| uniform | +0.0013 | [-0.0331, +0.0375] | 0.552 | 1.0009 | 1420 | 3 |

## Prior sensitivity (rule 20b)

The slope moves from +0.0009 under the Beta(1,10) prior to +0.0013 under uniform, a shift of 0.0004. P(slope > 0) moves from 0.535 to 0.552. Reported whether or not it is small, as the rule requires.

## Limitation, stated rather than buried

This is a two-stage analysis. Each bin's posterior is summarised by its median and standard deviation, and those summaries are carried into the trend model, so trend uncertainty is propagated through a normal approximation to each bin's posterior rather than exactly. A single joint model of bins and trend would be better and is not what this is. The improvement over what it replaces is nonetheless categorical: a posterior standard deviation has generative status, and a bootstrap standard deviation does not.
