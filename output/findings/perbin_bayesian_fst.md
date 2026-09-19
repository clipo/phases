# Per-bin Bayesian F_ST, and the trend as a posterior

Produced by `analyses/50_perbin_bayesian_fst.py` (full). St. Francis basin, 29 assemblages, 3 spatial clusters, 6 seriation bins.

Replaces two frequentist objects with one generative one (`docs/FREQUENTIST_INVENTORY.md` items B and H): the rarefaction resampling interval on the reported trajectory, and the 800-draw assemblage-bootstrap standard errors that `analyses/40_hierarchical_convergence.py` feeds into its likelihood as if known.

## Per-bin Gini-Simpson F_ST posteriors

| bin | prior | median | 95% CI | posterior SD | plug-in | sherds | R-hat | min ESS | div |
|---|---|---|---|---|---|---|---|---|---|
| 0 | Beta(1,10) | 0.0053 | [0.0021, 0.0102] | 0.0021 | 0.0056 | 1533 | 1.0020 | 4061 | 0 |
| 1 | Beta(1,10) | 0.0030 | [0.0007, 0.0073] | 0.0017 | 0.0028 | 1737 | 1.0040 | 1997 | 0 |
| 2 | Beta(1,10) | 0.0283 | [0.0168, 0.0419] | 0.0064 | 0.0292 | 1519 | 1.0013 | 3833 | 0 |
| 3 | Beta(1,10) | 0.1330 | [0.1174, 0.1484] | 0.0079 | 0.1355 | 3300 | 1.0013 | 7376 | 0 |
| 4 | Beta(1,10) | 0.0074 | [0.0052, 0.0102] | 0.0013 | 0.0075 | 14224 | 1.0009 | 2384 | 0 |
| 5 | Beta(1,10) | (fewer than two clusters) | - | - | - | - | - | - | - |
| 0 | uniform | 0.0056 | [0.0023, 0.0105] | 0.0021 | 0.0056 | 1533 | 1.0015 | 3776 | 0 |
| 1 | uniform | 0.0030 | [0.0007, 0.0074] | 0.0017 | 0.0028 | 1737 | 1.0010 | 2714 | 0 |
| 2 | uniform | 0.0286 | [0.0172, 0.0428] | 0.0065 | 0.0292 | 1519 | 1.0018 | 4368 | 0 |
| 3 | uniform | 0.1336 | [0.1180, 0.1491] | 0.0079 | 0.1355 | 3300 | 1.0006 | 6629 | 0 |
| 4 | uniform | 0.0075 | [0.0051, 0.0102] | 0.0013 | 0.0075 | 14224 | 1.0017 | 2532 | 0 |
| 5 | uniform | (fewer than two clusters) | - | - | - | - | - | - | - |

## The trend, as a posterior slope

Slope of per-bin F_ST against standardized bin position. This is the quantity the manuscript currently reports as a Spearman rank correlation with a resampling interval.

| prior | slope median | 95% CI | P(slope > 0) | R-hat | min ESS | div |
|---|---|---|---|---|---|---|
| Beta(1,10) | +0.0189 | [-0.0392, +0.0770] | 0.773 | 1.0003 | 4051 | 0 |
| uniform | +0.0185 | [-0.0401, +0.0760] | 0.764 | 1.0008 | 3599 | 0 |

## Prior sensitivity (rule 20b)

The slope moves from +0.0189 under the Beta(1,10) prior to +0.0185 under uniform, a shift of 0.0004. P(slope > 0) moves from 0.773 to 0.764. Reported whether or not it is small, as the rule requires.

## Limitation, stated rather than buried

This is a two-stage analysis. Each bin's posterior is summarised by its median and standard deviation, and those summaries are carried into the trend model, so trend uncertainty is propagated through a normal approximation to each bin's posterior rather than exactly. A single joint model of bins and trend would be better and is not what this is. The improvement over what it replaces is nonetheless categorical: a posterior standard deviation has generative status, and a bootstrap standard deviation does not.
