# A0 prior- and posterior-predictive resemblance summaries (rule 20a)

Produced by `analyses/probe_latent_groups/a0_predictive_summaries.py` from the counts written by `refit_a0_basin.R` (200 draws each). Summaries computed with the operating-point screen's own functions (`47_revision_analysis.diversity`, plug-in `cultural_fst` on the screen's spatial clusters). Exploratory design validation, not evidence about the record.

Observed counts: 28 assemblages, 14101 sherds, counts sha256 `41bd753fe1c75396b09df992947a760de16739a9dcae0c7f5827f7adf351c7a7`. Spatial clusters: 2 (sizes 17, 11), labels sha256 `d574c257c8823ada15783c3a5733035c94fe488429a1ed7256863920dbde4777`, in `output/probe_latent_groups/spatial_clusters.json`.

| summary | observed | prior predictive median [5%, 95%] | fraction of prior draws below observed | posterior predictive median [5%, 95%] | fraction of posterior draws below observed |
|---|---|---|---|---|---|
| mean within-assemblage Gini-Simpson | 0.5633 | 0.1309 [0.0000, 0.4328] | 1.00 | 0.5637 [0.5532, 0.5733] | 0.47 |
| mean class richness | 6.2500 | 2.1786 [1.0000, 3.9714] | 1.00 | 6.4286 [6.1768, 6.7500] | 0.15 |
| pooled Gini-Simpson | 0.5916 | 0.2693 [0.0000, 0.6643] | 0.85 | 0.5911 [0.5861, 0.5983] | 0.56 |
| plug-in F_ST between the spatial clusters | 0.0062 | 0.0135 [0.0001, 0.2110] | 0.38 | 0.0063 [0.0051, 0.0077] | 0.47 |

Undefined draws, excluded from the medians and intervals above: 13 of 200 prior-predictive draws for plug-in F_ST between the spatial clusters (the pooled record holds a single class, so there is no diversity to partition).

Reading: the prior-predictive column says what the priors alone generate on this geometry and these sherd totals; a fraction near 0 or 1 means the observed value sits in a tail of the prior predictive. The posterior-predictive column is the ordinary predictive check of the fitted model on the same four summaries.
