# A0 prior- and posterior-predictive resemblance summaries (rule 20a)

Produced by `analyses/probe_latent_groups/a0_predictive_summaries.py` from the counts written by `refit_a0_basin.R` (200 draws each). Summaries computed with the operating-point screen's own functions (`47_revision_analysis.diversity`, plug-in `cultural_fst` on the screen's spatial clusters). Exploratory design validation, not evidence about the record.

Observed counts: 28 assemblages, 14101 sherds, counts sha256 `41bd753fe1c75396b09df992947a760de16739a9dcae0c7f5827f7adf351c7a7`. Spatial clusters: 2 (sizes 17, 11), labels sha256 `d574c257c8823ada15783c3a5733035c94fe488429a1ed7256863920dbde4777`, in `output/probe_latent_groups/spatial_clusters.json`.

| summary | observed | prior predictive median [5%, 95%] | fraction of prior draws below observed | posterior predictive median [5%, 95%] | fraction of posterior draws below observed |
|---|---|---|---|---|---|
| mean within-assemblage Gini-Simpson | 0.5633 | 0.1234 [0.0000, 0.3616] | 0.99 | 0.5646 [0.5536, 0.5738] | 0.43 |
| mean class richness | 6.2500 | 2.1250 [1.0000, 4.3250] | 1.00 | 6.4643 [6.1768, 6.7518] | 0.12 |
| pooled Gini-Simpson | 0.5916 | 0.2258 [0.0000, 0.6794] | 0.83 | 0.5913 [0.5858, 0.5981] | 0.56 |
| plug-in F_ST between the spatial clusters | 0.0062 | 0.0146 [0.0001, 0.1455] | 0.37 | 0.0062 [0.0050, 0.0078] | 0.54 |

Undefined draws, excluded from the medians and intervals above: 20 of 200 prior-predictive draws for plug-in F_ST between the spatial clusters (the pooled record holds a single class, so there is no diversity to partition).

Reading: the prior-predictive column says what the priors alone generate on this geometry and these sherd totals; a fraction near 0 or 1 means the observed value sits in a tail of the prior predictive. The posterior-predictive column is the ordinary predictive check of the fitted model on the same four summaries.
