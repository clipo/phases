# A0 prior- and posterior-predictive resemblance summaries (rule 20a)

Produced by `analyses/probe_latent_groups/a0_predictive_summaries.py` from the counts written by `refit_a0_basin.R` (200 draws each). Summaries computed with the operating-point screen's own functions (`47_revision_analysis.diversity`, plug-in `cultural_fst` on the screen's spatial clusters). Exploratory design validation, not evidence about the record.

Observed counts: 28 assemblages, 14101 sherds, counts sha256 `41bd753fe1c75396b09df992947a760de16739a9dcae0c7f5827f7adf351c7a7`. Spatial clusters: 2 (sizes 17, 11), labels sha256 `d574c257c8823ada15783c3a5733035c94fe488429a1ed7256863920dbde4777`, in `output/probe_latent_groups/spatial_clusters.json`.

| summary | observed | prior predictive median [5%, 95%] | fraction of prior draws below observed | posterior predictive median [5%, 95%] | fraction of posterior draws below observed |
|---|---|---|---|---|---|
| mean within-assemblage Gini-Simpson | 0.5633 | 0.4098 [0.2221, 0.6686] | 0.83 | 0.5651 [0.5547, 0.5748] | 0.38 |
| mean class richness | 6.2500 | 6.2143 [3.9286, 8.9286] | 0.51 | 6.5000 [6.2500, 6.8571] | 0.04 |
| pooled Gini-Simpson | 0.5916 | 0.7503 [0.3456, 0.8449] | 0.16 | 0.5921 [0.5858, 0.5979] | 0.46 |
| plug-in F_ST between the spatial clusters | 0.0062 | 0.0280 [0.0019, 0.1756] | 0.16 | 0.0064 [0.0049, 0.0078] | 0.43 |

Reading: the prior-predictive column says what the priors alone generate on this geometry and these sherd totals; a fraction near 0 or 1 means the observed value sits in a tail of the prior predictive. The posterior-predictive column is the ordinary predictive check of the fitted model on the same four summaries.
