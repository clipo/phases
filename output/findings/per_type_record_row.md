# The record's per-class row

Produced by `analyses/per_type/drift_sensitivity_figure.py` from the curated 28-assemblage matrix and the two-cluster spatial partition (`a0_predictive_summaries.spatial_clusters`). F is the site-weighted method-of-moments F with binomial sampling variance removed; sites is the number of the 28 sites at which the class is present. Bands in the sensitivity findings are the 10th to 90th percentile of a cell's finite values under NumPy's linear quantile rule, endpoints inclusive, eligible when at least 45 of the records give a finite value.

| class | F | sites present | status | Spearman share vs size |
|---|---|---|---|---|
| Parkin_Punctated | 0.1646 | 28 | finite | +0.42 |
| Barton/Kent/MPI | 0.1172 | 28 | finite | +0.26 |
| Painted | 0.1986 | 28 | finite | -0.57 |
| Fortune_Noded | 0.0255 | 15 | finite | +0.56 |
| Ranch_Incised | 0.0396 | 22 | finite | -0.03 |
| Walls_Engraved | 0.0485 | 14 | finite | -0.30 |
| Wallace_Incised | 0.0054 | 7 | finite | +0.56 |
| Rhodes_Incised | 0.0267 | 15 | finite | -0.15 |
| Vernon_Paul_Applique | 0.0047 | 10 | finite | +0.43 |
| Hull_Engraved | 0.0057 | 8 | finite | -0.09 |

Site F_ST (sherd-weighted pool) 0.0479; site F_ST (site-mean pool) 0.1426; two-cluster F_ST 0.0062; richness 6.25.
