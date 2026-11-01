# Softmax fields of variance F/p_c against Dirichlet drift

Produced by `analyses/probe_latent_groups/softmax_drift_check.py`: 400,000 draws (seed 1) of K = 10 independent Gaussian fields with mean log p_c and variance F/p_c at the record's pooled shares (`data/stan/basin_composition.json`), passed through the softmax; realized per-class moment F = var(p_c) / (mean(p_c) (1 - mean(p_c))) and the ratio of the realized mean share to the pooled share. To first order both would equal F and 1. Descriptive; the basis of D27.

| class | pooled share | F set 0.01: realized F, mean ratio | F set 0.05: realized F, mean ratio | F set 0.12: realized F, mean ratio |
|---|---|---|---|---|
| Parkin_Punctated | 0.5462 | 0.020, 0.93 | 0.175, 0.57 | 0.274, 0.31 |
| Barton/Kent/MPI | 0.3194 | 0.014, 0.93 | 0.108, 0.59 | 0.191, 0.33 |
| Painted | 0.0868 | 0.011, 0.97 | 0.069, 0.71 | 0.160, 0.50 |
| Fortune_Noded | 0.0106 | 0.021, 1.45 | 0.251, 3.24 | 0.489, 4.99 |
| Ranch_Incised | 0.0152 | 0.017, 1.27 | 0.180, 2.08 | 0.409, 2.90 |
| Walls_Engraved | 0.0098 | 0.023, 1.51 | 0.264, 3.60 | 0.509, 5.69 |
| Wallace_Incised | 0.0035 | 0.084, 3.51 | 0.522, 19.04 | 0.706, 31.28 |
| Rhodes_Incised | 0.0048 | 0.052, 2.49 | 0.443, 11.09 | 0.652, 18.48 |
| Vernon_Paul_Applique | 0.0024 | 0.138, 5.75 | 0.602, 35.60 | 0.758, 56.70 |
| Hull_Engraved | 0.0013 | 0.291, 17.41 | 0.716, 103.51 | 0.824, 152.25 |

The first-order matching fails at every F tried, for the common classes as well: the implied variances F/p_c run from 0.02 (Parkin_Punctated at F 0.01) to 94 logits^2 (Hull_Engraved at F 0.12), and the rare classes' dominate the normalization. The figures 1.48 and 1.21 in the item 1.1 entry (1/sqrt(1 - p_c)) are first-order values that no draw reproduces and are withdrawn.
