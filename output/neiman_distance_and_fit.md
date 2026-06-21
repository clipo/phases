# Neiman interassemblage-distance model and Frequency Increment Test

## B. Neiman interassemblage-distance neutral model

Squared Euclidean distance d_ij^2 = sum_k (p_ik - p_jk)^2 over 10 decorated types, 29 basin assemblages.

- **Diversity-distance relationship** (Neiman's convergence-of-two-measures precedent): within-assemblage diversity (t_E) vs mean interassemblage distance, Spearman rho = -0.28 (p = 0.145). Neiman (1995, Fig. 7) found rho = -0.62 under drift+innovation in one interacting population; a negative relationship here is consistent with a single drifting field rather than assorted, bounded groups.

- **Divergence trajectory**: mean within-bin pairwise distance across the 6 CA bins trends with seriation position at Spearman rho = +0.37 (p = 0.468). A rising trend would indicate growing between-assemblage divergence (fragmentation); a flat or falling trend indicates no growing divergence. Per-bin mean d^2: 0:0.014, 1:0.013, 2:0.408, 3:0.186, 4:0.169, 5:0.055.

- **Between- vs within-cluster distance** (k = 3 spatial clusters): within-cluster mean d^2 = 0.219, between-cluster = 0.248, ratio = 1.13. A ratio near 1 means spatial clusters are no more internally similar than the basin at large (no bounded ceramic groups); a large ratio would mark bounded groups.

## D. Frequency Increment Test (Feder, Kryazhimskiy & Plotkin 2014)

Per decorated type: rescaled increments Y_i = (v_i - v_{i-1}) / sqrt(2 v_{i-1}(1-v_{i-1}) dt) along the ordered bins (dt = 1); a one-sample t-test of mean(Y) = 0 tests neutral drift. Increments at v in {0,1} are dropped.

| type | n increments | mean Y | t p-value | departs neutral? |
|---|---|---|---|---|
| Parkin_Punctated | 5 | -0.13 | 0.021 | YES |
| Barton/Kent/MPI | 5 | +0.03 | 0.552 | no |
| Painted | 5 | +0.18 | 0.005 | YES |
| Fortune_Noded | 5 | -0.01 | 0.369 | no |
| Ranch_Incised | 5 | +0.03 | 0.606 | no |
| Walls_Engraved | 4 | +0.03 | 0.187 | no |
| Wallace_Incised | 2 | - | - | (too few) |
| Rhodes_Incised | 4 | +0.00 | 0.902 | no |
| Vernon_Paul_Applique | 3 | +0.00 | 0.960 | no |
| Hull_Engraved | 4 | +0.01 | 0.571 | no |

**2 of 9 testable types depart from neutral drift at p < 0.05.** A predominantly neutral result corroborates the static neutrality signature: the decorated types behave as drifting neutral variants along the sequence, not as markers under conformist or anti-conformist bias. (With six bins per type the per-type test has low power; read the count, not individual p-values.)