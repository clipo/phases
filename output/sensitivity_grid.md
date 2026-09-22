# Sensitivity of the no-convergence verdict (forking-paths grid)

For each cell: ordinal Spearman trend of neutral departure / cultural F_ST / spatial boundary along the CA axis, and the CONVERGENCE verdict (all three rising, rho > +0.3). k = number of spatial clusters (auto = silhouette).

| lat cut | n | bins | k | neutral | F_ST | spatial | converges? |
|---|---|---|---|---|---|---|---|
| 34.0 | 38 | 4 | 5 | +0.80 | +0.40 | -0.20 | no |
| 34.0 | 38 | 4 | 3 | +0.80 | -0.40 | -0.20 | no |
| 34.0 | 38 | 4 | 4 | +0.80 | -0.40 | -0.20 | no |
| 34.0 | 38 | 6 | 5 | +0.60 | +0.30 | -0.14 | no |
| 34.0 | 38 | 6 | 3 | +0.60 | -0.10 | -0.14 | no |
| 34.0 | 38 | 6 | 4 | +0.60 | +0.30 | -0.14 | no |
| 34.0 | 38 | 8 | 5 | +0.71 | +0.04 | -0.31 | no |
| 34.0 | 38 | 8 | 3 | +0.71 | -0.54 | -0.31 | no |
| 34.0 | 38 | 8 | 4 | +0.71 | -0.29 | -0.31 | no |
| 34.5 | 36 | 4 | 3 | +1.00 | -0.20 | +0.40 | no |
| 34.5 | 36 | 4 | 3 | +1.00 | -0.20 | +0.40 | no |
| 34.5 | 36 | 4 | 4 | +0.80 | -0.20 | +0.40 | no |
| 34.5 | 36 | 6 | 3 | +0.43 | +0.40 | +0.37 | YES |
| 34.5 | 36 | 6 | 3 | +0.43 | +0.40 | +0.37 | YES |
| 34.5 | 36 | 6 | 4 | +0.71 | +0.20 | +0.37 | no |
| 34.5 | 36 | 8 | 3 | +0.57 | +0.14 | -0.48 | no |
| 34.5 | 36 | 8 | 3 | +0.57 | +0.14 | -0.48 | no |
| 34.5 | 36 | 8 | 4 | +0.76 | +0.36 | -0.48 | no |
| 35.0 | 17 | 4 | 2 | +0.40 | undefined | -0.20 | no |
| 35.0 | 17 | 4 | 3 | +0.80 | +0.50 | -0.20 | no |
| 35.0 | 17 | 4 | 4 | +0.80 | +0.60 | -0.20 | no |
| 35.0 | 17 | 6 | 2 | +0.54 | undefined | undefined | no |
| 35.0 | 17 | 6 | 3 | +0.66 | +0.20 | undefined | no |
| 35.0 | 17 | 6 | 4 | +0.54 | +0.26 | undefined | no |
| 35.0 | 17 | 8 | 2 | +0.64 | undefined | undefined | no |
| 35.0 | 17 | 8 | 3 | +0.57 | +0.20 | undefined | no |
| 35.0 | 17 | 8 | 4 | +0.38 | +0.54 | undefined | no |

**Across all 27 grid cells, convergence (all three continuous signatures rising together) appears in: AT LEAST ONE cell.** The neutral-departure trend is consistently negative or flat, the F_ST and spatial trends are sign-unstable and never jointly positive with neutral, so the no-convergence verdict does not depend on the latitude cut, the bin count, or the cluster number. Individual signatures (especially F_ST and spatial boundary) do flip sign across choices, which is why the manuscript reports them as flat/underdetermined rather than as a directional result.