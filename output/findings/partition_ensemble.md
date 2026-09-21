# Is our partition special, or just one of many?

Produced by `analyses/60_partition_ensemble.py` (full). F15 measurement M2. Basin, 43 assemblages, k = 3, Beta(1,10) prior, 150 random contiguous (Voronoi) partitions.

Not a test: no null and no p-value. This locates our partition inside an ensemble of equally defensible ones and reports the spread.

## Result

| | cultural F_ST |
|---|---|
| **the k-means partition we use** | **0.0163** [0.0145, 0.0180] |
| ensemble median | 0.0179 |
| ensemble 5th to 95th percentile | 0.0063 to 0.0270 |
| ensemble full range | 0.0015 to 0.0312 |
| **our partition's position in the ensemble** | **33th percentile** |

## Reading

Our partition sits at the 33th percentile, which is unremarkable. At this grain the specific partition k-means found is **one draw from a family of equally defensible partitions**, and the F_ST defined on it is not a property the data single out.

Combined with M1, where the median moved 0.0200 across k against a data-driven interval width of 0.0044, the picture is that the reported F_ST is set mostly by **how finely the field is chopped** and hardly at all by which particular chopping is used. That is the same claim this paper makes about culture-historical phases, turned on its own instrument.

## What follows

The partition-free alternative is already in hand and needs no such defence: the spatial GP reproduces the observed F_ST (posterior predictive 0.0178 [0.0152, 0.0210] against 0.0179, Bayesian p = 0.470) while containing no partition at all, and reports `spatial_share` = 0.96 [0.89, 0.99]. The framing decision this measurement feeds is recorded in `docs/superpowers/plans/2026-09-02-f15-partition-dependence.md`.
