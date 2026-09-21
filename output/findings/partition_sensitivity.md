# How much of the reported F_ST is the partition we drew?

Produced by `analyses/59_partition_sensitivity.py` (full). F15 measurements M1 and M3. Basin, 43 assemblages, Beta(1,10) prior throughout. Not a test: this is the spread of the posterior across choices we made, reported as a spread.

**Pre-stated threshold.** The reported F_ST is 0.0179 with a 95 percent interval of width 0.0044. A choice that moves the median by more than that matters more than the data's own uncertainty. This was written down before the numbers were seen.

## M1. Sensitivity to k

| k | silhouette | F_ST median | 95% CI | R-hat | min ESS | divergences |
|---|---|---|---|---|---|---|
| 2 | 0.4555 | **0.0075** | [0.0065, 0.0086] | 1.0009 | 1753 | 0/6000 |
| 3 | 0.4173 | **0.0162** | [0.0144, 0.0181] | 1.0016 | 1709 | 0/6000 |
| 4 | 0.5022 | **0.0205** | [0.0185, 0.0226] | 1.0019 | 3075 | 0/6000 |
| 5 **(selected)** | 0.5665 | **0.0277** | [0.0255, 0.0300] | 1.0018 | 3568 | 0/6000 |
| 6 | 0.5138 | **0.0415** | [0.0386, 0.0445] | 1.0023 | 3253 | 0/6000 |

Across k = 2 to 6 the posterior median spans 0.0075 to 0.0415, a spread of **0.0340**, against a reported interval width of 0.0044. **That is 7.7 times the data's own uncertainty**, so the choice of k matters more than the evidence the data carry about F_ST.

## M3. Sensitivity to the k-means seed at k = 5

20 seeds. Posterior medians span 0.0277 to 0.0277, a spread of **0.0001** (0.01 times the reported interval width).

The seed is not a meaningful degree of freedom: the k-means helper takes the best of several initialisations, and at this k the partition is stable across seeds. One of the three degrees of freedom named in F15 can be struck.

## What this means for M2

k-dependence is the dominant term, so the partition-at-fixed-grain ensemble (M2) is worth running: the question becomes whether the specific partition is special among partitions of the same grain, or just one draw.

## Diagnostics (rule 16)

- k = 2: R-hat 1.0009, min ESS 1753, divergences 0/6000
- k = 3: R-hat 1.0016, min ESS 1709, divergences 0/6000
- k = 4: R-hat 1.0019, min ESS 3075, divergences 0/6000
- k = 5: R-hat 1.0018, min ESS 3568, divergences 0/6000
- k = 6: R-hat 1.0023, min ESS 3253, divergences 0/6000
- seed sweep: worst R-hat 1.0018, total divergences 0
