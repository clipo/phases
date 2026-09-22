# How much of the reported F_ST is the partition we drew?

Produced by `analyses/59_partition_sensitivity.py` (full). F15 measurements M1 and M3. Basin, 28 assemblages, Beta(1,10) prior throughout. Not a test: this is the spread of the posterior across choices we made, reported as a spread.

**Pre-stated threshold.** The reported F_ST is 0.0063 with a 95 percent interval of width 0.0025. A choice that moves the median by more than that matters more than the data's own uncertainty. This was written down before the numbers were seen.

## M1. Sensitivity to k

| k | silhouette | F_ST median | 95% CI | R-hat | min ESS | divergences |
|---|---|---|---|---|---|---|
| 2 **(selected)** | 0.5288 | **0.0062** | [0.0051, 0.0076] | 1.0015 | 2059 | 0/6000 |
| 3 | 0.5016 | **0.0155** | [0.0133, 0.0180] | 1.0049 | 2085 | 0/6000 |
| 4 | 0.4692 | **0.0185** | [0.0160, 0.0212] | 1.0030 | 2521 | 0/6000 |
| 5 | 0.4701 | **0.0302** | [0.0267, 0.0342] | 1.0014 | 3671 | 0/6000 |
| 6 | 0.3830 | **0.0248** | [0.0217, 0.0283] | 1.0024 | 3616 | 0/6000 |

Across k = 2 to 6 the posterior median spans 0.0062 to 0.0302, a spread of **0.0240**, against a reported interval width of 0.0025. **That is 9.6 times the data's own uncertainty**, so the choice of k matters more than the evidence the data carry about F_ST.

## M3. Sensitivity to the k-means seed at k = 2

20 seeds. Posterior medians span 0.0062 to 0.0062, a spread of **0.0000** (0.00 times the reported interval width).

The seed is not a meaningful degree of freedom: the k-means helper takes the best of several initialisations, and at this k the partition is stable across seeds. One of the three degrees of freedom named in F15 can be struck.

## What this means for M2

k-dependence is the dominant term, so the partition-at-fixed-grain ensemble (M2) is worth running: the question becomes whether the specific partition is special among partitions of the same grain, or just one draw.

## Diagnostics (rule 16)

- k = 2: R-hat 1.0015, min ESS 2059, divergences 0/6000
- k = 3: R-hat 1.0049, min ESS 2085, divergences 0/6000
- k = 4: R-hat 1.0030, min ESS 2521, divergences 0/6000
- k = 5: R-hat 1.0014, min ESS 3671, divergences 0/6000
- k = 6: R-hat 1.0024, min ESS 3616, divergences 0/6000
- seed sweep: worst R-hat 1.0015, total divergences 0
