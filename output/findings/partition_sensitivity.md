# How much of the reported F_ST is the partition we drew?

Produced by `analyses/59_partition_sensitivity.py` (full). F15 measurements M1 and M3. Basin, 29 assemblages, Beta(1,10) prior throughout. Not a test: this is the spread of the posterior across choices we made, reported as a spread.

**Pre-stated threshold.** The reported F_ST is 0.0179 with a 95 percent interval of width 0.0044. A choice that moves the median by more than that matters more than the data's own uncertainty. This was written down before the numbers were seen.

## M1. Sensitivity to k

| k | silhouette | F_ST median | 95% CI | R-hat | min ESS | divergences |
|---|---|---|---|---|---|---|
| 2 | 0.5528 | **0.0147** | [0.0128, 0.0167] | 1.0026 | 1886 | 0/6000 |
| 3 **(selected)** | 0.6008 | **0.0178** | [0.0157, 0.0200] | 1.0015 | 2444 | 0/6000 |
| 4 | 0.5051 | **0.0283** | [0.0256, 0.0313] | 1.0018 | 2883 | 0/6000 |
| 5 | 0.4935 | **0.0315** | [0.0287, 0.0345] | 1.0021 | 3423 | 0/6000 |
| 6 | 0.4774 | **0.0347** | [0.0319, 0.0376] | 1.0014 | 4200 | 0/6000 |

Across k = 2 to 6 the posterior median spans 0.0147 to 0.0347, a spread of **0.0201**, against a reported interval width of 0.0044. **That is 4.6 times the data's own uncertainty**, so the choice of k matters more than the evidence the data carry about F_ST.

## M3. Sensitivity to the k-means seed at k = 3

20 seeds. Posterior medians span 0.0178 to 0.0178, a spread of **0.0000** (0.00 times the reported interval width).

The seed is not a meaningful degree of freedom: the k-means helper takes the best of several initialisations, and at this k the partition is stable across seeds. One of the three degrees of freedom named in F15 can be struck.

## What this means for M2

k-dependence is the dominant term, so the partition-at-fixed-grain ensemble (M2) is worth running: the question becomes whether the specific partition is special among partitions of the same grain, or just one draw.

## Diagnostics (rule 16)

- k = 2: R-hat 1.0026, min ESS 1886, divergences 0/6000
- k = 3: R-hat 1.0015, min ESS 2444, divergences 0/6000
- k = 4: R-hat 1.0018, min ESS 2883, divergences 0/6000
- k = 5: R-hat 1.0021, min ESS 3423, divergences 0/6000
- k = 6: R-hat 1.0014, min ESS 4200, divergences 0/6000
- seed sweep: worst R-hat 1.0015, total divergences 0
