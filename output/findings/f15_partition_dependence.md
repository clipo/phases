# F15: the reported F_ST is set by the grain, not by the partition

Combined verdict on the three measurements planned in
`docs/superpowers/plans/2026-09-02-f15-partition-dependence.md`. Sources:
`output/findings/partition_sensitivity.md` (M1, M3) and
`output/findings/partition_ensemble.md` (M2). Basin, 29 assemblages, Beta(1,10)
prior throughout, all fits R-hat below 1.003 with zero divergences.

**Pre-stated yardstick**, written before any of these numbers were seen: the
reported F_ST is 0.0179 with a 95 percent credible interval of width 0.0044. A
researcher choice that moves the posterior median by more than that matters more
than the data's own uncertainty.

## The three degrees of freedom, measured

| choice | effect on the posterior median | verdict |
|---|---|---|
| **k-means seed** (M3, 20 seeds) | spread **0.0000** | irrelevant |
| **which partition at fixed k** (M2, 150 random contiguous partitions) | ours sits at the **51st percentile**, ensemble median 0.0178 against our 0.0179 | irrelevant |
| **k, the number of clusters** (M1) | spread **0.0200**, monotone | **dominant, 4.5x the data's own uncertainty** |

### M1 in full

| k | silhouette | F_ST median | 95% CI |
|---|---|---|---|
| 2 | 0.5528 | 0.0147 | [0.0128, 0.0167] |
| **3 (selected)** | **0.6008** | **0.0178** | [0.0157, 0.0200] |
| 4 | 0.5051 | 0.0283 | [0.0256, 0.0313] |
| 5 | 0.4935 | 0.0315 | [0.0287, 0.0345] |
| 6 | 0.4774 | 0.0347 | [0.0319, 0.0376] |

F_ST rises **monotonically** with k, and the intervals at k = 2 and k = 6 are
**entirely disjoint**. The number the paper reports is a statement about k = 3.

## What this means

Two of the three researcher degrees of freedom named in F15 can be struck: the
seed does nothing, and the particular partition does nothing. What remains is
the one that matters most.

**The reported F_ST is set almost entirely by how finely the field is chopped,
and not at all by where it is chopped.** A quantity with that behaviour is
measuring the chopping. This is the same claim the paper makes about
culture-historical phases, arriving at its own instrument, and it reproduces
independently what `../mataa` found when it ran the equivalent check on its own
imposed partitions, which is why that project retired the partition as its
method.

## What it does NOT mean

**The paper's conclusion is not threatened.** Three things stand untouched:

1. **The partition-free result says the same thing, more strongly.** The spatial
   GP contains no partition at all, reproduces the observed value in posterior
   predictive check (0.0178 [0.0152, 0.0210] against 0.0179, Bayesian p = 0.470),
   and reports `spatial_share` = 0.96 [0.89, 0.99]: almost all assemblage-level
   compositional variance is spatially structured, leaving almost none of the
   site-specific character bounded groups produce.
2. **The drift-versus-groups comparison scores simulated and real data through
   the same partition**, so grain-dependence cancels in the comparison even
   though it does not cancel in the absolute value. That needs stating in the
   text rather than assuming, and it is why the positive control is unaffected.
3. **The seriation and settlement lines involve no partition** and are untouched.

## Disposition

The finding supports leading with the partition-free result and demoting the
partition-based F_ST to corroboration reported with its grain-dependence stated.
That is a framing decision for the author, not a code change, and it is the last
open item in `docs/CODE_REVIEW_2026-08-31.md`.

What must change regardless of that decision: the manuscript currently reports
F_ST without saying that the value is a function of k. Under rule 6 the
operating point has to ride along with the quantity, and "k = 3, chosen by
silhouette" is part of that operating point.
