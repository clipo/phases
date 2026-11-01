# Are the phase boundaries where the differences are?

Basin phase set, 28 assemblages, 3 phases (Kent 11, Parkin 12, Walls 5).
Phase partition cultural F_ST **0.0123**, within-group inertia **197.5 km^2** per assemblage.

| partition ensemble | median F_ST | median inertia (km^2) | the phases sit at |
|---|---|---|---|
| labels shuffled, n = 2000 | 0.0048 | (geography destroyed) | above 1,882 of 2,000 draws |
| same sizes, boundaries at random, n = 1000 | 0.0116 | 278 | **57th pct** (54th-60th) |
| same sizes, compactness matched, n = 1000 | 0.0108 | 212 | **62nd pct** (59th-65th) |

Percentiles are one-sided empirical percentiles of the ensemble, with the
2.5th-97.5th range over 2000 resamples of the same draw count in brackets.

## Did the compactness match hold?

The phases' inertia is 197.5 km^2 per assemblage. The random-boundary
ensemble's median is 278 (157-443); the compactness-matched ensemble's is 212 (157-322).
41 percent of the compactness-matched partitions and 15 percent of the
random-boundary ones are tighter than the phases.

Across the random-boundary draws, inertia and F_ST correlate at Spearman -0.02, and
across the compactness-matched draws at -0.63. So differentiation here is NOT
a simple function of how tight the groups are, and the looseness above does not
by itself bias the comparison. The superseded greedy generator did show such a
coupling (-0.35), which is a fact about that generator's fractured partitions
rather than about this statistic. Panel B plots both clouds.

## What the compactness search recovers

Of 1000 restarts, **0** (0.0 percent) returned the phase partition
exactly (adjusted Rand index 1.0); the mean ARI with the phases is 0.53.
The lowest inertia found anywhere in the search is 157.1 km^2 against the phases'
197.5, so the phases are within 25.7 percent of the most compact 3-way division
of these assemblages at these group sizes.

## Reading

Beating shuffled labels says only that nearby assemblages resemble one another,
which isolation by distance produces on its own. The comparison that bears on
the phase lines is their position among the size-matched ensembles. A partition
marking interaction communities should sit high against them, because its
boundaries would be where the differences are and not merely where the gaps
between sites are; a partition that sits near the middle separates the pottery
about as well as any division of the map into groups of these sizes. The
compactness search says how far the phase division is from the tightest division
of these sizes, and how often a search that never sees a potsherd returns it.
This file states those positions and draws no conclusion from them.

The drift comparison in panel A is a separate statement: where the observed
differentiation lies against calibrated spatial drift at each number of
clusters. It is tabulated, with the share of runs reaching the observed value,
in scale_sweep.md (analysis 71), which this figure panel plots.

Figure written to fig12_phase_partition.png and its siblings; ensembles in phase_partition_ensembles.csv.
