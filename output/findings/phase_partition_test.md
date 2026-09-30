# Are the phase boundaries where the differences are?

Basin phase set, 28 assemblages, 3 phases (Kent 8, Parkin 11, Walls 9).
Phase partition cultural F_ST **0.0110**, within-group inertia **189.9 km^2** per assemblage.

| partition ensemble | median F_ST | median inertia (km^2) | the phases sit at |
|---|---|---|---|
| labels shuffled, n = 2000 | 0.0047 | (geography destroyed) | above 1,832 of 2,000 draws |
| same sizes, boundaries at random, n = 1000 | 0.0154 | 269 | **18th pct** (15th-20th) |
| same sizes, compactness matched, n = 1000 | 0.0171 | 223 | **5th pct** (4th-6th) |

Percentiles are one-sided empirical percentiles of the ensemble, with the
2.5th-97.5th range over 2000 resamples of the same draw count in brackets.

## Did the compactness match hold?

The phases' inertia is 189.9 km^2 per assemblage. The random-boundary
ensemble's median is 269 (194-487); the compactness-matched ensemble's is 223 (190-267).
The match is one-sided: the comparison partitions are looser than the phases,
not tighter, because the phases sit essentially at the compactness optimum for
these group sizes.

Across the random-boundary draws, inertia and F_ST correlate at Spearman -0.02, and
across the compactness-matched draws at +0.45. So differentiation here is NOT
a simple function of how tight the groups are, and the looseness above does not
by itself bias the comparison. The superseded greedy generator did show such a
coupling (-0.35), which is a fact about that generator's fractured partitions
rather than about this statistic. Panel B plots both clouds.

## What the compactness search recovers

Of 1000 restarts, **156** (15.6 percent) returned the phase partition
exactly (adjusted Rand index 1.0); the mean ARI with the phases is 0.73.
The lowest inertia found anywhere in the search is 189.9 km^2 against the phases'
189.9, so the phases are within 0.0 percent of the most compact 3-way division
of these assemblages at these group sizes.

## Reading

The phases are not arbitrary lines; they are where generations of workers
thought they saw structure. What the comparison shows is what that perception
was tracking. Beating shuffled labels says only that nearby assemblages
resemble one another, which isolation by distance produces on its own. Sitting
inside the size- and compactness-matched distribution says the particular
placement of these boundaries carries little information about where ceramic
differences lie beyond the fact that they enclose compact groups. And a search
that minimises inertia and never sees a potsherd rediscovers the phase
partition itself, from random starts, at the rate reported above.

A partition marking interaction communities should behave differently: it
should sit high against compactness-matched cuts, because the boundaries would
be where the differences are, not merely where the gaps between sites are.

The drift comparison in panel A is a separate statement: where the observed
differentiation lies against calibrated spatial drift at each number of
clusters. It is tabulated, with the share of runs reaching the observed value,
in scale_sweep.md (analysis 71), which this figure panel plots.

Figure written to fig12_phase_partition.png and its siblings; ensembles in phase_partition_ensembles.csv.
