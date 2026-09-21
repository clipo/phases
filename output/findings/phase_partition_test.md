# Are the phase boundaries where the differences are?

Basin phase set, 43 assemblages, 5 phases (Kent 10, Nodena 3, Parchman 4, Parkin 11, Walls 15).
Phase partition cultural F_ST **0.0259**, within-group inertia **161.8 km^2** per assemblage.

| partition ensemble | median F_ST | median inertia (km^2) | the phases sit at |
|---|---|---|---|
| labels shuffled, n = 2000 | 0.0091 | (geography destroyed) | above 1,994 of 2,000 draws |
| same sizes, boundaries at random, n = 1000 | 0.0277 | 424 | **33rd pct** (30th-36th) |
| same sizes, compactness matched, n = 1000 | 0.0290 | 335 | **22nd pct** (20th-25th) |

Percentiles are one-sided empirical percentiles of the ensemble, with the
2.5th-97.5th range over 2000 resamples of the same draw count in brackets.

## Did the compactness match hold?

The phases' inertia is 161.8 km^2 per assemblage. The random-boundary
ensemble's median is 424 (291-588); the compactness-matched ensemble's is 335 (226-461).
The match is one-sided: the comparison partitions are looser than the phases,
not tighter, because the phases sit essentially at the compactness optimum for
these group sizes.

Across the random-boundary draws, inertia and F_ST correlate at Spearman -0.05, and
across the compactness-matched draws at +0.23. So differentiation here is NOT
a simple function of how tight the groups are, and the looseness above does not
by itself bias the comparison. The superseded greedy generator did show such a
coupling (-0.35), which is a fact about that generator's fractured partitions
rather than about this statistic. Panel B plots both clouds.

## What the compactness search recovers

Of 1000 restarts, **6** (0.6 percent) returned the phase partition
exactly (adjusted Rand index 1.0); the mean ARI with the phases is 0.50.
The lowest inertia found anywhere in the search is 159.8 km^2 against the phases'
161.8, so the phases are within 1.2 percent of the most compact 5-way division
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

The drift comparison in panel A is a separate and stronger statement, and it is
the one the paper's argument rests on: differentiation exceeds what calibrated
spatial drift produces on this geography at every scale the test can resolve,
whatever partition is used to measure it.

Figure written to fig12_phase_partition.png and its siblings; ensembles in phase_partition_ensembles.csv.
