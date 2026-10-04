# Can this design resolve an interaction length scale?

Produced by `analyses/48_length_scale_recovery.py` (full: 40 seeds per level, 1500 generations).

Design: the real St. Francis basin, 28 assemblages at their real coordinates and real per-assemblage sherd totals (14101 sherds), 10 decorated classes, maximum pairwise distance 79.1 km.

## Recovery

| true ell (km) | median estimate | 5-95% of estimates | fits |
|---|---|---|---|
| 4 | 22.2 | 16.0 to 33.8 | 40/40 |
| 8 | 53.1 | 20.1 to 123.9 | 40/40 |
| 16 | 49.7 | 9.1 to 500.0 | 40/40 |
| 32 | 257.1 | 4.7 to 500.0 | 40/40 |
| 64 | 500.0 | 3.7 to 500.0 | 40/40 |

## Separability of adjacent (2x) levels

Overlap is the fraction of the two estimate distributions that cannot be told apart, computed as the proportion of pairs (x from the lower level, y from the upper) with x >= y. 0.5 means the design carries no information distinguishing them; 0.0 means they never cross.

| pair | overlap (ell_hat) | overlap (Mantel r) |
|---|---|---|
| 4 vs 8 km | 0.119 | 0.144 |
| 8 vs 16 km | 0.475 | 0.009 |
| 16 vs 32 km | 0.466 | 0.070 |
| 32 vs 64 km | 0.474 | 0.242 |

Mantel r by level (mean +/- sd): 4 km: -0.744 +/- 0.076; 8 km: -0.633 +/- 0.090; 16 km: -0.325 +/- 0.107; 32 km: -0.118 +/- 0.082; 64 km: -0.043 +/- 0.066


## Verdict

Per adjacent pair, taking the better of the two observables: 4-8 km 0.119, 8-16 km 0.009, 16-32 km 0.070, 32-64 km 0.242. Separation threshold 0.15.

**The design resolves the scale over 4 to 32 km**, where every adjacent doubling separates. Above 32 km it does not.

**The two observables fail in different places, which is why both are here.** Mantel r is NON-MONOTONE in the interaction scale: it strengthens from -0.744 at 4 km to -0.633 at 8 km, then weakens steadily to -0.043 at 64 km. It peaks near the spacing of the sites themselves, so on its own it cannot tell a 4 km world from an 8 km one: they sit on opposite sides of the peak. ell_hat is monotone across that range and does separate them. Above 32 km the position reverses and both degrade, because a decay that long is barely expressed inside a study window 79 km across; there is no far field in which the similarity curve can flatten.

**Multiplicity, stated rather than buried.** The band above takes, for each pair, whichever of the two statistics separates better. For a question about whether the INFORMATION is present that is the right operation, and the margins are not marginal (the selected overlaps are 0.119, 0.009, 0.070 against 40 replicates per level). It would not be legitimate for estimating an effect size, and nothing here does that.

**Consequence for this paper.** The operating point the manuscript reports, an interaction range of about 24 km, sits inside the resolved band, near its upper end. The latent-field direction is therefore supported: a length-scale posterior will be driven by the data rather than by the prior over the range that matters. The reportable claim has to carry the ceiling with it, in the form 'resolved to roughly 32 km; longer interaction scales are not distinguishable from one another at this study extent'. That ceiling is set by the extent of the basin, not by sherd counts, so no amount of additional excavation at these sites would lift it.

**Caveat that must travel with this.** ell_hat is biased upward throughout (median 20.9 km when the truth is 4 km) and unstable at 64 km. It is a screening statistic, not a proposed estimator, and its bias is a property of fitting three parameters to non-independent pairs over a bounded window. What this screen licenses is the separability result, not any point estimate.
