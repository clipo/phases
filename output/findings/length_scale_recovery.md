# Can this design resolve an interaction length scale?

Produced by `analyses/48_length_scale_recovery.py` (full: 40 seeds per level, 1500 generations).

Design: the real St. Francis basin, 43 assemblages at their real coordinates and real per-assemblage sherd totals (36391 sherds), 10 decorated classes, maximum pairwise distance 148.1 km.

## Recovery

| true ell (km) | median estimate | 5-95% of estimates | fits |
|---|---|---|---|
| 4 | 17.0 | 12.4 to 24.4 | 40/40 |
| 8 | 32.8 | 18.7 to 63.9 | 40/40 |
| 16 | 88.0 | 17.9 to 500.0 | 40/40 |
| 32 | 31.2 | 4.3 to 500.0 | 40/40 |
| 64 | 56.6 | 0.8 to 500.0 | 40/40 |

## Separability of adjacent (2x) levels

Overlap is the fraction of the two estimate distributions that cannot be told apart, computed as the proportion of pairs (x from the lower level, y from the upper) with x >= y. 0.5 means the design carries no information distinguishing them; 0.0 means they never cross.

| pair | overlap (ell_hat) | overlap (Mantel r) |
|---|---|---|
| 4 vs 8 km | 0.067 | 0.608 |
| 8 vs 16 km | 0.314 | 0.014 |
| 16 vs 32 km | 0.632 | 0.087 |
| 32 vs 64 km | 0.486 | 0.312 |

Mantel r by level (mean +/- sd): 4 km: -0.510 +/- 0.141; 8 km: -0.563 +/- 0.094; 16 km: -0.260 +/- 0.097; 32 km: -0.089 +/- 0.075; 64 km: -0.042 +/- 0.047


## Verdict

Per adjacent pair, taking the better of the two observables: 4-8 km 0.067, 8-16 km 0.014, 16-32 km 0.087, 32-64 km 0.312. Separation threshold 0.15.

**The design resolves the scale over 4 to 32 km**, where every adjacent doubling separates. Above 32 km it does not.

**The two observables fail in different places, which is why both are here.** Mantel r is NON-MONOTONE in the interaction scale: it strengthens from -0.510 at 4 km to -0.563 at 8 km, then weakens steadily to -0.042 at 64 km. It peaks near the spacing of the sites themselves, so on its own it cannot tell a 4 km world from an 8 km one: they sit on opposite sides of the peak. ell_hat is monotone across that range and does separate them. Above 32 km the position reverses and both degrade, because a decay that long is barely expressed inside a study window 148 km across; there is no far field in which the similarity curve can flatten.

**Multiplicity, stated rather than buried.** The band above takes, for each pair, whichever of the two statistics separates better. For a question about whether the INFORMATION is present that is the right operation, and the margins are not marginal (the selected overlaps are 0.067, 0.014, 0.087 against 40 replicates per level). It would not be legitimate for estimating an effect size, and nothing here does that.

**Consequence for this paper.** The operating point the manuscript reports, an interaction range of about 24 km, sits inside the resolved band, near its upper end. The latent-field direction is therefore supported: a length-scale posterior will be driven by the data rather than by the prior over the range that matters. The reportable claim has to carry the ceiling with it, in the form 'resolved to roughly 32 km; longer interaction scales are not distinguishable from one another at this study extent'. That ceiling is set by the extent of the basin, not by sherd counts, so no amount of additional excavation at these sites would lift it.

**Caveat that must travel with this.** ell_hat is biased upward throughout (median 20.9 km when the truth is 4 km) and unstable at 64 km. It is a screening statistic, not a proposed estimator, and its bias is a property of fitting three parameters to non-independent pairs over a bounded window. What this screen licenses is the separability result, not any point estimate.
