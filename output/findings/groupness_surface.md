# Is groupness a gradient or an edge?

Basin phase set, 28 assemblages, 10 decorated classes. Turnover is the root sum of squares of the
weighted-least-squares gradient of the class proportions, Gaussian bandwidth 15 km,
on a 220 x 220 grid, masked where the local weight sum falls below 2.

## Do the phase boundaries lie on ridges of turnover?

Median turnover along the internal phase boundaries is **0.01037** per km,
from 87 points sampled at 1.5 km spacing.

**Against boundaries of the same construction drawn elsewhere** (200 partitions
carrying the phases' group sizes, seeds at random, assigned by exact minimum-cost
matching, each dissolved into territories the same way): median 0.01239 per km,
and the phase boundaries sit at the **0th percentile**. This is the number to read.

**Against the whole mapped area** (median 0.01509 over 27,892 unmasked grid points): the 24th
percentile. **Do not read this one.** Turnover falls where local support rises
(Spearman -0.40; median 0.0159 per km where the local weight sum is 2 to 4, against
0.0108 where it is 6 to 9), because a regression fitted to few nearby assemblages
extrapolates steeply. Voronoi boundaries lie in the interior of the site distribution,
where support is highest, so ANY boundaries built this way score low against the whole
field. The matched comparison above removes that.

A boundary between interaction communities is a place where composition changes fast
over a short distance, so it would sit high against the matched comparison.

## At what spatial scale does drift stop accounting for it?

| group radius (km) | k | observed F_ST | drift 95% upper | above? | smallest group |
|---|---|---|---|---|---|
| 17.5 | 2 | 0.0062 | 0.0098 | no | 11 |
| 13.7 | 3 | 0.0156 | 0.0134 | yes | 9 |
| 11.7 | 4 | 0.0186 | 0.0157 | yes | 5 |
| 10.6 | 5 | 0.0307 | 0.0185 | yes | 4 |
| 9.1 | 6 | 0.0307 | 0.0187 | yes | 3 |
| 8.6 | 7 | 0.0308 | 0.0215 | yes | 1 |
| 7.3 | 8 | 0.0330 | 0.0241 | yes | 1 |
| 6.9 | 9 | 0.0470 | 0.0260 | yes | 1 |
| 6.0 | 10 | 0.0435 | 0.0279 | yes | 1 |
| 5.8 | 11 | 0.0446 | 0.0302 | yes | 1 |
| 4.9 | 12 | 0.0480 | 0.0283 | yes | 1 |

The record sits inside the drift band at a radius of 18 km and leaves it by 14 km.

**The smallest-group column is the limit on this table.** Once a partition isolates a
single assemblage, that group matches its own profile exactly and contributes to
between-group variance for a reason that is arithmetic rather than archaeological.
Rows whose smallest group is 1 or 2 are reported but should not be read as evidence
of structure at that scale.

Figure written to fig13_groupness_surface.png and its siblings.
