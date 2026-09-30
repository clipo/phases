# Is groupness a gradient or an edge?

Basin phase set, 28 assemblages, 10 decorated classes. Turnover is the root sum of squares of the
weighted-least-squares gradient of the class proportions, Gaussian bandwidth 15 km,
on a 220 x 220 grid, masked where the local weight sum falls below 2.

## Do the phase boundaries lie on ridges of turnover?

Median turnover along the internal phase boundaries is **0.01052** per km,
from 87 points sampled at 1.5 km spacing.

**Against boundaries of the same construction drawn elsewhere** (200 partitions
carrying the phases' group sizes, seeds at random, assigned by exact minimum-cost
matching, each dissolved into territories the same way): median 0.01242 per km,
and the phase boundaries sit at the **0th percentile**. This is the number to read.

**Against the whole mapped area** (median 0.01509 over 27,885 unmasked grid points): the 24th
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
| 17.6 | 2 | 0.0062 | 0.0124 | no | 11 |
| 13.6 | 3 | 0.0156 | 0.0135 | yes | 9 |
| 11.7 | 4 | 0.0186 | 0.0146 | yes | 5 |
| 10.0 | 5 | 0.0206 | 0.0148 | yes | 1 |
| 8.7 | 6 | 0.0327 | 0.0152 | yes | 1 |
| 8.1 | 7 | 0.0345 | 0.0161 | yes | 1 |
| 7.1 | 8 | 0.0346 | 0.0165 | yes | 1 |
| 6.5 | 9 | 0.0426 | 0.0183 | yes | 1 |
| 5.8 | 10 | 0.0441 | 0.0185 | yes | 1 |
| 5.1 | 11 | 0.0488 | 0.0180 | yes | 1 |
| 4.5 | 12 | 0.0529 | 0.0193 | yes | 1 |

The record sits inside the drift band at a radius of 18 km and leaves it by 14 km.

**The smallest-group column is the limit on this table.** Once a partition isolates a
single assemblage, that group matches its own profile exactly and contributes to
between-group variance for a reason that is arithmetic rather than archaeological.
Rows whose smallest group is 1 or 2 are reported but should not be read as evidence
of structure at that scale.

Figure written to fig13_groupness_surface.png and its siblings.
