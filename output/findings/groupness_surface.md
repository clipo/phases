# Is groupness a gradient or an edge?

Basin phase set, 43 assemblages, 10 decorated classes. Turnover is the root sum of squares of the
weighted-least-squares gradient of the class proportions, Gaussian bandwidth 15 km,
on a 160 x 160 grid, masked where the local weight sum falls below 2.

## Do the phase boundaries lie on ridges of turnover?

Median turnover along the internal phase boundaries is **0.00866** per km,
from 109 points sampled at 1.5 km spacing.

**Against boundaries of the same construction drawn elsewhere** (200 partitions
carrying the phases' group sizes, seeds at random, assigned by exact minimum-cost
matching, each dissolved into territories the same way): median 0.01117 per km,
and the phase boundaries sit at the **3rd percentile**. This is the number to read.

**Against the whole mapped area** (median 0.01438 over 13,526 unmasked grid points): the 19th
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
| 26.3 | 2 | 0.0076 | 0.0136 | no | 13 |
| 21.3 | 3 | 0.0163 | 0.0135 | yes | 14 |
| 17.3 | 4 | 0.0207 | 0.0143 | yes | 3 |
| 12.7 | 5 | 0.0280 | 0.0151 | yes | 3 |
| 11.4 | 6 | 0.0420 | 0.0168 | yes | 3 |
| 10.1 | 7 | 0.0447 | 0.0169 | yes | 3 |
| 9.1 | 8 | 0.0412 | 0.0168 | yes | 1 |
| 8.6 | 9 | 0.0430 | 0.0172 | yes | 2 |
| 8.0 | 10 | 0.0473 | 0.0173 | yes | 1 |
| 7.7 | 11 | 0.0507 | 0.0178 | yes | 2 |
| 6.7 | 12 | 0.0512 | 0.0175 | yes | 1 |

The record sits inside the drift band at a radius of 26 km and leaves it by 21 km.

**The smallest-group column is the limit on this table.** Once a partition isolates a
single assemblage, that group matches its own profile exactly and contributes to
between-group variance for a reason that is arithmetic rather than archaeological.
Rows whose smallest group is 1 or 2 are reported but should not be read as evidence
of structure at that scale.

Figure written to fig13_groupness_surface.png and its siblings.
