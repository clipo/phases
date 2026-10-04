# What recovers the phase scheme: the pottery, or the map?

Basin phase set, 28 assemblages, 3 phases, 10 decorated classes.
2 assemblages lie inside no Phillips (1970) phase area and take the nearest, a geographic
rule; they are EXCLUDED from the column the conclusion rests on, leaving 26
assemblages inside a phase area.

Agreement is the adjusted Rand index: 0 in expectation for unrelated partitions, 1 for
identical ones. Weight 0 is the site map alone; the last row is the pottery alone.

| composition transform | weight on composition | ARI, all 28 | ARI, the 26 inside an area | k-means seed range |
|---|---|---|---|---|
| raw | **map alone** | 0.509 | **0.480** | 0.480-0.480 |
| raw | 0.1 | 0.509 | **0.480** | 0.480-0.480 |
| raw | 0.25 | 0.607 | **0.580** | 0.580-0.580 |
| raw | 0.5 | 0.492 | **0.508** | 0.508-0.508 |
| raw | 1 | 0.492 | **0.508** | 0.508-0.508 |
| raw | 2 | 0.115 | **0.167** | 0.167-0.167 |
| raw | 4 | 0.115 | **0.167** | 0.167-0.167 |
| raw | 8 | 0.115 | **0.167** | 0.167-0.167 |
| raw | 32 | 0.077 | **0.128** | 0.128-0.128 |
| raw | pottery alone | 0.077 | **0.128** | 0.128-0.128 |
| clr | **map alone** | 0.509 | **0.480** | 0.480-0.480 |
| clr | 0.1 | 0.509 | **0.480** | 0.480-0.480 |
| clr | 0.25 | 0.607 | **0.580** | 0.580-0.580 |
| clr | 0.5 | 0.509 | **0.522** | 0.522-0.522 |
| clr | 1 | 0.492 | **0.508** | 0.508-0.508 |
| clr | 2 | 0.115 | **0.167** | 0.167-0.167 |
| clr | 4 | 0.115 | **0.167** | 0.167-0.167 |
| clr | 8 | 0.115 | **0.167** | 0.167-0.167 |
| clr | 32 | 0.077 | **0.128** | 0.128-0.128 |
| clr | pottery alone | 0.077 | **0.128** | 0.128-0.128 |
| chisq | **map alone** | 0.509 | **0.480** | 0.480-0.480 |
| chisq | 0.1 | 0.509 | **0.480** | 0.480-0.480 |
| chisq | 0.25 | 0.756 | **0.740** | 0.740-0.740 |
| chisq | 0.5 | 0.756 | **0.740** | 0.740-0.740 |
| chisq | 1 | 0.756 | **0.740** | 0.740-0.740 |
| chisq | 2 | 0.485 | **0.539** | 0.539-0.539 |
| chisq | 4 | 0.485 | **0.539** | 0.539-0.539 |
| chisq | 8 | 0.389 | **0.435** | 0.435-0.435 |
| chisq | 32 | 0.326 | **0.358** | 0.358-0.358 |
| chisq | pottery alone | 0.326 | **0.358** | 0.358-0.358 |

## Does the verdict depend on the algorithm?

| transform | algorithm | map alone | pottery alone |
|---|---|---|---|
| raw | kmeans | 0.480 | 0.128 |
| raw | ward | 0.740 | 0.128 |
| raw | average | 0.740 | 0.128 |
| clr | kmeans | 0.480 | 0.128 |
| clr | ward | 0.740 | 0.131 |
| clr | average | 0.740 | 0.131 |
| chisq | kmeans | 0.480 | 0.358 |
| chisq | ward | 0.740 | 0.358 |
| chisq | average | 0.740 | 0.084 |

## How much of the log-ratio answer is the zero convention?

105 of 280 cells (38 percent) are zero, so the log-ratio transform
cannot be computed without deciding what a zero is worth. Map plus log-ratio composition
at half weight, on the assemblages inside a phase area, under five conventions:

| zero convention | ARI |
|---|---|
| pseudocount 0.5 (the default here) | 0.522 |
| pseudocount 5.0 | 0.508 |
| floor 2.0x min positive | 0.740 |
| floor 0.5x min positive | 0.740 |
| floor 0.1x min positive | 0.863 |

A convention that amplifies absence harder gives a better match, up to a
perfect one. That is a property of the convention, not of the pottery, and it
is why the chi-square transform -- which needs no such choice -- is the one
the reading below uses.

## Which assemblages each clustering misplaces

Each cluster is given the phase most of its members carry; an assemblage is misplaced
when its own phase is another. All assemblages are clustered and vote; the count is over the
26 inside a phase area, the ones the index is scored on, and an assemblage placed by the
nearest area is listed in the last column and not counted. Chi-square composition; k-means
is the run from seed 0 (the seed range in the first table shows whether the start matters).

| clustered on | method | adjusted Rand index (inside an area) | misplaced, of those scored | misplaced, placed by nearest area (not scored) |
|---|---|---|---|---|
| site map alone | kmeans | 0.480 | 6: Beck (Kent, with Walls); Belle_Meade (Kent, with Walls); Big_Eddy (Parkin, with Kent); Castile_Landing (Parkin, with Kent); Commerce (Kent, with Walls); Hollywood (Kent, with Walls) | none |
| site map alone | ward | 0.740 | 2: Big_Eddy (Parkin, with Kent); Castile_Landing (Parkin, with Kent) | none |
| site map alone | average | 0.740 | 2: Big_Eddy (Parkin, with Kent); Castile_Landing (Parkin, with Kent) | none |
| map plus composition, weight 0.5 | kmeans | 0.740 | 2: Big_Eddy (Parkin, with Kent); Castile_Landing (Parkin, with Kent) | none |
| map plus composition, weight 0.5 | ward | 0.740 | 2: Big_Eddy (Parkin, with Kent); Castile_Landing (Parkin, with Kent) | none |
| map plus composition, weight 0.5 | average | 0.740 | 2: Big_Eddy (Parkin, with Kent); Castile_Landing (Parkin, with Kent) | none |
| composition alone | kmeans | 0.358 | 6: Beck (Kent, with Parkin); Belle_Meade (Kent, with Parkin); Big_Eddy (Parkin, with Kent); Commerce (Kent, with Parkin); Mound_Place (Walls, with Parkin); Nickel (Kent, with Parkin) | Cramor_Place (Kent, with Parkin) |
| composition alone | ward | 0.358 | 6: Beck (Kent, with Parkin); Belle_Meade (Kent, with Parkin); Big_Eddy (Parkin, with Kent); Commerce (Kent, with Parkin); Mound_Place (Walls, with Parkin); Nickel (Kent, with Parkin) | Cramor_Place (Kent, with Parkin) |
| composition alone | average | 0.084 | 13: Barton_Ranch (Parkin, with Kent); Big_Eddy (Parkin, with Kent); Castile_Landing (Parkin, with Kent); Cummins (Parkin, with Kent); Fortune (Parkin, with Kent); Mound_Place (Walls, with Kent); Neeleys_Ferry (Parkin, with Kent); Parkin (Parkin, with Kent); Rose_Mound (Parkin, with Kent); Vernon_Paul (Parkin, with Kent); Walls (Walls, with Kent); Williamson (Parkin, with Kent); Woodlyn (Walls, with Kent) | none |

## Reading

Partitions carrying the phases' own group sizes with boundaries placed at random agree with the phases at a median ARI of 0.485 (400 draws), which is what
the group sizes manufacture on their own.

Under the chi-square transform, which carries no free parameter, k-means, on the 26 assemblages
assigned directly: the site map alone, 0.480; the best agreement, 0.740, at a composition weight of 0.25; the pottery alone, 0.358.

How to read the three values. A scheme that divides the map would be recovered by
the map alone, well above what the group sizes manufacture, and adding the pottery
would not help. A scheme that reads the pottery would be recovered better as the
pottery enters. This file states the values and draws no conclusion from them.

**What this does not show.** Agreement with a published scheme says how the scheme
was drawn, not how its units behaved: workers who sorted these same collections by
resemblance and by proximity are matched best by a clustering that uses both.

Figure written to fig14_phase_recovery.png and its siblings; full grid in phase_recovery.csv.
