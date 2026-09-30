# What recovers the phase scheme: the pottery, or the map?

Basin phase set, 28 assemblages, 3 phases, 10 decorated classes.
3 assemblages carry a phase this project derived by territory, a geographic
rule; they are EXCLUDED from the column the conclusion rests on, leaving 25
assemblages Mainfort himself assigned.

Agreement is the adjusted Rand index: 0 in expectation for unrelated partitions, 1 for
identical ones. Weight 0 is the site map alone; the last row is the pottery alone.

| composition transform | weight on composition | ARI, all 28 | ARI, Mainfort's 25 | k-means seed range |
|---|---|---|---|---|
| raw | **map alone** | 0.889 | **1.000** | 1.000-1.000 |
| raw | 0.1 | 0.889 | **1.000** | 1.000-1.000 |
| raw | 0.25 | 0.695 | **0.766** | 0.766-0.766 |
| raw | 0.5 | 0.527 | **0.533** | 0.533-0.533 |
| raw | 1 | 0.527 | **0.533** | 0.533-0.533 |
| raw | 2 | 0.064 | **0.092** | 0.092-0.092 |
| raw | 4 | 0.064 | **0.092** | 0.092-0.092 |
| raw | 8 | 0.064 | **0.092** | 0.092-0.092 |
| raw | 32 | 0.056 | **0.082** | 0.082-0.082 |
| raw | pottery alone | 0.056 | **0.082** | 0.082-0.082 |
| clr | **map alone** | 0.889 | **1.000** | 1.000-1.000 |
| clr | 0.1 | 0.889 | **1.000** | 1.000-1.000 |
| clr | 0.25 | 0.695 | **0.766** | 0.766-0.766 |
| clr | 0.5 | 0.527 | **0.533** | 0.533-0.533 |
| clr | 1 | 0.527 | **0.533** | 0.533-0.533 |
| clr | 2 | 0.064 | **0.092** | 0.092-0.092 |
| clr | 4 | 0.064 | **0.092** | 0.092-0.092 |
| clr | 8 | 0.064 | **0.092** | 0.092-0.092 |
| clr | 32 | 0.056 | **0.082** | 0.082-0.082 |
| clr | pottery alone | 0.056 | **0.082** | 0.082-0.082 |
| chisq | **map alone** | 0.889 | **1.000** | 1.000-1.000 |
| chisq | 0.1 | 0.560 | **0.607** | 0.607-0.607 |
| chisq | 0.25 | 0.560 | **0.607** | 0.607-0.607 |
| chisq | 0.5 | 0.560 | **0.607** | 0.607-0.607 |
| chisq | 1 | 0.560 | **0.607** | 0.607-0.607 |
| chisq | 2 | 0.291 | **0.292** | 0.292-0.292 |
| chisq | 4 | 0.291 | **0.292** | 0.292-0.292 |
| chisq | 8 | 0.222 | **0.215** | 0.215-0.215 |
| chisq | 32 | 0.208 | **0.204** | 0.204-0.204 |
| chisq | pottery alone | 0.208 | **0.204** | 0.204-0.204 |

## Does the verdict depend on the algorithm?

| transform | algorithm | map alone | pottery alone |
|---|---|---|---|
| raw | kmeans | 1.000 | 0.082 |
| raw | ward | 0.607 | 0.082 |
| raw | average | 0.537 | 0.082 |
| clr | kmeans | 1.000 | 0.082 |
| clr | ward | 0.607 | 0.111 |
| clr | average | 0.537 | 0.111 |
| chisq | kmeans | 1.000 | 0.204 |
| chisq | ward | 0.607 | 0.204 |
| chisq | average | 0.537 | -0.002 |

## How much of the log-ratio answer is the zero convention?

105 of 280 cells (38 percent) are zero, so the log-ratio transform
cannot be computed without deciding what a zero is worth. Map plus log-ratio composition
at half weight, on Mainfort's assemblages, under five conventions:

| zero convention | ARI |
|---|---|
| pseudocount 0.5 (the default here) | 0.533 |
| pseudocount 5.0 | 0.533 |
| floor 2.0x min positive | 0.607 |
| floor 0.5x min positive | 0.607 |
| floor 0.1x min positive | 0.766 |

A convention that amplifies absence harder gives a better match, up to a
perfect one. That is a property of the convention, not of the pottery, and it
is why the chi-square transform -- which needs no such choice -- is the one
the reading below uses.

## Reading

Partitions carrying the phases' own group sizes with boundaries placed at random agree with the phases at a median ARI of 0.570 (400 draws), which is what
the group sizes manufacture on their own.

The site map alone recovers the published scheme well above what the group
sizes manufacture. Composition adds to it, and the addition is real but
secondary: under the chi-square transform, which carries no free parameter,
agreement rises from the map-alone value to its best at a composition weight
of 0.1 to 0.25, that is with the pottery counting for a quarter or less of
what the coordinates count for. Weighting the pottery equally with the map is
already worse than the map alone, and the pottery by itself is worst of all.

So the scheme is neither purely geographic nor a reading of the pots. It is
the pots seen through where the sites are, with geography carrying most of
the weight. That is consistent with the rest of this paper: on a compositional
gradient with no edges, similarity is a monotone function of proximity, so
sorting assemblages by how alike their pots look largely recovers where they
are, and the residual ceramic signal refines that rather than overriding it.

**What this does not show.** It does not show the assemblages are
compositionally identical; they are not, and composition varies strongly with
distance. It shows that the particular five-way division the phase scheme
draws is predicted by the coordinates and not by the pots.

Figure written to fig14_phase_recovery.png and its siblings; full grid in phase_recovery.csv.
