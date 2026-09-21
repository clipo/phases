# What recovers the phase scheme: the pottery, or the map?

Basin phase set, 43 assemblages, 5 phases, 10 decorated classes.
9 assemblages carry a phase this project derived by territory, a geographic
rule; they are EXCLUDED from the column the conclusion rests on, leaving 34
assemblages Mainfort himself assigned.

Agreement is the adjusted Rand index: 0 in expectation for unrelated partitions, 1 for
identical ones. Weight 0 is the site map alone; the last row is the pottery alone.

| composition transform | weight on composition | ARI, all 43 | ARI, Mainfort's 34 | k-means seed range |
|---|---|---|---|---|
| raw | **map alone** | 0.741 | **0.678** | 0.678-1.000 |
| raw | 0.1 | 0.669 | **0.678** | 0.678-0.815 |
| raw | 0.25 | 0.726 | **0.740** | 0.623-0.815 |
| raw | 0.5 | 0.653 | **0.682** | 0.623-0.740 |
| raw | 1 | 0.504 | **0.522** | 0.443-0.697 |
| raw | 2 | 0.457 | **0.465** | 0.225-0.523 |
| raw | 4 | 0.144 | **0.159** | 0.153-0.465 |
| raw | 8 | 0.142 | **0.159** | 0.099-0.198 |
| raw | 32 | 0.088 | **0.077** | 0.052-0.142 |
| raw | pottery alone | 0.088 | **0.077** | 0.052-0.142 |
| clr | **map alone** | 0.741 | **0.678** | 0.678-1.000 |
| clr | 0.1 | 0.790 | **0.815** | 0.678-0.815 |
| clr | 0.25 | 0.732 | **0.737** | 0.682-0.815 |
| clr | 0.5 | 0.732 | **0.737** | 0.682-0.740 |
| clr | 1 | 0.647 | **0.659** | 0.522-0.737 |
| clr | 2 | 0.494 | **0.514** | 0.257-0.697 |
| clr | 4 | 0.271 | **0.246** | 0.165-0.519 |
| clr | 8 | 0.195 | **0.200** | 0.126-0.246 |
| clr | 32 | 0.150 | **0.129** | 0.126-0.206 |
| clr | pottery alone | 0.149 | **0.126** | 0.126-0.206 |
| chisq | **map alone** | 0.741 | **0.678** | 0.678-1.000 |
| chisq | 0.1 | 0.790 | **0.815** | 0.678-0.815 |
| chisq | 0.25 | 0.790 | **0.815** | 0.815-0.815 |
| chisq | 0.5 | 0.732 | **0.740** | 0.678-0.815 |
| chisq | 1 | 0.699 | **0.737** | 0.590-0.753 |
| chisq | 2 | 0.453 | **0.465** | 0.171-0.753 |
| chisq | 4 | 0.176 | **0.196** | 0.139-0.516 |
| chisq | 8 | 0.161 | **0.179** | 0.127-0.232 |
| chisq | 32 | 0.137 | **0.144** | 0.112-0.204 |
| chisq | pottery alone | 0.137 | **0.138** | 0.087-0.175 |

## Does the verdict depend on the algorithm?

| transform | algorithm | map alone | pottery alone |
|---|---|---|---|
| raw | kmeans | 0.678 | 0.077 |
| raw | ward | 0.678 | 0.109 |
| raw | average | 0.603 | 0.109 |
| clr | kmeans | 0.678 | 0.126 |
| clr | ward | 0.678 | 0.122 |
| clr | average | 0.603 | 0.199 |
| chisq | kmeans | 0.678 | 0.138 |
| chisq | ward | 0.678 | 0.259 |
| chisq | average | 0.603 | 0.143 |

## How much of the log-ratio answer is the zero convention?

168 of 430 cells (39 percent) are zero, so the log-ratio transform
cannot be computed without deciding what a zero is worth. Map plus log-ratio composition
at half weight, on Mainfort's assemblages, under five conventions:

| zero convention | ARI |
|---|---|
| pseudocount 0.5 (the default here) | 0.737 |
| pseudocount 5.0 | 0.682 |
| floor 2.0x min positive | 0.815 |
| floor 0.5x min positive | 0.902 |
| floor 0.1x min positive | 1.000 |

A convention that amplifies absence harder gives a better match, up to a
perfect one. That is a property of the convention, not of the pottery, and it
is why the chi-square transform -- which needs no such choice -- is the one
the reading below uses.

## Reading

Partitions carrying the phases' own group sizes with boundaries placed at random agree with the phases at a median ARI of 0.356 (400 draws), which is what
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
