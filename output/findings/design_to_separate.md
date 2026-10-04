# What record would separate local innovation from restricted copying?

Produced by `analyses/97_design_to_separate.py`. 300 runs per setting on analysis 84's sixteen-setting grid, at the calibrated combination (N 2000, innovation 0.001, mixing 0.02), for each design. Each run made with **restricted copying** (copying across the phase lines cut to 0.03, one regional pool) or with **local innovation** (copying unrestricted, 0.2 of classes reordered by phase) is put through analysis 84's rejection ABC (5% of runs kept, six summaries), left out of its own reference table, and the posterior probability of a restriction at the phase lines is read off. The prior is 0.75.

Added sites: 53 late-period sites of the settlement compilation lie more than 1 km from every analyzed assemblage (analysis 80's selection); 43 of them fall in or nearest the Kent, Parkin or Walls areas and are used. Phases of the 71-site design: Kent 28, Parkin 29, Walls 14. Sherd counts for added sites are drawn from the observed counts (seed 97). The excluded collections are those of `output/findings/small_collections.md`, each at its own location, Phillips phase and decorated-sherd count.

| design | sites | classes | P(restriction), made with one: median (10th to 90th) | made without one: median (10th to 90th) | pairs ordered correctly | between-phase F_ST, medians (with / without) | boundary excess, medians (with / without) |
|---|---|---|---|---|---|---|---|
| the record as it is | 28 | 10 | 0.82 (0.70 to 0.94) | 0.58 (0.44 to 0.75) | 0.93 | 0.0062 / 0.0056 | +6.5 / +2.9 |
| the record as it is, sherds x4 | 28 | 10 | 0.84 (0.70 to 0.95) | 0.58 (0.44 to 0.74) | 0.94 | 0.0063 / 0.0057 | +6.4 / +2.5 |
| classes split in two | 28 | 20 | 0.91 (0.73 to 0.99) | 0.56 (0.46 to 0.73) | 0.97 | 0.0070 / 0.0049 | +10.6 / +3.2 |
| classes split in two, sherds x4 | 28 | 20 | 0.94 (0.78 to 1.00) | 0.54 (0.46 to 0.71) | 0.99 | 0.0071 / 0.0049 | +11.0 / +2.7 |
| classes split in four | 28 | 40 | 0.99 (0.85 to 1.00) | 0.53 (0.42 to 0.66) | 1.00 | 0.0082 / 0.0045 | +14.4 / +3.9 |
| classes split in four, sherds x4 | 28 | 40 | 1.00 (0.87 to 1.00) | 0.52 (0.41 to 0.66) | 1.00 | 0.0081 / 0.0043 | +15.2 / +3.5 |
| 21 sites added | 49 | 10 | 0.83 (0.71 to 0.90) | 0.56 (0.45 to 0.77) | 0.94 | 0.0039 / 0.0027 | +2.7 / +1.3 |
| 21 sites added, sherds x4 | 49 | 10 | 0.85 (0.70 to 0.93) | 0.57 (0.44 to 0.75) | 0.94 | 0.0039 / 0.0024 | +3.2 / +1.4 |
| 43 sites added | 71 | 10 | 0.81 (0.68 to 0.88) | 0.60 (0.47 to 0.76) | 0.90 | 0.0027 / 0.0024 | +2.5 / +1.9 |
| 43 sites added, sherds x4 | 71 | 10 | 0.83 (0.70 to 0.91) | 0.58 (0.46 to 0.75) | 0.92 | 0.0027 / 0.0025 | +2.7 / +2.1 |
| 43 sites added, classes split in four | 71 | 40 | 0.92 (0.80 to 0.99) | 0.53 (0.43 to 0.67) | 0.99 | 0.0032 / 0.0016 | +5.1 / +2.1 |
| 43 sites added, classes split in four, sherds x4 | 71 | 40 | 0.95 (0.83 to 1.00) | 0.51 (0.42 to 0.65) | 0.99 | 0.0032 / 0.0016 | +5.8 / +2.1 |
| the 7 excluded collections with at least 50 decorated sherds added | 35 | 10 | 0.83 (0.70 to 0.95) | 0.59 (0.47 to 0.73) | 0.94 | 0.0059 / 0.0049 | +5.4 / +2.3 |
| the 15 excluded collections with at least 25 decorated sherds added | 43 | 10 | 0.79 (0.67 to 0.88) | 0.62 (0.50 to 0.76) | 0.89 | 0.0046 / 0.0048 | +3.3 / +1.8 |
| the 24 excluded collections with at least 10 decorated sherds added | 52 | 10 | 0.78 (0.69 to 0.88) | 0.61 (0.47 to 0.76) | 0.89 | 0.0038 / 0.0032 | +3.2 / +2.3 |
| the 24 excluded collections with at least 10 added, classes split in four | 52 | 40 | 0.92 (0.75 to 0.99) | 0.54 (0.45 to 0.74) | 0.98 | 0.0047 / 0.0023 | +6.8 / +3.1 |

## Reading

"Pairs ordered correctly" is the share of (restricted, local-innovation) pairs of records in which the restricted one gets the higher posterior probability of a restriction: 0.5 means the procedure cannot tell the two explanations apart, 1 that it always can. It is a description of the simulated records, not a test (rule 18).

On the record as it is, the share is 0.93, and the posterior is 0.82 for records made with a restriction against 0.58 for records made without one (prior 0.75). The highest share among the designs is 1.00, for "classes split in four, sherds x4".

The model is not recalibrated to the finer classes or the larger site set, so compare the designs with each other and not with the record. Finer classes are a geometric split of the ten classes' frequencies; real attribute classes need not behave that way. The added sites are where late-period sites are recorded, which is not where collections exist.
