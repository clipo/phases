# Is the drift shortfall an artefact of merging data sources?

Basin phase set, 28 assemblages: 9 rows carrying the survey's counts alone (PFGData.xlsx) and 19 carrying Lipo's (2001) compilation, which holds the Phillips-Ford-Griffin counts plus his 1996-97 field collections. The two are tallies of the same material, so each assemblage uses one of them, the larger, and none is a sum (scripts/build_analysis_matrix.py).

## 1. Is source confounded with place?

| cluster | assemblages | survey alone | Lipo compilation | sherds, survey alone | sherds, Lipo compilation |
|---|---|---|---|---|---|
| 0 (Barton_Ranch, Fortune, Holden_Lake...) | 9 | 3 | 6 | 986 | 5,308 |
| 1 (Irby, Lake_Cormorant, Mound_Place...) | 5 | 2 | 3 | 275 | 643 |
| 2 (Cummins...) | 1 | 1 | 0 | 261 | 0 |
| 3 (Big_Eddy, Castile_Landing, Clay_Hill...) | 5 | 2 | 3 | 410 | 2,338 |
| 4 (Beck, Belle_Meade, Commerce...) | 8 | 1 | 7 | 93 | 3,787 |

## 2. Do the sources differ once place is held?

Stated as a posterior probability (rule 18, 2026-09-23): each assemblage's class proportions are drawn from Dirichlet(counts + 1/2), and for each draw the F_ST between the two kinds of row is compared with the F_ST of a random split of the same cluster into groups of the same sizes.

| cluster | F_ST, survey-alone rows against compilation rows, posterior median [95%] | random splits, posterior median | P(source split differs more than a random split) |
|---|---|---|---|
| 0 | 0.0003 [0.0001, 0.0009] | 0.0042 | 0.04 |
| 1 | 0.0033 [0.0011, 0.0076] | 0.0075 | 0.15 |
| 2 | fewer than two rows of one kind | | not testable |
| 3 | 0.0320 [0.0232, 0.0418] | 0.0110 | 0.89 |
| 4 | fewer than two rows of one kind | | not testable |

## 3. Does the excess over drift survive in a single source?

Between-cluster F_ST at k = 5, the same clusters throughout, against 200 calibrated drift realizations scored on the same subset.

| subset | assemblages | clusters represented | median sherds per assemblage | observed F_ST | drift median (95 percent) | shortfall |
|---|---|---|---|---|---|---|
| all 28 | 28 | 5 | 251 | 0.0206 | 0.0044 (0.0011-0.0132) | **4.6x** |
| rows carrying the survey's counts alone | 9 | 5 | 181 | 0.1138 | 0.0075 (0.0021-0.0221) | **15.2x** |
| rows carrying Lipo's compilation | 19 | 4 | 424 | 0.0163 | 0.0046 (0.0010-0.0153) | **3.5x** |
| rows carrying the survey's counts alone, shared clusters only | 8 | 4 | 174 | 0.1014 | 0.0074 (0.0019-0.0241) | **13.8x** |
| rows carrying Lipo's compilation, shared clusters only | 19 | 4 | 424 | 0.0163 | 0.0046 (0.0010-0.0153) | **3.5x** |

**Source and sample size are not separable here.** The rows carrying the survey's counts alone are also the small
ones, so a larger shortfall among them may be an analyst or collection-regime effect, or
overdispersion in small surface collections, or both. Either is a property of the record
rather than of past interaction, which is the distinction that matters for the residual.
