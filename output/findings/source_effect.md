# Is the drift shortfall an artefact of merging data sources?

Basin phase set, 28 assemblages: 9 rows carrying the survey's counts alone (PFGData.xlsx) and 19 carrying Lipo's (2001) compilation, which holds the Phillips-Ford-Griffin counts plus his 1996-97 field collections. The two are tallies of the same material, so each assemblage uses one of them, the larger, and none is a sum (scripts/build_analysis_matrix.py).

## 1. Is source confounded with place?

| cluster | assemblages | survey alone | Lipo compilation | sherds, survey alone | sherds, Lipo compilation |
|---|---|---|---|---|---|
| 0 (Barton_Ranch, Cummins, Fortune...) | 10 | 4 | 6 | 1,247 | 5,308 |
| 1 (Clay_Hill, Davis, Grant...) | 5 | 2 | 3 | 274 | 491 |
| 2 (Irby, Lake_Cormorant, Mound_Place...) | 5 | 2 | 3 | 275 | 643 |
| 3 (Big_Eddy, Castile_Landing, Cramor_Place...) | 4 | 1 | 3 | 229 | 2,338 |
| 4 (Beck, Belle_Meade, Commerce...) | 4 | 0 | 4 | 0 | 3,296 |

## 2. Do the sources differ once place is held?

| cluster | F_ST, survey-alone rows against compilation rows | random splits of the same sizes, median (5th-95th) | position |
|---|---|---|---|
| 0 | 0.0009 | 0.0040 (0.0003-0.0200) | above 17 percent of 2,000 |
| 1 | 0.0001 | 0.0033 (0.0001-0.0118) | above 0 percent of 2,000 |
| 2 | 0.0033 | 0.0069 (0.0021-0.0296) | above 18 percent of 2,000 |
| 3 | fewer than two rows of one kind | | not testable |
| 4 | fewer than two rows of one kind | | not testable |

## 3. Does the excess over drift survive in a single source?

Between-cluster F_ST at k = 5, the same clusters throughout, against 200 calibrated drift realizations scored on the same subset.

| subset | assemblages | clusters represented | median sherds per assemblage | observed F_ST | drift median (95 percent) | shortfall |
|---|---|---|---|---|---|---|
| all 43 | 28 | 5 | 251 | 0.0307 | 0.0065 (0.0021-0.0210) | **4.8x** |
| rows carrying the survey's counts alone | 9 | 4 | 181 | 0.0998 | 0.0107 (0.0040-0.0388) | **9.3x** |
| rows carrying Lipo's compilation | 19 | 5 | 424 | 0.0261 | 0.0077 (0.0024-0.0215) | **3.4x** |
| rows carrying the survey's counts alone, shared clusters only | 9 | 4 | 181 | 0.0998 | 0.0107 (0.0040-0.0388) | **9.3x** |
| rows carrying Lipo's compilation, shared clusters only | 15 | 4 | 424 | 0.0295 | 0.0056 (0.0017-0.0217) | **5.2x** |

**Source and sample size are not separable here.** The rows carrying the survey's counts alone are also the small
ones, so a larger shortfall among them may be an analyst or collection-regime effect, or
overdispersion in small surface collections, or both. Either is a property of the record
rather than of past interaction, which is the distinction that matters for the residual.
