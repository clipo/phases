# What Mainfort's rim-attribute file holds, as counts

Written by `scripts/build_mainfort_rim_counts.py`. The file (`data/raw/mainfort_miss_rims.xls`, sent by R. C. Mainfort, Jr.) is the data behind Mainfort (2003b), "Late Period Ceramic Rim Attribute Variation in the Central Mississippi Valley", *Southeastern Archaeology* 22(1):33-46. That paper's Table 1 prints the counts of six attributes for the same 55 sites and 7,195 rims; it is transcribed in `data/raw/mainfort2003_rim_table1_counts.csv` and compared with the file here.

## The columns

| column(s) | what it is | how the values relate |
|---|---|---|
| `total` | N, the number of rims analyzed | equals the published N at all 55 sites |
| `beadrim`, `noderim`, `widenode`, `notchtop`, `notchext` | beaded applique strip; band of nodes; widely spaced nodes; notched top of the lip; notched exterior | each is the percentage of the N rims with that attribute. Separate yes-or-no observations: they do not sum to 100 with one another |
| `sharpbev`, `modbev`, `slitebev`, `nobev` | degree of interior bevel | parts of one whole: sum 77 to 110 (100 within rounding at 52 of 55 sites). The first three together are the paper's "Beveled" |
| `flatlip`, `roundlip` | lip form | parts of one whole: sum 98 to 101 |
| `nf`, `bell` | Neeley's Ferry Plain and Bell Plain | NOT parts of one whole: sum 20 to 102, median 89. The paste of the rim sherds (author), so percentages of N; what the remaining rims are is not established. Not in the paper |
| `latitude`, `longitude` | rough locations | some are wrong (Davis is at latitude 24.86); the project's own coordinates are used for matched sites |

## The check against the published counts

For each of the six published attributes, the published count divided by N is compared with the file's percentage. They agree within 1.5 points except where listed, so the file's percentages are percentages of N.

- Bishop, bevel degrees: the file's four bevel percentages sum to 77, not 100; the published beveled count (20 of 32) is used and divided in the file's proportions
- Denton Mounds, bevel degrees: the file's four bevel percentages sum to 110, not 100; the published beveled count (49 of 145) is used and divided in the file's proportions
- Gant, beveled: published 9 of 84 (10.7 percent) against 9 percent in the file; the published count is used
- Graves Lake, beaded: printed count 12 of 203 (5.9 percent) against 9 percent in the file (which would be 18); the printed count is used
- Graves Lake, bevel degrees: the file's four bevel percentages sum to 110, not 100; the published beveled count (91 of 203) is used and divided in the file's proportions
- Pecan Point, beaded: printed count 10 of 30 (33.3 percent) against 10 percent in the file (which would be 3); the printed count is used
- Walls, beveled: published 14 of 107 (13.1 percent) against 16 percent in the file; the published count is used

The paper's Table 2 (percentages) lists Parkin with N = 111, where Table 1 and the file have 67; the file's Parkin percentages are those of 67 rims.

## The table as counts

`data/processed/mainfort_rim_counts.csv`, 55 rows. The six attributes are the published counts as printed; bevel degrees divide the published beveled count in the file's proportions; lip form and paste are the file's percentages times N and are not in the paper.

| site | N | beaded | noded | wide nodes | notched top | notched exterior | bevel: sharp / moderate / slight / none | lip: flat / round | Neeley's Ferry / Bell / other |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| 3CT245 | 57 | 2 | 1 | 1 | 0 | 15 | 5 / 5 / 5 / 42 | 11 / 46 | 22 / 28 / 7 |
| 40LA17 | 27 | 7 | 1 | 0 | 0 | 6 | 5 / 10 / 1 / 11 | 8 / 19 | 13 / 14 / 0 |
| Armorel | 141 | 29 | 4 | 0 | 6 | 24 | 26 / 9 / 1 / 105 | 55 / 86 | 34 / 85 / 22 |
| Barton Ranch | 102 | 0 | 0 | 0 | 1 | 6 | 1 / 0 / 1 / 100 | 35 / 67 | 66 / 16 / 20 |
| Beck | 183 | 4 | 4 | 0 | 2 | 22 | 16 / 20 / 6 / 141 | 79 / 104 | 70 / 73 / 40 |
| Bell | 35 | 2 | 0 | 0 | 0 | 0 | 4 / 2 / 2 / 27 | 6 / 29 | 6 / 18 / 11 |
| Belle Meade | 143 | 4 | 3 | 1 | 0 | 29 | 12 / 21 / 11 / 99 | 74 / 69 | 56 / 34 / 53 |
| Berry | 473 | 186 | 0 | 0 | 4 | 77 | 152 / 152 / 76 / 93 | 151 / 322 | 166 / 218 / 89 |
| Big Eddy | 88 | 9 | 0 | 0 | 0 | 1 | 8 / 4 / 2 / 74 | 20 / 68 | 40 / 23 / 25 |
| Bishop | 32 | 9 | 1 | 0 | 0 | 1 | 5 / 10 / 5 / 12 | 13 / 19 | 23 / 8 / 1 |
| Bradley (3CT7) | 337 | 3 | 1 | 2 | 13 | 77 | 8 / 7 / 0 / 322 | 135 / 202 | 148 / 189 / 0 |
| Cagle Lake | 273 | 54 | 10 | 0 | 2 | 31 | 37 / 43 / 24 / 169 | 71 / 202 | 134 / 126 / 13 |
| Campbell | 115 | 43 | 0 | 0 | 0 | 37 | 30 / 43 / 19 / 23 | 40 / 75 | 49 / 58 / 8 |
| Carson Lake | 156 | 7 | 3 | 0 | 2 | 41 | 28 / 27 / 11 / 90 | 48 / 108 | 105 / 48 / 3 |
| Castile Landing | 72 | 0 | 1 | 2 | 1 | 0 | 0 / 1 / 0 / 71 | 12 / 60 | 45 / 15 / 12 |
| Cheatam | 18 | 0 | 1 | 0 | 1 | 1 | 3 / 1 / 0 / 14 | 9 / 9 | 4 / 14 / 0 |
| Chucalissa | 90 | 8 | 1 | 0 | 2 | 11 | 3 / 7 / 0 / 80 | 36 / 54 | 32 / 56 / 2 |
| Clay Hill | 159 | 0 | 0 | 0 | 4 | 10 | 5 / 9 / 5 / 140 | 59 / 100 | 92 / 11 / 56 |
| Davis | 96 | 0 | 0 | 0 | 1 | 2 | 0 / 1 / 1 / 94 | 31 / 65 | 54 / 15 / 27 |
| Denton Mounds | 145 | 13 | 2 | 0 | 2 | 13 | 13 / 26 / 10 / 96 | 30 / 115 | 100 / 22 / 23 |
| Dry Arm | 30 | 4 | 0 | 0 | 0 | 3 | 7 / 5 / 1 / 17 | 6 / 24 | 21 / 7 / 2 |
| Fortune | 84 | 0 | 0 | 0 | 0 | 0 | 0 / 0 / 0 / 84 | 16 / 68 | 8 / 8 / 68 |
| Fullen | 111 | 7 | 3 | 0 | 3 | 19 | 12 / 9 / 13 / 77 | 60 / 51 | 48 / 60 / 3 |
| Gant | 84 | 5 | 0 | 0 | 2 | 15 | 1 / 6 / 2 / 75 | 56 / 28 | 56 / 18 / 10 |
| Graves Lake | 203 | 12 | 3 | 1 | 5 | 31 | 20 / 36 / 35 / 112 | 87 / 116 | 112 / 87 / 4 |
| Hatchie | 140 | 12 | 9 | 5 | 3 | 28 | 22 / 17 / 8 / 93 | 41 / 99 | 49 / 85 / 6 |
| Irby | 277 | 3 | 2 | 2 | 3 | 57 | 17 / 17 / 14 / 229 | 101 / 176 | 47 / 224 / 6 |
| Jeter | 61 | 2 | 0 | 0 | 0 | 10 | 3 / 2 / 0 / 56 | 27 / 34 | 26 / 30 / 5 |
| Jones Bayou | 96 | 14 | 6 | 1 | 1 | 9 | 16 / 11 / 9 / 60 | 59 / 37 | 45 / 49 / 2 |
| Kent | 153 | 1 | 0 | 1 | 3 | 4 | 3 / 8 / 3 / 139 | 54 / 99 | 73 / 35 / 45 |
| Knappenberger | 219 | 36 | 9 | 1 | 8 | 40 | 32 / 32 / 17 / 138 | 64 / 155 | 88 / 107 / 24 |
| Middle Nodena | 117 | 14 | 0 | 2 | 1 | 41 | 22 / 30 / 3 / 62 | 40 / 77 | 50 / 62 / 5 |
| Mound Place | 48 | 2 | 0 | 0 | 1 | 0 | 1 / 3 / 0 / 44 | 17 / 31 | 17 / 31 / 0 |
| Neeleys Ferry | 108 | 1 | 1 | 0 | 2 | 1 | 0 / 0 / 1 / 107 | 23 / 85 | 54 / 14 / 40 |
| Neely (3CT40) | 253 | 15 | 3 | 3 | 5 | 82 | 8 / 26 / 10 / 209 | 129 / 124 | 76 / 177 / 0 |
| Notgrass | 190 | 57 | 7 | 0 | 3 | 23 | 2 / 4 / 7 / 177 | 65 / 125 | 57 / 106 / 27 |
| Otto Sharpe | 402 | 93 | 3 | 7 | 0 | 21 | 48 / 97 / 92 / 165 | 84 / 318 | 189 / 133 / 80 |
| Parkin | 67 | 0 | 0 | 1 | 0 | 4 | 0 / 1 / 4 / 62 | 30 / 37 | 37 / 12 / 18 |
| Pecan Point | 30 | 10 | 1 | 0 | 0 | 6 | 7 / 1 / 0 / 22 | 14 / 16 | 13 / 17 / 0 |
| Porter | 111 | 12 | 11 | 2 | 3 | 8 | 5 / 2 / 9 / 95 | 83 / 28 | 83 / 28 / 0 |
| Rast | 102 | 6 | 0 | 1 | 2 | 9 | 11 / 14 / 5 / 72 | 41 / 61 | 31 / 71 / 0 |
| Richardson's Landing | 115 | 9 | 1 | 0 | 3 | 14 | 8 / 8 / 11 / 88 | 66 / 49 | 60 / 52 / 3 |
| Rose Mound | 76 | 0 | 1 | 0 | 0 | 0 | 1 / 1 / 0 / 74 | 24 / 52 | 57 / 0 / 19 |
| Shawnee Village | 133 | 2 | 2 | 1 | 2 | 21 | 0 / 3 / 3 / 127 | 57 / 76 | 64 / 51 / 18 |
| Soudan | 154 | 0 | 0 | 0 | 9 | 1 | 0 / 1 / 1 / 152 | 46 / 108 | 108 / 9 / 37 |
| Starkley | 108 | 2 | 1 | 1 | 2 | 7 | 13 / 6 / 2 / 87 | 26 / 82 | 42 / 26 / 40 |
| Sweat | 132 | 29 | 2 | 2 | 0 | 12 | 29 / 36 / 8 / 59 | 70 / 62 | 69 / 61 / 2 |
| Upper Nodena | 103 | 26 | 6 | 1 | 1 | 27 | 8 / 15 / 6 / 74 | 49 / 54 | 14 / 85 / 4 |
| Vernon Paul | 99 | 0 | 0 | 1 | 1 | 5 | 2 / 0 / 1 / 96 | 37 / 62 | 66 / 10 / 23 |
| Walls | 107 | 1 | 0 | 1 | 1 | 8 | 7 / 5 / 2 / 93 | 45 / 62 | 18 / 66 / 23 |
| Wilder | 141 | 11 | 4 | 1 | 2 | 26 | 23 / 20 / 8 / 90 | 70 / 71 | 58 / 82 / 1 |
| Wildy | 179 | 3 | 0 | 1 | 8 | 22 | 2 / 0 / 1 / 176 | 90 / 89 | 134 / 39 / 6 |
| Williamson | 68 | 0 | 0 | 0 | 0 | 0 | 1 / 0 / 1 / 66 | 19 / 49 | 24 / 9 / 35 |
| Woodlyn | 27 | 1 | 0 | 0 | 1 | 4 | 3 / 6 / 1 / 17 | 6 / 21 | 8 / 19 / 0 |
| Young | 125 | 0 | 2 | 1 | 1 | 20 | 21 / 14 / 3 / 87 | 48 / 77 | 21 / 88 / 16 |

## What the paste percentages are percentages of

The author reads `nf` and `bell` as the paste of the rim sherds on which the attributes were recorded. Two checks agree that the base is the N rims and not the whole sherd assemblage.

- Each value is a whole number of rims out of N at 53 of 55 sites for `nf` and 54 for `bell`.
- At the 32 sites also in Mainfort (2003a: Table 1), Bell Plain's share of all sherds there differs from this file's `bell` by a mean of 11 points. Parkin is 18 percent Bell here and 1 percent there. Rims and body sherds are different samples of the pottery, so the two need not agree.

## To confirm with Mainfort

1. What the rims in neither paste are. `nf` and `bell` are the paste of the rim sherds (see above), and the remainder is large where decorated jars are common (Fortune 20 percent in total).
2. Whether lip form and bevel degree were recorded on all N rims. A few percentages cannot be a whole number of rims out of N (for example Bishop's lip values), which suggests a slightly smaller base at some sites.
3. The rows where the file and the published counts disagree, listed above.
