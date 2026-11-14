# Composition by collection: the 1940 PFG surface collections against Lipo's later additions

Produced by `analyses/per_type/composition_by_collection.py`. Of the 28 basin assemblages, 18 are the 1940 Phillips-Ford-Griffin collections alone; 6 are those plus a large-count addition from Lipo's 1996 to 1997 collections (Belle_Meade +1455, Beck +1154, Nickel +1120, Rose_Mound +738, Cramor_Place +286, Castile_Landing +220); 3 carry a small one (Barton_Ranch +9, Hollywood +25, Neeleys_Ferry +4); Holden_Lake is Lipo's alone. Large and small are count-size labels at 100 sherds, not a judgment of compositional weight: Hollywood's 25 added sherds are 24 Painted, the largest change in any site's Painted share among the nine. Decomposition checked on this run: one cpl-pfg id per assemblage, the matrix row equal to that id's row, disjoint collection letters within and across assemblages, every addition non-negative (the addition is defined as the matrix row minus the PFG part, so that identity is not a check). Descriptive readings under rule 18. Intervals in section 2 are 5th to 95th percentiles of the contrast under independent multinomial sampling of each part with a symmetric Dirichlet(0.5) prior over the ten classes (20,000 draws, seed 2026101960); they describe the latent class proportions of each collection under that model and nothing about why the two collections differ (source attribution, vessel clustering, which parts of a site were walked, recovery or classification differences, or the fifty years between the visits are all outside them). Prior (rule 20): it leans toward equal class shares (1/K) in each part and moves a small collection further than a large one; its effect on the contrast between the parts is not a shrinkage toward zero, since pulling one part toward 1/K can attenuate, amplify or reverse a contrast, and the per-site alpha columns show all three (Neeleys Ferry's Parkin_Punctated contrast changes sign across the alphas); prior predictive share per class median 0.052, 90 percent interval [0.0004, 0.365] (a distribution of latent shares, not of counts at these collection sizes); sensitivity under alpha 0.1 and 1.0 is shown beside every per-site contrast; recovery on constructed counts, including ten classes at 25 against 200 sherds, is in the tests.

## 1. Share against collection size

Spearman correlation of each class's site share with the site's decorated-sherd total. The last two columns hold the site set fixed (the 27 sites whose PFG part has at least 75 sherds) and change only the collection: the current counts against the PFG parts alone. An undefined entry says why.

| class | all 28 rows | 18 rows, no addition | 21 rows, no large addition | 27 matched sites, current counts | 27 matched sites, PFG parts |
|---|---|---|---|---|---|
| Parkin_Punctated | +0.42 (28) | +0.48 (18) | +0.48 (21) | +0.46 (27) | +0.49 (27) |
| Barton/Kent/MPI | +0.26 (28) | +0.03 (18) | +0.14 (21) | +0.26 (27) | -0.00 (27) |
| Painted | -0.57 (28) | -0.56 (18) | -0.60 (21) | -0.58 (27) | -0.51 (27) |
| Fortune_Noded | +0.56 (15) | +0.34 (6) | +0.42 (8) | +0.57 (14) | +0.35 (12) |
| Ranch_Incised | -0.03 (22) | -0.03 (12) | -0.04 (15) | -0.03 (21) | -0.04 (20) |
| Walls_Engraved | -0.30 (14) | -0.27 (9) | -0.35 (10) | -0.30 (14) | -0.34 (12) |
| Wallace_Incised | +0.56 (7) | undefined (class at fewer than 2 rows) (0) | undefined (class at fewer than 2 rows) (0) | +0.56 (6) | undefined (class at fewer than 2 rows) (0) |
| Rhodes_Incised | -0.15 (15) | -0.27 (7) | -0.32 (8) | -0.15 (14) | -0.17 (11) |
| Vernon_Paul_Applique | +0.43 (10) | +0.24 (4) | +0.19 (4) | +0.45 (9) | +0.29 (6) |
| Hull_Engraved | -0.09 (8) | -0.23 (5) | -0.30 (5) | -0.08 (8) | -0.25 (5) |

Rows present in parentheses. Totals of the 18 rows with no addition: 79 to 1303 sherds; of the six large-addition rows: 424 to 1709; PFG parts of the matched sites: 79 to 1303.

## 2. The nine mixed sites: PFG part against Lipo addition

Log odds ratio is the Lipo addition minus the PFG part, so a positive value means the class is a larger share of Lipo's collection than of the 1940 collection at the same site. The per-site rows are the paired reading. The pooled rows contrast two collection-size-weighted mixtures (the pooled PFG parts against the pooled additions), which can differ even when no site differs within itself; they are not a site-adjusted effect.

### Pooled over the six large-addition sites (size-weighted)

| class | PFG parts: count, share | Lipo additions: count, share | log odds ratio, Lipo minus PFG, median [5%, 95%] |
|---|---|---|---|
| Parkin_Punctated | 787, 0.610 | 2633, 0.529 | -0.32 [-0.43, -0.22] |
| Barton/Kent/MPI | 339, 0.263 | 1713, 0.344 | +0.39 [+0.28, +0.51] |
| Painted | 134, 0.104 | 359, 0.072 | -0.40 [-0.57, -0.22] |
| Fortune_Noded | 8, 0.006 | 53, 0.011 | +0.53 [-0.04, +1.20] |
| Ranch_Incised | 9, 0.007 | 66, 0.013 | +0.64 [+0.10, +1.27] |
| Walls_Engraved | 3, 0.002 | 55, 0.011 | +1.51 [+0.69, +2.60] |
| Wallace_Incised | 0, 0.000 | 47, 0.009 | +3.99 [+1.84, +8.72] |
| Rhodes_Incised | 8, 0.006 | 22, 0.004 | -0.35 [-1.00, +0.36] |
| Vernon_Paul_Applique | 2, 0.002 | 20, 0.004 | +0.88 [-0.14, +2.23] |
| Hull_Engraved | 0, 0.000 | 5, 0.001 | +1.78 [-0.52, +6.54] |

### Pooled over all nine mixed sites (size-weighted)

| class | PFG parts: count, share | Lipo additions: count, share | log odds ratio, Lipo minus PFG, median [5%, 95%] |
|---|---|---|---|
| Parkin_Punctated | 2003, 0.592 | 2646, 0.528 | -0.26 [-0.33, -0.19] |
| Barton/Kent/MPI | 1090, 0.322 | 1714, 0.342 | +0.09 [+0.01, +0.17] |
| Painted | 213, 0.063 | 383, 0.076 | +0.21 [+0.07, +0.35] |
| Fortune_Noded | 29, 0.009 | 53, 0.011 | +0.21 [-0.17, +0.60] |
| Ranch_Incised | 24, 0.007 | 66, 0.013 | +0.62 [+0.24, +1.02] |
| Walls_Engraved | 11, 0.003 | 55, 0.011 | +1.22 [+0.70, +1.79] |
| Wallace_Incised | 0, 0.000 | 47, 0.009 | +4.94 [+2.80, +9.76] |
| Rhodes_Incised | 10, 0.003 | 22, 0.004 | +0.39 [-0.21, +1.05] |
| Vernon_Paul_Applique | 2, 0.001 | 20, 0.004 | +1.84 [+0.83, +3.19] |
| Hull_Engraved | 0, 0.000 | 5, 0.001 | +2.72 [+0.42, +7.42] |

### Per site, the three common classes, with the prior sensitivity

| site | addition | PFG part total | Lipo addition total | Parkin_Punctated: PFG share, Lipo share; log odds ratio at alpha 0.5 [5%, 95%]; medians at alpha 0.1, 1.0 | Barton/Kent/MPI: PFG share, Lipo share; log odds ratio at alpha 0.5 [5%, 95%]; medians at alpha 0.1, 1.0 | Painted: PFG share, Lipo share; log odds ratio at alpha 0.5 [5%, 95%]; medians at alpha 0.1, 1.0 |
|---|---|---|---|---|---|---|
| Belle_Meade | large | 254 | 1455 | 0.504, 0.546; +0.20 [-0.03, +0.41]; +0.18, +0.22 | 0.429, 0.363; -0.25 [-0.48, -0.03]; -0.27, -0.23 | 0.039, 0.022; -0.59 [-1.18, +0.04]; -0.58, -0.61 |
| Beck | large | 79 | 1154 | 0.405, 0.484; +0.39 [+0.01, +0.78]; +0.34, +0.46 | 0.544, 0.394; -0.50 [-0.88, -0.13]; -0.59, -0.42 | 0.025, 0.056; +0.80 [-0.19, +2.15]; +0.96, +0.65 |
| Nickel | large | 180 | 1120 | 0.633, 0.548; -0.30 [-0.58, -0.04]; -0.34, -0.25 | 0.228, 0.288; +0.34 [+0.04, +0.66]; +0.33, +0.35 | 0.122, 0.123; +0.02 [-0.37, +0.44]; +0.02, +0.03 |
| Rose_Mound | large | 245 | 738 | 0.629, 0.535; -0.35 [-0.60, -0.11]; -0.38, -0.32 | 0.265, 0.356; +0.44 [+0.18, +0.71]; +0.43, +0.45 | 0.094, 0.073; -0.26 [-0.68, +0.18]; -0.27, -0.26 |
| Cramor_Place | large | 138 | 286 | 0.659, 0.395; -1.03 [-1.38, -0.67]; -1.08, -0.97 | 0.261, 0.420; +0.73 [+0.36, +1.11]; +0.72, +0.73 | 0.051, 0.122; +0.96 [+0.31, +1.71]; +0.99, +0.92 |
| Castile_Landing | large | 394 | 220 | 0.680, 0.718; +0.15 [-0.15, +0.45]; +0.18, +0.12 | 0.114, 0.109; -0.06 [-0.51, +0.38]; -0.06, -0.05 | 0.178, 0.159; -0.14 [-0.52, +0.23]; -0.14, -0.14 |
| Hollywood | small | 201 | 25 | 0.303, 0.000; -3.98 [-8.83, -1.80]; -9.76, -3.01 | 0.363, 0.040; -2.57 [-4.53, -1.32]; -2.86, -2.37 | 0.274, 0.960; +2.54 [+1.77, +3.45]; +3.71, +1.97 |
| Barton_Ranch | small | 629 | 9 | 0.628, 1.000; +0.28 [-0.65, +1.35]; +2.17, -0.38 | 0.331, 0.000; -3.33 [-8.21, -1.14]; -9.02, -2.51 | 0.010, 0.000; +0.57 [-4.25, +2.91]; -5.01, +1.32 |
| Neeleys_Ferry | small | 1262 | 4 | 0.602, 1.000; -0.41 [-1.55, +0.77]; +1.43, -1.02 | 0.372, 0.000; -3.02 [-7.75, -0.79]; -8.45, -2.37 | 0.014, 0.000; +0.67 [-4.12, +2.95]; -4.65, +1.31 |

### Per site, all ten classes: counts in the PFG part and in the Lipo addition

| site | Parkin_Punctated | Barton/Kent/MPI | Painted | Fortune_Noded | Ranch_Incised | Walls_Engraved | Wallace_Incised | Rhodes_Incised | Vernon_Paul_Applique | Hull_Engraved |
|---|---|---|---|---|---|---|---|---|---|---|
| Belle_Meade | 128 + 795 | 109 + 528 | 10 + 32 | 0 + 12 | 3 + 30 | 2 + 25 | 0 + 15 | 2 + 11 | 0 + 5 | 0 + 2 |
| Beck | 32 + 558 | 43 + 455 | 2 + 65 | 0 + 10 | 1 + 20 | 1 + 18 | 0 + 12 | 0 + 8 | 0 + 7 | 0 + 1 |
| Nickel | 114 + 614 | 41 + 323 | 22 + 138 | 1 + 8 | 0 + 5 | 0 + 8 | 0 + 14 | 2 + 1 | 0 + 7 | 0 + 2 |
| Rose_Mound | 154 + 395 | 65 + 263 | 23 + 54 | 1 + 18 | 1 + 3 | 0 + 0 | 0 + 3 | 0 + 1 | 1 + 1 | 0 + 0 |
| Cramor_Place | 91 + 113 | 36 + 120 | 7 + 35 | 3 + 4 | 1 + 7 | 0 + 4 | 0 + 2 | 0 + 1 | 0 + 0 | 0 + 0 |
| Castile_Landing | 268 + 158 | 45 + 24 | 70 + 35 | 3 + 1 | 3 + 1 | 0 + 0 | 0 + 1 | 4 + 0 | 1 + 0 | 0 + 0 |
| Hollywood | 61 + 0 | 73 + 1 | 55 + 24 | 0 + 0 | 2 + 0 | 8 + 0 | 0 + 0 | 2 + 0 | 0 + 0 | 0 + 0 |
| Barton_Ranch | 395 + 9 | 208 + 0 | 6 + 0 | 16 + 0 | 4 + 0 | 0 + 0 | 0 + 0 | 0 + 0 | 0 + 0 | 0 + 0 |
| Neeleys_Ferry | 760 + 4 | 470 + 0 | 18 + 0 | 5 + 0 | 9 + 0 | 0 + 0 | 0 + 0 | 0 + 0 | 0 + 0 | 0 + 0 |

Cells are PFG part + Lipo addition.

## 3. The per-class table on PFG material alone

Moment F (common classes) and sites present (rare classes) as in the per-type table. The matched columns hold the 27 sites with a PFG part of at least 75 sherds and change only the collection (Holden_Lake has no PFG part; every PFG part clears the cutoff).

| class | 28 rows, current | 18 rows, no addition | 27 matched sites, current | 27 matched sites, PFG parts |
|---|---|---|---|---|
| Parkin_Punctated | F 0.165 | F 0.189 | F 0.145 | F 0.151 |
| Barton/Kent/MPI | F 0.117 | F 0.099 | F 0.076 | F 0.087 |
| Painted | F 0.199 | F 0.209 | F 0.195 | F 0.195 |
| Fortune_Noded | 15 of 28 sites | 6 of 18 sites | 14 of 27 sites | 12 of 27 sites |
| Ranch_Incised | 22 of 28 sites | 12 of 18 sites | 21 of 27 sites | 20 of 27 sites |
| Walls_Engraved | 14 of 28 sites | 9 of 18 sites | 14 of 27 sites | 12 of 27 sites |
| Wallace_Incised | 7 of 28 sites | 0 of 18 sites | 6 of 27 sites | 0 of 27 sites |
| Rhodes_Incised | 15 of 28 sites | 7 of 18 sites | 14 of 27 sites | 11 of 27 sites |
| Vernon_Paul_Applique | 10 of 28 sites | 4 of 18 sites | 9 of 27 sites | 6 of 27 sites |
| Hull_Engraved | 8 of 28 sites | 5 of 18 sites | 8 of 27 sites | 5 of 27 sites |

PFG parts of the nine mixed sites: Belle_Meade 254, Beck 79, Nickel 180, Rose_Mound 245, Cramor_Place 138, Castile_Landing 394, Hollywood 201, Barton_Ranch 629, Neeleys_Ferry 1262.

## Provenance: source rows used

| assemblage | cpl-pfg id | PFG rows | PFG part | addition |
|---|---|---|---|---|
| Barton_Ranch | 11-O-10 | 11-O-10 | 629 | 9 |
| Beck | 13-O-7 | 13-O-7 | 79 | 1154 |
| Belle_Meade | 13-O-5 | 13-O-5 | 254 | 1455 |
| Big_Eddy |  | matrix row (PFGData.xlsx) | 229 | 0 |
| Castile_Landing | 13-N-21 | 13-N-21 | 394 | 220 |
| Clay_Hill |  | matrix row (PFGData.xlsx) | 181 | 0 |
| Commerce | 13-O-11 | 13-O-11 | 128 | 0 |
| Cramor_Place | 12-O-5 | 12-O-5 | 138 | 286 |
| Cummins |  | matrix row (PFGData.xlsx) | 261 | 0 |
| Davis | 13-N-5 | 13-N-5 | 79 | 0 |
| Fortune |  | matrix row (PFGData.xlsx) | 200 | 0 |
| Grant |  | matrix row (PFGData.xlsx) | 93 | 0 |
| Holden_Lake | Holden_Lake | none | 0 | 360 |
| Hollywood | 13-O-10 | 13-O-10 | 201 | 25 |
| Irby | 13-P-10 | 13-P-10 | 91 | 0 |
| Kent_Place | 13-N-4 | 13-N-4/B,C,D,E | 241 | 0 |
| Lake_Cormorant | 13-P-8 | 13-P-8 | 192 | 0 |
| Mound_Place |  | matrix row (PFGData.xlsx) | 168 | 0 |
| Neeleys_Ferry | 11-N-4 | 11-N-4 | 1262 | 4 |
| Nickel | 13-N-15 | 13-N-15/B; 13-N-15/C,D,E | 180 | 1120 |
| Parkin | 11-N-1 | 11-N-1 | 1303 | 0 |
| Rose_Mound | 12-N-3 | 12-N-3/A&B; 12-N-3/C | 245 | 738 |
| Starkley | 13-N-16 | 13-N-16 | 171 | 0 |
| Turnbow |  | matrix row (PFGData.xlsx) | 128 | 0 |
| Vernon_Paul | 11-N-9 | 11-N-9 | 758 | 0 |
| Walls | 13-P-1 | 13-P-1 | 360 | 0 |
| Williamson |  | matrix row (PFGData.xlsx) | 658 | 0 |
| Woodlyn |  | matrix row (PFGData.xlsx) | 107 | 0 |
