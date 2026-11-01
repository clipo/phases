# The attribute scale: Mainfort's rim attributes (exploratory)

Produced by `analyses/100_rim_attributes.py`. Counts are Table 1 of Mainfort (2003b), *Southeastern Archaeology* 22(1):33-46, as printed: 55 assemblages, 7,195 rims, six attributes, each a separate yes-or-no observation on a rim. Paste, lip form and bevel degree come from the file Mainfort sent the author and are not in the paper; they are in the last part. Every statistic is a posterior median with a 95 percent credible interval over 2,000 draws of each assemblage's frequencies (Beta or Dirichlet, counts + 1/2).

## A. How each attribute varies across the region (48 assemblages with at least 45 rims)

Coordinates: 21 from the analysis matrix, 14 from the rim file (rough), 13 from the Mainfort 2003a site file. Distances are straight-line; the pairs span 0 to 220 km.

| attribute | percent of rims: median assemblage (lowest to highest) | rank correlation of difference with distance | with latitude (north positive) | with longitude (east positive) |
|---|---|---|---|---|
| interior bevel | 17.4 (0.0 to 80.3) | +0.34 (+0.30 to +0.39) | +0.54 (+0.49 to +0.58) | +0.73 (+0.69 to +0.77) |
| beaded applique strip | 4.3 (0.0 to 39.3) | +0.40 (+0.35 to +0.46) | +0.70 (+0.64 to +0.76) | +0.78 (+0.73 to +0.83) |
| exterior notching | 12.1 (0.0 to 35.0) | +0.11 (+0.05 to +0.17) | +0.35 (+0.27 to +0.42) | +0.53 (+0.46 to +0.59) |
| notched top | 1.3 (0.0 to 5.8) | +0.02 (-0.04 to +0.09) | -0.01 (-0.18 to +0.17) | +0.04 (-0.14 to +0.22) |
| band of nodes | 0.9 (0.0 to 9.9) | +0.03 (-0.04 to +0.10) | +0.26 (+0.11 to +0.41) | +0.38 (+0.22 to +0.52) |
| widely spaced nodes | 0.5 (0.0 to 3.6) | -0.02 (-0.08 to +0.06) | -0.04 (-0.23 to +0.16) | +0.04 (-0.16 to +0.24) |

A positive correlation with distance means assemblages farther apart differ more in that attribute. A positive correlation with latitude means the attribute is commoner to the north.

## B. Do the attributes vary together?

Rank correlation between two attributes' frequencies across the same assemblages.

| attribute pair | rank correlation (95% CrI) |
|---|---|
| interior bevel and beaded applique strip | +0.69 (+0.61 to +0.75) |
| interior bevel and exterior notching | +0.54 (+0.46 to +0.62) |
| interior bevel and notched top | -0.15 (-0.35 to +0.02) |
| interior bevel and band of nodes | +0.29 (+0.11 to +0.45) |
| interior bevel and widely spaced nodes | +0.08 (-0.14 to +0.28) |
| beaded applique strip and exterior notching | +0.42 (+0.31 to +0.52) |
| beaded applique strip and notched top | -0.04 (-0.25 to +0.15) |
| beaded applique strip and band of nodes | +0.34 (+0.17 to +0.49) |
| beaded applique strip and widely spaced nodes | -0.02 (-0.23 to +0.19) |
| exterior notching and notched top | +0.06 (-0.15 to +0.27) |
| exterior notching and band of nodes | +0.23 (+0.04 to +0.41) |
| exterior notching and widely spaced nodes | +0.07 (-0.16 to +0.28) |
| notched top and band of nodes | +0.05 (-0.17 to +0.28) |
| notched top and widely spaced nodes | -0.07 (-0.31 to +0.19) |
| band of nodes and widely spaced nodes | +0.11 (-0.12 to +0.34) |

## C. Do the phase lines divide the attributes? (19 assemblages of the primary analysis)

19 of the 28 analyzed assemblages have at least 25 rims: Kent 6, Parkin 9, Walls 4 (2,065 rims; 27 to 277 per assemblage, 1 below Mainfort's minimum of 45). Without rim data: Commerce, Cramor_Place, Cummins, Grant, Holden_Lake, Hollywood, Lake_Cormorant, Nickel, Turnbow.

Percent of rims, pooled by phase:

| attribute | class | Kent | Parkin | Walls |
|---|---|---:|---:|---:|
| interior bevel | beveled | 16.9 | 3.9 | 16.6 |
| beaded applique strip | present | 1.3 | 1.3 | 1.5 |
| exterior notching | present | 8.8 | 2.2 | 15.0 |
| notched top | present | 1.4 | 0.7 | 1.3 |
| band of nodes | present | 1.0 | 0.4 | 0.4 |
| widely spaced nodes | present | 0.4 | 0.5 | 0.7 |

The paper's phase-line comparison (analysis 86), attribute by attribute. A probability near one half means the phase lines divide the attribute no better than arbitrary lines.

| attribute | F_ST, posterior median (95% CrI) | median F_ST, random-center / compact | P(lines separate better than random-center / compact) | boundary excess | P(larger than random-center / compact) |
|---|---|---|---|---:|---|
| decorated types, same assemblages | 0.0148 (0.0117 to 0.0183) | 0.0150 / 0.0150 | 0.47 / 0.44 | +23.5 | 0.88 / 0.87 |
| interior bevel | 0.0345 (0.0235 to 0.0466) | 0.0304 / 0.0359 | 0.62 / 0.39 | +6.7 | 0.79 / 0.74 |
| beaded applique strip | 0.0007 (0.0000 to 0.0038) | 0.0024 / 0.0023 | 0.17 / 0.15 | +0.3 | 0.68 / 0.63 |
| exterior notching | 0.0306 (0.0181 to 0.0455) | 0.0262 / 0.0299 | 0.62 / 0.48 | +5.5 | 0.77 / 0.68 |
| notched top | 0.0010 (0.0000 to 0.0047) | 0.0010 / 0.0009 | 0.46 / 0.47 | +0.5 | 0.68 / 0.69 |
| band of nodes | 0.0008 (0.0000 to 0.0042) | 0.0011 / 0.0011 | 0.39 / 0.37 | +0.1 | 0.56 / 0.53 |
| widely spaced nodes | 0.0009 (0.0000 to 0.0046) | 0.0007 / 0.0007 | 0.57 / 0.57 | +0.1 | 0.53 / 0.49 |
| the six attributes together (mean F_ST) | 0.0117 (0.0089 to 0.0149) | 0.0103 / 0.0117 | 0.68 / 0.44 | | |

## D. Does any attribute go with the difference in decoration that distance leaves?

For the 19 assemblages of part C (171 pairs). The difference in decoration between two assemblages is half the summed absolute difference of their decorated-class proportions. It is correlated with river distance at +0.18 (+0.14 to +0.22) and with straight-line distance at +0.41 (+0.37 to +0.45). Each row asks whether pairs that differ more in an attribute also differ more in decoration, first as it stands and then with distance removed.

| attribute | rank correlation of decoration difference with attribute difference | with river distance removed | with straight-line distance removed |
|---|---|---|---|
| interior bevel (published) | +0.14 (+0.02 to +0.24) | +0.09 (-0.03 to +0.20) | +0.06 (-0.06 to +0.18) |
| exterior notching (published) | +0.07 (-0.02 to +0.16) | +0.04 (-0.05 to +0.13) | +0.02 (-0.07 to +0.11) |
| paste: Bell share of rims (file only) | +0.12 (+0.04 to +0.21) | +0.05 (-0.03 to +0.15) | +0.04 (-0.05 to +0.13) |

Interior bevel is read as increasing through time. A clear positive value in its row would say that assemblages of different age differ more in decoration than distance alone predicts, which is the reading that the unexplained difference is chronological. A value near zero says bevel, and so age as bevel measures it, does not go with that difference.


## The attributes in Mainfort's file that the paper does not print

Paste is the paste of the rim sherds (author); lip form and bevel degree are as recorded. Counts are the file's percentages times N. Use with that in mind.

Percent of rims, pooled by phase:

| attribute | class | Kent | Parkin | Walls |
|---|---|---:|---:|---:|
| paste: Bell share of rims in the two pastes | Bell | 33.4 | 21.2 | 79.1 |
| paste, three classes | Neeley's Ferry | 46.0 | 52.0 | 19.6 |
| paste, three classes | Bell | 23.0 | 14.0 | 74.1 |
| paste, three classes | other | 31.0 | 34.0 | 6.3 |
| lip form | flat | 38.4 | 28.3 | 36.8 |
| interior bevel, by degree | sharp | 5.8 | 1.7 | 6.1 |
| interior bevel, by degree | moderate | 7.7 | 0.9 | 6.8 |
| interior bevel, by degree | slight | 3.3 | 1.3 | 3.7 |
| interior bevel, by degree | none | 83.1 | 96.1 | 83.4 |

| attribute | F_ST, posterior median (95% CrI) | median F_ST, random-center / compact | P(lines separate better than random-center / compact) | boundary excess | P(larger than random-center / compact) |
|---|---|---|---|---:|---|
| paste: Bell share of rims in the two pastes | 0.2250 (0.1852 to 0.2643) | 0.1316 / 0.1629 | 0.95 / 0.91 | +16.7 | 0.85 / 0.78 |
| paste, three classes | 0.1232 (0.1039 to 0.1422) | 0.0632 / 0.0752 | 0.92 / 0.87 | +11.2 | 0.77 / 0.76 |
| lip form | 0.0098 (0.0034 to 0.0187) | 0.0069 / 0.0090 | 0.69 / 0.53 | +1.2 | 0.67 / 0.53 |
| interior bevel, by degree | 0.0207 (0.0136 to 0.0285) | 0.0195 / 0.0217 | 0.55 / 0.36 | +6.9 | 0.76 / 0.70 |

## Reading

**A.** interior bevel: difference grows with distance; beaded applique strip: difference grows with distance; exterior notching: difference grows with distance.

**B.** 6 of the 15 attribute pairs have a rank correlation whose interval excludes zero: interior bevel and beaded applique strip (+0.69); interior bevel and exterior notching (+0.54); interior bevel and band of nodes (+0.29); beaded applique strip and exterior notching (+0.42); beaded applique strip and band of nodes (+0.34); exterior notching and band of nodes (+0.23).

**C.** Against random-center divisions, the phase lines separate the six published attributes with probabilities from 0.17 to 0.62, and the six together with 0.68; the decorated types on the same assemblages give 0.47. Of the unpublished attributes, paste: Bell share of rims in the two pastes 0.95; paste, three classes 0.92; lip form 0.69; interior bevel, by degree 0.55.

A restriction on learning at the phase lines predicts high probabilities in every attribute and attributes that vary together. Variation that arises locally predicts probabilities scattered around one half and attributes that vary independently. Cautions: part C rests on 19 assemblages; rim form and paste may answer to vessel function, cost and time as well as to who learned from whom; and the table gives each attribute's total, not the attributes of each rim, so it cannot say which attributes occur on the same vessels.
