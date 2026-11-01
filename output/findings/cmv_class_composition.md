# What the southeast-Missouri comparison classes actually are

Produced by `analyses/66_cmv_class_composition.py`. 39 assemblages, 42 classes, 8958 sherds: the set analysis 47 calls `cmv`.

## Composition

| class | share of all sherds |
|---|---:|
| Mississippi Plain | 81.51% |
| Bell Plain | 7.51% |
| Varney Red | 3.82% |
| Old Town Red Filmed | 1.72% |
| Wickliffe Plain | 0.88% |
| Parkin Punctate | 0.87% |
| Thin ware | 0.55% |
| Wickliffe Incised | 0.44% |
| Matthews Incised | 0.38% |
| O'Byam Engraved | 0.37% |
| Fine ware | 0.33% |
| Kimmswick Plain | 0.23% |

Plain wares (Mississippi Plain, Bell Plain, Wickliffe Plain, Kimmswick Plain) hold **8075 of 8958 sherds, 90.1 percent**. The non-plain remainder is 883 sherds (9.9 percent), a median of 16 per assemblage over 38 classes; **12 of 39 assemblages** reach the 30-sherd minimum the study imposes when it is applied to that material.

### What the non-plain sherds are

| class | share of non-plain |
|---|---:|
| Varney Red | 38.7% |
| Old Town Red Filmed | 17.4% |
| Parkin Punctate | 8.8% |
| Thin ware | 5.5% |
| Wickliffe Incised | 4.4% |
| Matthews Incised | 3.9% |
| O'Byam Engraved | 3.7% |
| Fine ware | 3.4% |

## Where the between-cluster differentiation lives

Cultural $F_{ST}$ across the 4 spatial clusters, by subset. The null resamples every assemblage from the regional pooled profile at its own observed count, so it carries no spatial structure at all (400 draws, seed 3); it separates sampling from everything else, and is NOT the paper's drift null.

| subset | classes | sherds | $F_{ST}$ | pooled-profile null (median, 95%) |
|---|---:|---:|---:|---|
| all classes, as analyzed | 42 | 8958 | 0.0337 | 0.0003 [0.0001, 0.0007] |
| plain wares only | 4 | 8075 | 0.0652 | 0.0003 [0.0001, 0.0010] |
| Mississippi Plain vs Bell Plain only | 2 | 7975 | 0.0755 | 0.0003 [0.0000, 0.0011] |
| non-plain only | 38 | 883 | 0.3245 | 0.0032 [0.0018, 0.0060] |
| non-plain, named types only | 24 | 734 | 0.3506 | 0.0038 [0.0017, 0.0077] |
| non-plain, named, minus Varney Red | 23 | 392 | 0.1796 | 0.0072 [0.0036, 0.0139] |

## Reading

The reported comparison value of 0.0337 is a blend in which 90 percent of the weight is undecorated utility ware. The plain component on its own differentiates more strongly (0.0652), and the coarse-versus-fine paste ratio alone more strongly still (0.0755). That is a manufacturing variable, so a between-group difference in it is the signature Table 1 files under environmental patchiness rather than under bounded interaction.

The non-plain subset differentiates far more (0.3245), and that is not pure sampling: the pooled-profile null for the same subset sits two orders of magnitude below it. But it is not a clean stylistic measurement either. Varney Red, a chronologically diagnostic red-filmed ware that the Woodland exclusion does not catch, is 38.7 percent of the non-plain sherds, and removing it moves the value to 0.1796. 29 of the 38 non-plain classes occur in only one of the four spatial clusters while holding 302 sherds between them, so rarity and spatial restriction are confounded at this sample size and cannot be separated.

The conclusion is not that the comparison set is differentiated or undifferentiated. It is that the set cannot resolve the question the basin set is being asked, which is the paper's own empirical-sufficiency criterion applied to itself.

## A defect in the source table

These class names collide under case normalization, so one class is carried as two. No result here turns on them; they should be merged in the source.

| names | sherds |
|---|---|
| 'Brown Slip', 'Brown slip' | 2, 2 |
| 'Fine Ware', 'Fine ware' | 4, 30 |
| 'Thin Ware', 'Thin ware' | 5, 49 |
