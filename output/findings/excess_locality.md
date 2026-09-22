# Where does the drift shortfall live?

Basin phase set, 28 assemblages, k = 5, 300 drift realizations, calibrated pooled model.
Observed F_ST 0.0307 against a drift median of 0.0067, an overall shortfall of 4.6x.

## Per-cluster contribution

| cluster | n | sherds | observed term | share of observed F_ST | drift median term | excess |
|---|---|---|---|---|---|---|
| 1 | 5 | 765 | 0.0150 | 49% | 0.0005 | **30.1x** |
| 2 | 5 | 918 | 0.0044 | 14% | 0.0006 | **7.0x** |
| 0 | 10 | 6,555 | 0.0052 | 17% | 0.0013 | **4.0x** |
| 4 | 4 | 3,296 | 0.0034 | 11% | 0.0014 | **2.3x** |
| 3 | 4 | 2,567 | 0.0028 | 9% | 0.0013 | **2.2x** |

Excess ranges from 2.2x to 30.1x, a ratio of 13.8.

**CONCENTRATED: one or two clusters carry the excess.**

## Leave one cluster out

| dropped | observed F_ST | drift median | shortfall |
|---|---|---|---|
| 0.0 | 0.0363 | 0.0063 | 5.7x |
| 1.0 | 0.0161 | 0.0062 | 2.6x |
| 2.0 | 0.0281 | 0.0064 | 4.4x |
| 3.0 | 0.0338 | 0.0057 | 5.9x |
| 4.0 | 0.0349 | 0.0051 | 6.8x |

If the excess were one cluster's closure, dropping that cluster would collapse the shortfall toward 1. Read the column for that.

## Cluster membership

- **cluster 1** (5): Clay_Hill, Davis, Grant, Kent_Place, Starkley
- **cluster 2** (5): Irby, Lake_Cormorant, Mound_Place, Walls, Woodlyn
- **cluster 0** (10): Barton_Ranch, Cummins, Fortune, Holden_Lake, Neeleys_Ferry, Parkin, Rose_Mound, Turnbow, Vernon_Paul, Williamson
- **cluster 4** (4): Beck, Belle_Meade, Commerce, Hollywood
- **cluster 3** (4): Big_Eddy, Castile_Landing, Cramor_Place, Nickel