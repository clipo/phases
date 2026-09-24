# Where does the drift shortfall live?

Basin phase set, 28 assemblages, k = 5, 300 drift realizations, calibrated pooled model.
Observed F_ST 0.0206 against a drift median of 0.0041, an overall shortfall of 5.1x.

## Per-cluster contribution

| cluster | n | sherds | observed term | share of observed F_ST | drift median term | excess |
|---|---|---|---|---|---|---|
| 2 | 1 | 261 | 0.0026 | 13% | 0.0002 | **14.4x** |
| 1 | 5 | 918 | 0.0044 | 21% | 0.0005 | **8.5x** |
| 4 | 8 | 3,880 | 0.0055 | 27% | 0.0009 | **6.5x** |
| 0 | 9 | 6,294 | 0.0045 | 22% | 0.0010 | **4.3x** |
| 3 | 5 | 2,748 | 0.0035 | 17% | 0.0009 | **4.1x** |

Excess ranges from 4.1x to 14.4x, a ratio of 3.5.

**CONCENTRATED: one or two clusters carry the excess.**

## Leave one cluster out

| dropped | observed F_ST | drift median | shortfall |
|---|---|---|---|
| 0.0 | 0.0209 | 0.0033 | 6.3x |
| 1.0 | 0.0172 | 0.0039 | 4.5x |
| 2.0 | 0.0182 | 0.0041 | 4.4x |
| 3.0 | 0.0205 | 0.0037 | 5.5x |
| 4.0 | 0.0187 | 0.0041 | 4.5x |

If the excess were one cluster's closure, dropping that cluster would collapse the shortfall toward 1. Read the column for that.

## Cluster membership

- **cluster 2** (1): Cummins
- **cluster 1** (5): Irby, Lake_Cormorant, Mound_Place, Walls, Woodlyn
- **cluster 4** (8): Beck, Belle_Meade, Commerce, Davis, Grant, Hollywood, Kent_Place, Starkley
- **cluster 0** (9): Barton_Ranch, Fortune, Holden_Lake, Neeleys_Ferry, Parkin, Rose_Mound, Turnbow, Vernon_Paul, Williamson
- **cluster 3** (5): Big_Eddy, Castile_Landing, Clay_Hill, Cramor_Place, Nickel