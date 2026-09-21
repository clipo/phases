# Where does the drift shortfall live?

Basin phase set, 43 assemblages, k = 5, 200 drift realizations, calibrated pooled model.
Observed F_ST 0.0280 against a drift median of 0.0037, an overall shortfall of 7.6x.

## Per-cluster contribution

| cluster | n | sherds | observed term | share of observed F_ST | drift median term | excess |
|---|---|---|---|---|---|---|
| 1 | 4 | 762 | 0.0121 | 43% | 0.0001 | **82.1x** |
| 2 | 3 | 433 | 0.0032 | 11% | 0.0002 | **17.9x** |
| 3 | 11 | 6,756 | 0.0067 | 24% | 0.0005 | **13.1x** |
| 0 | 10 | 15,580 | 0.0056 | 20% | 0.0015 | **3.9x** |
| 4 | 15 | 12,860 | 0.0005 | 2% | 0.0010 | **0.5x** |

Excess ranges from 0.5x to 82.1x, a ratio of 180.8.

**CONCENTRATED: one or two clusters carry the excess.**

## Leave one cluster out

| dropped | observed F_ST | drift median | shortfall |
|---|---|---|---|
| 0.0 | 0.0289 | 0.0014 | 20.2x |
| 1.0 | 0.0163 | 0.0037 | 4.5x |
| 2.0 | 0.0253 | 0.0036 | 7.0x |
| 3.0 | 0.0252 | 0.0035 | 7.1x |
| 4.0 | 0.0429 | 0.0028 | 15.3x |

If the excess were one cluster's closure, dropping that cluster would collapse the shortfall toward 1. Read the column for that.

## Cluster membership

- **cluster 1** (4): Dundee, Parchman, Salomon, West_Mounds
- **cluster 2** (3): Carson_Lake, Notgrass, Upper_Nodena
- **cluster 3** (11): Big_Eddy, Castile_Landing, Clay_Hill, Connor, Cramor_Place, Davis, Grant, Kent_Place, Nickel, Soudan, Starkley
- **cluster 0** (10): Barton_Ranch, Cummins, Fortune, Holden_Lake, Neeleys_Ferry, Parkin, Rose_Mound, Turnbow, Vernon_Paul, Williamson
- **cluster 4** (15): Beck, Belle_Meade, Cheatham, Chuccalissa, Commerce, Hollywood, Irby, Lake_Cormorant, Mound_Place, Norfolk, Pouncey, Wall, Walls, Woodlyn, Young