# Is the finer-scale excess a sampling effect?

Produced by `analyses/85_excess_robustness.py`. Basin, 28 assemblages, k = 5 spatial clusters (the division of analysis 72), calibrated pooled-profile cell (2000 learners, innovation 0.002, mixing 0.005), 300 drift runs, 2000 draws for (a) and (b).

| cluster | assemblages | sherds | observed term | (a) sampling alone, median [95%] | drift median | (b) excess over drift, posterior median [95%] | P(excess > 1) |
|---|---|---|---|---|---|---|---|
| 1 (Clay_Hill, Davis, Grant, Kent_Place, Starkley) | 5 | 765 | 0.0150 | 0.00004 [0.00001, 0.00028] | 0.0005 | 29.0 [23.6, 34.9] | 1.00 |
| 2 (Irby, Lake_Cormorant, Mound_Place, Walls, Woodlyn) | 5 | 918 | 0.0044 | 0.00004 [0.00001, 0.00024] | 0.0006 | 6.9 [5.4, 8.4] | 1.00 |
| 0 (Barton_Ranch, Cummins, Fortune, Holden_Lake, Neeleys_Ferry, Parkin, Rose_Mound, Turnbow, Vernon_Paul, Williamson) | 10 | 6,555 | 0.0052 | 0.00002 [0.00000, 0.00014] | 0.0013 | 3.9 [3.4, 4.6] | 1.00 |
| 4 (Beck, Belle_Meade, Commerce, Hollywood) | 4 | 3,296 | 0.0034 | 0.00004 [0.00001, 0.00022] | 0.0014 | 2.3 [1.5, 3.3] | 1.00 |
| 3 (Big_Eddy, Castile_Landing, Cramor_Place, Nickel) | 4 | 2,567 | 0.0028 | 0.00004 [0.00001, 0.00022] | 0.0013 | 2.2 [1.5, 2.9] | 1.00 |

## (c) Dropping one assemblage at a time

| cluster | assemblage dropped | observed term | drift median | excess |
|---|---|---|---|---|
| 1 | Clay_Hill | 0.0140 | 0.0005 | 29.6 |
| 1 | Davis | 0.0142 | 0.0005 | 28.0 |
| 1 | Grant | 0.0116 | 0.0005 | 22.7 |
| 1 | Kent_Place | 0.0109 | 0.0004 | 29.0 |
| 1 | Starkley | 0.0117 | 0.0004 | 29.2 |
| 2 | Irby | 0.0042 | 0.0006 | 6.8 |
| 2 | Lake_Cormorant | 0.0038 | 0.0006 | 6.7 |
| 2 | Mound_Place | 0.0045 | 0.0006 | 7.6 |
| 2 | Walls | 0.0028 | 0.0004 | 6.5 |
| 2 | Woodlyn | 0.0038 | 0.0006 | 6.5 |
| 0 | Barton_Ranch | 0.0051 | 0.0014 | 3.6 |
| 0 | Cummins | 0.0050 | 0.0014 | 3.6 |
| 0 | Fortune | 0.0051 | 0.0014 | 3.8 |
| 0 | Holden_Lake | 0.0066 | 0.0015 | 4.4 |
| 0 | Neeleys_Ferry | 0.0051 | 0.0013 | 3.9 |
| 0 | Parkin | 0.0055 | 0.0015 | 3.7 |
| 0 | Rose_Mound | 0.0063 | 0.0014 | 4.5 |
| 0 | Turnbow | 0.0056 | 0.0013 | 4.2 |
| 0 | Vernon_Paul | 0.0049 | 0.0013 | 3.7 |
| 0 | Williamson | 0.0050 | 0.0014 | 3.5 |
| 4 | Beck | 0.0022 | 0.0018 | 1.3 |
| 4 | Belle_Meade | 0.0046 | 0.0017 | 2.8 |
| 4 | Commerce | 0.0029 | 0.0014 | 2.1 |
| 4 | Hollywood | 0.0034 | 0.0015 | 2.3 |
| 3 | Big_Eddy | 0.0023 | 0.0014 | 1.7 |
| 3 | Castile_Landing | 0.0018 | 0.0014 | 1.3 |
| 3 | Cramor_Place | 0.0037 | 0.0014 | 2.6 |
| 3 | Nickel | 0.0034 | 0.0009 | 3.7 |

## Reading

(a) says how large each cluster's term would be from finite samples alone. (b) carries the uncertainty in the observed counts into the excess; a posterior that stays above 1 says the excess is not an accident of which sherds were collected. (c) says whether the excess belongs to the cluster or to one collection in it.
