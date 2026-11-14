# Is the finer-scale excess a sampling effect?

Produced by `analyses/85_excess_robustness.py`. Basin, 28 assemblages, k = 5 spatial clusters (the division of analysis 72), calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), 300 drift runs, 2000 draws for (a) and (b).

| cluster | assemblages | sherds | observed term | (a) sampling alone, median [95%] | drift median | (b) excess over drift median, posterior median [95%] | share of drift runs reaching the observed term |
|---|---|---|---|---|---|---|---|
| 1 (Clay_Hill, Davis, Grant, Kent_Place, Starkley) | 5 | 765 | 0.0150 | 0.00004 [0.00001, 0.00028] | 0.0004 | 32.6 [26.5, 39.1] | 0.0% |
| 2 (Irby, Lake_Cormorant, Mound_Place, Walls, Woodlyn) | 5 | 918 | 0.0044 | 0.00004 [0.00001, 0.00024] | 0.0005 | 8.0 [6.3, 9.8] | 0.0% |
| 0 (Barton_Ranch, Cummins, Fortune, Holden_Lake, Neeleys_Ferry, Parkin, Rose_Mound, Turnbow, Vernon_Paul, Williamson) | 10 | 6,555 | 0.0052 | 0.00002 [0.00000, 0.00014] | 0.0012 | 4.3 [3.6, 5.0] | 4.7% |
| 4 (Beck, Belle_Meade, Commerce, Hollywood) | 4 | 3,296 | 0.0034 | 0.00004 [0.00001, 0.00022] | 0.0008 | 4.0 [2.6, 5.8] | 7.3% |
| 3 (Big_Eddy, Castile_Landing, Cramor_Place, Nickel) | 4 | 2,567 | 0.0028 | 0.00004 [0.00001, 0.00022] | 0.0010 | 2.7 [1.9, 3.6] | 11.0% |

## (c) Dropping one assemblage at a time

| cluster | assemblage dropped | observed term | drift median | excess |
|---|---|---|---|---|
| 1 | Clay_Hill | 0.0140 | 0.0004 | 33.7 |
| 1 | Davis | 0.0142 | 0.0004 | 34.9 |
| 1 | Grant | 0.0116 | 0.0004 | 29.1 |
| 1 | Kent_Place | 0.0109 | 0.0003 | 31.7 |
| 1 | Starkley | 0.0117 | 0.0004 | 33.1 |
| 2 | Irby | 0.0042 | 0.0005 | 8.0 |
| 2 | Lake_Cormorant | 0.0038 | 0.0005 | 7.5 |
| 2 | Mound_Place | 0.0045 | 0.0005 | 8.7 |
| 2 | Walls | 0.0028 | 0.0004 | 7.8 |
| 2 | Woodlyn | 0.0038 | 0.0005 | 7.8 |
| 0 | Barton_Ranch | 0.0051 | 0.0013 | 4.0 |
| 0 | Cummins | 0.0050 | 0.0012 | 4.2 |
| 0 | Fortune | 0.0051 | 0.0012 | 4.3 |
| 0 | Holden_Lake | 0.0066 | 0.0013 | 5.1 |
| 0 | Neeleys_Ferry | 0.0051 | 0.0012 | 4.5 |
| 0 | Parkin | 0.0055 | 0.0012 | 4.7 |
| 0 | Rose_Mound | 0.0063 | 0.0013 | 4.8 |
| 0 | Turnbow | 0.0056 | 0.0012 | 4.7 |
| 0 | Vernon_Paul | 0.0049 | 0.0012 | 4.1 |
| 0 | Williamson | 0.0050 | 0.0012 | 4.1 |
| 4 | Beck | 0.0022 | 0.0007 | 3.0 |
| 4 | Belle_Meade | 0.0046 | 0.0009 | 4.9 |
| 4 | Commerce | 0.0029 | 0.0008 | 3.6 |
| 4 | Hollywood | 0.0034 | 0.0008 | 4.0 |
| 3 | Big_Eddy | 0.0023 | 0.0011 | 2.1 |
| 3 | Castile_Landing | 0.0018 | 0.0010 | 1.8 |
| 3 | Cramor_Place | 0.0037 | 0.0009 | 4.0 |
| 3 | Nickel | 0.0034 | 0.0006 | 5.6 |

## Reading

(a) says how large each cluster's term would be from finite samples alone. (b) carries the uncertainty in the observed counts into the excess; a posterior that stays above 1 says the excess is not an accident of which sherds were collected. (c) says whether the excess belongs to the cluster or to one collection in it.
