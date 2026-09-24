# Is the finer-scale excess a sampling effect?

Produced by `analyses/85_excess_robustness.py`. Basin, 28 assemblages, k = 5 spatial clusters (the division of analysis 72), calibrated pooled-profile cell (2000 learners, innovation 0.001, mixing 0.02), 300 drift runs, 2000 draws for (a) and (b).

| cluster | assemblages | sherds | observed term | (a) sampling alone, median [95%] | drift median | (b) excess over drift median, posterior median [95%] | share of drift runs reaching the observed term |
|---|---|---|---|---|---|---|---|
| 2 (Cummins) | 1 | 261 | 0.0026 | 0.00005 [0.00001, 0.00026] | 0.0002 | 13.6 [8.2, 19.9] | 0.0% |
| 1 (Irby, Lake_Cormorant, Mound_Place, Walls, Woodlyn) | 5 | 918 | 0.0044 | 0.00004 [0.00001, 0.00024] | 0.0005 | 8.4 [6.6, 10.3] | 1.0% |
| 4 (Beck, Belle_Meade, Commerce, Davis, Grant, Hollywood, Kent_Place, Starkley) | 8 | 3,880 | 0.0055 | 0.00003 [0.00000, 0.00019] | 0.0009 | 6.2 [4.6, 8.3] | 0.3% |
| 0 (Barton_Ranch, Fortune, Holden_Lake, Neeleys_Ferry, Parkin, Rose_Mound, Turnbow, Vernon_Paul, Williamson) | 9 | 6,294 | 0.0045 | 0.00003 [0.00000, 0.00014] | 0.0010 | 4.3 [3.7, 5.0] | 6.3% |
| 3 (Big_Eddy, Castile_Landing, Clay_Hill, Cramor_Place, Nickel) | 5 | 2,748 | 0.0035 | 0.00004 [0.00000, 0.00022] | 0.0009 | 4.1 [3.1, 5.1] | 5.0% |

## (c) Dropping one assemblage at a time

| cluster | assemblage dropped | observed term | drift median | excess |
|---|---|---|---|---|
| 1 | Irby | 0.0042 | 0.0005 | 8.6 |
| 1 | Lake_Cormorant | 0.0038 | 0.0005 | 8.1 |
| 1 | Mound_Place | 0.0045 | 0.0005 | 9.3 |
| 1 | Walls | 0.0028 | 0.0004 | 7.1 |
| 1 | Woodlyn | 0.0038 | 0.0005 | 7.9 |
| 4 | Beck | 0.0060 | 0.0008 | 7.4 |
| 4 | Belle_Meade | 0.0111 | 0.0010 | 11.6 |
| 4 | Commerce | 0.0049 | 0.0009 | 5.6 |
| 4 | Davis | 0.0055 | 0.0009 | 6.4 |
| 4 | Grant | 0.0049 | 0.0008 | 5.8 |
| 4 | Hollywood | 0.0047 | 0.0009 | 5.4 |
| 4 | Kent_Place | 0.0048 | 0.0009 | 5.4 |
| 4 | Starkley | 0.0048 | 0.0008 | 5.6 |
| 0 | Barton_Ranch | 0.0044 | 0.0011 | 4.2 |
| 0 | Fortune | 0.0045 | 0.0010 | 4.3 |
| 0 | Holden_Lake | 0.0056 | 0.0012 | 4.7 |
| 0 | Neeleys_Ferry | 0.0043 | 0.0011 | 3.8 |
| 0 | Parkin | 0.0050 | 0.0012 | 4.3 |
| 0 | Rose_Mound | 0.0055 | 0.0011 | 4.8 |
| 0 | Turnbow | 0.0049 | 0.0011 | 4.5 |
| 0 | Vernon_Paul | 0.0044 | 0.0010 | 4.2 |
| 0 | Williamson | 0.0043 | 0.0011 | 3.8 |
| 3 | Big_Eddy | 0.0027 | 0.0009 | 3.2 |
| 3 | Castile_Landing | 0.0028 | 0.0008 | 3.4 |
| 3 | Clay_Hill | 0.0030 | 0.0008 | 3.5 |
| 3 | Cramor_Place | 0.0045 | 0.0008 | 5.4 |
| 3 | Nickel | 0.0046 | 0.0005 | 9.3 |

## Reading

(a) says how large each cluster's term would be from finite samples alone. (b) carries the uncertainty in the observed counts into the excess; a posterior that stays above 1 says the excess is not an accident of which sherds were collected. (c) says whether the excess belongs to the cluster or to one collection in it.
