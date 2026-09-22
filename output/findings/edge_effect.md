# Is the drift shortfall an edge effect of a closed basin?

Basin phase set, 28 assemblages scored; the open arm simulates 38 (10 outside the set).
k = 5, 300 realizations per arm on shared seeds, calibrated pooled cell (2000 learners, innovation 0.002, mixing 0.005), interaction length 24 river-km.

| arm | observed F_ST | drift median | shortfall | H_S | richness | H_T |
|---|---|---|---|---|---|---|
| observed | 0.0307 | | | 0.563 | 6.25 | 0.592 |
| closed | 0.0307 | 0.0072 | **4.3x** | 0.570 | 5.73 | 0.588 |
| open | 0.0307 | 0.0067 | **4.6x** | 0.572 | 5.79 | 0.590 |
| open, distinct outside | 0.0307 | 0.0073 | **4.2x** | 0.593 | 6.00 | 0.609 |

## Per cluster

| cluster | n | outside assemblages within one interaction length | observed term | excess, closed | excess, open | excess, open with distinct outside |
|---|---|---|---|---|---|---|
| 0 (Barton_Ranch, Cummins, Fortune...) | 10 | 0 | 0.0052 | 3.8x | 3.8x | **3.3x** |
| 1 (Clay_Hill, Davis, Grant...) | 5 | 1 | 0.0150 | 32.4x | 33.0x | **29.7x** |
| 2 (Irby, Lake_Cormorant, Mound_Place...) | 5 | 2 | 0.0044 | 7.4x | 7.3x | **7.6x** |
| 3 (Big_Eddy, Castile_Landing, Cramor_Place...) | 4 | 2 | 0.0028 | 1.9x | 2.0x | **1.9x** |
| 4 (Beck, Belle_Meade, Commerce...) | 4 | 1 | 0.0034 | 2.1x | 2.3x | **2.2x** |

The outside assemblages' pooled profile differs from the basin's by an L1 distance of 0.596
(0 identical, 2 disjoint).

A cluster with no outside assemblage within one interaction length has an edge this test
cannot open; its row says nothing about whether an edge effect operates there.
