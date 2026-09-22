# Is the drift shortfall an edge effect of a closed basin?

Basin phase set, 28 assemblages scored; the open arm simulates 38 (10 outside the set).
k = 5, 300 realizations per arm on shared seeds, calibrated pooled cell (120 learners, innovation 0.024, mixing 0.01), interaction length 24 river-km.

| arm | observed F_ST | drift median | shortfall | H_S | richness | H_T |
|---|---|---|---|---|---|---|
| observed | 0.0307 | | | 0.563 | 6.25 | 0.592 |
| closed | 0.0307 | 0.0062 | **5.0x** | 0.565 | 5.64 | 0.590 |
| open | 0.0307 | 0.0068 | **4.5x** | 0.569 | 5.71 | 0.591 |
| open, distinct outside | 0.0307 | 0.0064 | **4.8x** | 0.579 | 5.75 | 0.600 |

## Per cluster

| cluster | n | outside assemblages within one interaction length | observed term | excess, closed | excess, open | excess, open with distinct outside |
|---|---|---|---|---|---|---|
| 0 (Barton_Ranch, Cummins, Fortune...) | 10 | 0 | 0.0052 | 5.7x | 4.6x | **5.2x** |
| 1 (Clay_Hill, Davis, Grant...) | 5 | 1 | 0.0150 | 37.9x | 40.9x | **38.1x** |
| 2 (Irby, Lake_Cormorant, Mound_Place...) | 5 | 2 | 0.0044 | 7.8x | 9.1x | **8.2x** |
| 3 (Big_Eddy, Castile_Landing, Cramor_Place...) | 4 | 2 | 0.0028 | 1.9x | 1.8x | **1.7x** |
| 4 (Beck, Belle_Meade, Commerce...) | 4 | 1 | 0.0034 | 2.2x | 1.8x | **2.2x** |

The outside assemblages' pooled profile differs from the basin's by an L1 distance of 0.596
(0 identical, 2 disjoint).

A cluster with no outside assemblage within one interaction length has an edge this test
cannot open; its row says nothing about whether an edge effect operates there.
