# Is the drift shortfall an edge effect of a closed basin?

Basin phase set, 28 assemblages scored; the open arm simulates 38 (10 outside the set).
k = 5, 300 realizations per arm on shared seeds, calibrated pooled cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 river-km.

| arm | observed F_ST | drift median | shortfall | H_S | richness | H_T |
|---|---|---|---|---|---|---|
| observed | 0.0307 | | | 0.563 | 6.25 | 0.592 |
| closed | 0.0307 | 0.0052 | **5.9x** | 0.578 | 6.00 | 0.588 |
| open | 0.0307 | 0.0042 | **7.3x** | 0.580 | 6.11 | 0.589 |
| open, distinct outside | 0.0307 | 0.0051 | **6.0x** | 0.612 | 6.39 | 0.617 |

## Per cluster

| cluster | n | outside assemblages within one interaction length | observed term | excess, closed | excess, open | excess, open with distinct outside |
|---|---|---|---|---|---|---|
| 0 (Barton_Ranch, Cummins, Fortune...) | 10 | 0 | 0.0052 | 3.9x | 4.6x | **3.8x** |
| 1 (Clay_Hill, Davis, Grant...) | 5 | 1 | 0.0150 | 39.7x | 37.7x | **34.1x** |
| 2 (Irby, Lake_Cormorant, Mound_Place...) | 5 | 2 | 0.0044 | 8.6x | 11.5x | **9.0x** |
| 3 (Big_Eddy, Castile_Landing, Cramor_Place...) | 4 | 2 | 0.0028 | 3.0x | 3.5x | **3.5x** |
| 4 (Beck, Belle_Meade, Commerce...) | 4 | 1 | 0.0034 | 3.1x | 4.5x | **3.7x** |

The outside assemblages' pooled profile differs from the basin's by an L1 distance of 0.596
(0 identical, 2 disjoint).

A cluster with no outside assemblage within one interaction length has an edge this test
cannot open; its row says nothing about whether an edge effect operates there.
