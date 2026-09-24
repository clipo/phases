# Is the drift shortfall an edge effect of a closed basin?

Basin phase set, 28 assemblages scored; the open arm simulates 38 (10 outside the set).
k = 5, 300 realizations per arm on shared seeds, calibrated pooled cell (2000 learners, innovation 0.001, mixing 0.02), interaction length 24 river-km.

| arm | observed F_ST | drift median | shortfall | H_S | richness | H_T |
|---|---|---|---|---|---|---|
| observed | 0.0206 | | | 0.563 | 6.25 | 0.592 |
| closed | 0.0206 | 0.0040 | **5.2x** | 0.575 | 6.04 | 0.587 |
| open | 0.0206 | 0.0037 | **5.6x** | 0.577 | 6.11 | 0.590 |
| open, distinct outside | 0.0206 | 0.0044 | **4.7x** | 0.615 | 6.46 | 0.622 |

## Per cluster

| cluster | n | outside assemblages within one interaction length | observed term | excess, closed | excess, open | excess, open with distinct outside |
|---|---|---|---|---|---|---|
| 0 (Barton_Ranch, Fortune, Holden_Lake...) | 9 | 0 | 0.0045 | 4.6x | 4.3x | **3.5x** |
| 1 (Irby, Lake_Cormorant, Mound_Place...) | 5 | 2 | 0.0044 | 8.9x | 10.6x | **8.9x** |
| 2 (Cummins...) | 1 | 0 | 0.0026 | 15.5x | 17.2x | **14.7x** |
| 3 (Big_Eddy, Castile_Landing, Clay_Hill...) | 5 | 2 | 0.0035 | 4.0x | 4.3x | **4.1x** |
| 4 (Beck, Belle_Meade, Commerce...) | 8 | 1 | 0.0055 | 8.3x | 7.1x | **6.9x** |

The outside assemblages' pooled profile differs from the basin's by an L1 distance of 0.596
(0 identical, 2 disjoint).

A cluster with no outside assemblage within one interaction length has an edge this test
cannot open; its row says nothing about whether an edge effect operates there.
