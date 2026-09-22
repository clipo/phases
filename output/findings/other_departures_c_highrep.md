# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 1200 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0062; calibrated cell N 120, innovation 0.024, mixing 0.01.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | innovation boundary | boundary strength 0.0 | 0.0016 | [0.0002, 0.0096] | 3.8x | 8.8% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0031 | [0.0003, 0.1924] | 2.0x | 33.9% | matched |
| basin | innovation boundary | boundary strength 0.18 | 0.0042 | [0.0003, 0.2022] | 1.5x | 41.6% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0054 | [0.0003, 0.2068] | 1.2x | 47.2% | matched |
| basin | innovation boundary | boundary strength 0.22 | 0.0071 | [0.0004, 0.2085] | 0.9x | 52.4% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0106 | [0.0004, 0.2245] | 0.6x | 59.2% | matched |
| basin | innovation boundary | boundary strength 0.28 | 0.0527 | [0.0004, 0.2322] | 0.1x | 63.8% | matched |
| basin | innovation boundary | boundary strength 0.3 | 0.0608 | [0.0004, 0.2292] | 0.1x | 68.2% | matched |
| basin | innovation boundary | boundary strength 0.35 | 0.0727 | [0.0006, 0.2330] | 0.1x | 76.8% | matched |
| basin | innovation boundary | boundary strength 0.4 | 0.0861 | [0.0007, 0.2333] | 0.1x | 81.0% | broken |

## Reading


**basin**

- **innovation boundary**: best shortfall 0.1x at boundary strength 0.4; the best cell that keeps diversity matched is boundary strength 0.35 at 0.1x.
