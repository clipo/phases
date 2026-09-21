# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 1200 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0207; calibrated cell N 10000, innovation 0.0002, mixing 0.005.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | innovation boundary | boundary strength 0.0 | 0.0029 | [0.0007, 0.0136] | 7.2x | 0.8% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0047 | [0.0009, 0.0854] | 4.4x | 15.8% | matched |
| basin | innovation boundary | boundary strength 0.18 | 0.0055 | [0.0009, 0.0843] | 3.8x | 20.2% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0064 | [0.0010, 0.0851] | 3.2x | 24.5% | matched |
| basin | innovation boundary | boundary strength 0.22 | 0.0078 | [0.0010, 0.0867] | 2.6x | 26.9% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0107 | [0.0011, 0.0942] | 1.9x | 30.9% | matched |
| basin | innovation boundary | boundary strength 0.28 | 0.0120 | [0.0013, 0.0933] | 1.7x | 34.5% | broken |
| basin | innovation boundary | boundary strength 0.3 | 0.0145 | [0.0013, 0.0982] | 1.4x | 38.8% | broken |
| basin | innovation boundary | boundary strength 0.35 | 0.0190 | [0.0018, 0.1005] | 1.1x | 48.2% | broken |
| basin | innovation boundary | boundary strength 0.4 | 0.0267 | [0.0019, 0.0997] | 0.8x | 56.1% | broken |

## Reading


**basin**

- **innovation boundary**: best shortfall 0.8x at boundary strength 0.4; the best cell that keeps diversity matched is boundary strength 0.25 at 1.9x.
