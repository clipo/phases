# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 1200 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0062; calibrated cell N 2000, innovation 0.002, mixing 0.005.
- **cmv**: observed $F_{ST}$ 0.0337; calibrated cell N 120, innovation 0.004, mixing 0.02.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | innovation boundary | boundary strength 0.0 | 0.0023 | [0.0003, 0.0136] | 2.7x | 16.4% | matched |
| basin | innovation boundary | boundary strength 0.1 | 0.0032 | [0.0003, 0.1331] | 2.0x | 29.5% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0040 | [0.0004, 0.1414] | 1.6x | 37.3% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0059 | [0.0004, 0.1558] | 1.1x | 48.9% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0125 | [0.0004, 0.1618] | 0.5x | 59.1% | matched |
| basin | innovation boundary | boundary strength 0.3 | 0.0447 | [0.0006, 0.1672] | 0.1x | 67.2% | matched |
| basin | innovation boundary | boundary strength 0.4 | 0.0614 | [0.0009, 0.1777] | 0.1x | 79.8% | broken |
| basin | innovation boundary | boundary strength 0.6 | 0.1194 | [0.0033, 0.1806] | 0.1x | 94.8% | broken |
| basin | innovation boundary | boundary strength 1.0 | 0.1386 | [0.0361, 0.1896] | 0.0x | 99.4% | broken |
| cmv | innovation boundary | boundary strength 0.0 | 0.0124 | [0.0041, 0.0362] | 2.7x | 3.2% | matched |
| cmv | innovation boundary | boundary strength 0.1 | 0.0166 | [0.0047, 0.1176] | 2.0x | 26.9% | matched |
| cmv | innovation boundary | boundary strength 0.15 | 0.0216 | [0.0048, 0.1256] | 1.6x | 37.5% | broken |
| cmv | innovation boundary | boundary strength 0.2 | 0.0329 | [0.0058, 0.1321] | 1.0x | 49.5% | broken |
| cmv | innovation boundary | boundary strength 0.25 | 0.0452 | [0.0064, 0.1406] | 0.7x | 57.6% | broken |
| cmv | innovation boundary | boundary strength 0.3 | 0.0545 | [0.0063, 0.1344] | 0.6x | 66.2% | broken |
| cmv | innovation boundary | boundary strength 0.4 | 0.0695 | [0.0090, 0.1379] | 0.5x | 81.2% | broken |
| cmv | innovation boundary | boundary strength 0.6 | 0.0828 | [0.0202, 0.1427] | 0.4x | 93.6% | broken |
| cmv | innovation boundary | boundary strength 1.0 | 0.0866 | [0.0480, 0.1391] | 0.4x | 99.8% | broken |

## Reading


**basin**

- **innovation boundary**: best shortfall 0.0x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.3 at 0.1x.

**cmv**

- **innovation boundary**: best shortfall 0.4x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.1 at 2.0x.

Run through `analyses/66_innovation_boundary_highrep.py` (65 --only C --reps 1200 --tag c_highrep) on 2026-09-22; the JSON beside it is `output/other_departures_c_highrep.json`.
