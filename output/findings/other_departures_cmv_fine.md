# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 250 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **cmv**: observed $F_{ST}$ 0.0337; calibrated cell N 120, innovation 0.012, mixing 0.005.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| cmv | innovation boundary | boundary strength 0.1 | 0.0186 | [0.0050, 0.4261] | 1.8x | 32.4% | matched |
| cmv | innovation boundary | boundary strength 0.11 | 0.0199 | [0.0049, 0.4333] | 1.7x | 35.6% | broken |
| cmv | innovation boundary | boundary strength 0.12 | 0.0194 | [0.0051, 0.4146] | 1.7x | 38.8% | matched |
| cmv | innovation boundary | boundary strength 0.13 | 0.0222 | [0.0054, 0.4329] | 1.5x | 42.4% | broken |
| cmv | innovation boundary | boundary strength 0.14 | 0.0252 | [0.0054, 0.4357] | 1.3x | 43.2% | broken |
| cmv | innovation boundary | boundary strength 0.15 | 0.0304 | [0.0068, 0.4340] | 1.1x | 48.4% | broken |

## Reading


**cmv**

- **innovation boundary**: best shortfall 1.1x at boundary strength 0.15; the best cell that keeps diversity matched is boundary strength 0.12 at 1.7x.
