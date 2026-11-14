# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 250 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0179; calibrated cell N 10000, innovation 0.0002, mixing 0.005.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | innovation boundary | boundary strength 0.0 | 0.0031 | [0.0007, 0.0126] | 5.8x | 0.8% | matched |
| basin | innovation boundary | boundary strength 0.1 | 0.0040 | [0.0007, 0.0572] | 4.5x | 12.0% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0054 | [0.0009, 0.0705] | 3.4x | 20.0% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0094 | [0.0008, 0.0823] | 1.9x | 35.2% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0164 | [0.0011, 0.0805] | 1.1x | 46.0% | matched |
| basin | innovation boundary | boundary strength 0.3 | 0.0194 | [0.0011, 0.0879] | 0.9x | 52.0% | matched |
| basin | innovation boundary | boundary strength 0.4 | 0.0286 | [0.0021, 0.0893] | 0.6x | 68.8% | broken |
| basin | innovation boundary | boundary strength 0.6 | 0.0529 | [0.0072, 0.0909] | 0.3x | 87.6% | broken |
| basin | innovation boundary | boundary strength 1.0 | 0.0610 | [0.0195, 0.1023] | 0.3x | 98.0% | broken |

## Reading


**basin**

- **innovation boundary**: best shortfall 0.3x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.3 at 0.9x.
