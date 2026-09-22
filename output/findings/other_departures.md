# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 250 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0062; calibrated cell N 120, innovation 0.024, mixing 0.01.
- **cmv**: observed $F_{ST}$ 0.0337; calibrated cell N 120, innovation 0.012, mixing 0.005.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | transport geography | river, exponential, 6 km | 0.0017 | [0.0002, 0.0101] | 3.7x | 9.2% | matched |
| basin | transport geography | river, exponential, 12 km | 0.0018 | [0.0003, 0.0096] | 3.4x | 8.8% | matched |
| basin | transport geography | river, exponential, 24 km | 0.0015 | [0.0002, 0.0116] | 4.0x | 10.0% | matched |
| basin | transport geography | river, exponential, 48 km | 0.0015 | [0.0002, 0.0108] | 4.2x | 8.4% | matched |
| basin | transport geography | river, exponential, 96 km | 0.0018 | [0.0002, 0.0095] | 3.5x | 11.2% | matched |
| basin | transport geography | straight-line, exponential, 6 km | 0.0017 | [0.0002, 0.0093] | 3.6x | 8.4% | matched |
| basin | transport geography | straight-line, exponential, 12 km | 0.0016 | [0.0001, 0.0089] | 4.0x | 5.6% | matched |
| basin | transport geography | straight-line, exponential, 24 km | 0.0016 | [0.0002, 0.0103] | 3.8x | 6.4% | matched |
| basin | transport geography | straight-line, exponential, 48 km | 0.0019 | [0.0001, 0.0101] | 3.3x | 10.0% | matched |
| basin | transport geography | straight-line, exponential, 96 km | 0.0017 | [0.0002, 0.0078] | 3.8x | 4.8% | matched |
| basin | transport geography | river, gaussian, 6 km | 0.0020 | [0.0003, 0.0106] | 3.2x | 11.6% | broken |
| basin | transport geography | river, gaussian, 12 km | 0.0016 | [0.0002, 0.0103] | 3.8x | 12.0% | broken |
| basin | transport geography | river, gaussian, 24 km | 0.0018 | [0.0002, 0.0107] | 3.4x | 10.4% | matched |
| basin | transport geography | river, gaussian, 48 km | 0.0015 | [0.0002, 0.0086] | 4.2x | 6.8% | matched |
| basin | transport geography | river, gaussian, 96 km | 0.0016 | [0.0002, 0.0081] | 3.8x | 5.6% | matched |
| basin | transport geography | river, power, 6 km | 0.0017 | [0.0002, 0.0112] | 3.6x | 12.4% | matched |
| basin | transport geography | river, power, 12 km | 0.0020 | [0.0003, 0.0097] | 3.1x | 8.4% | matched |
| basin | transport geography | river, power, 24 km | 0.0018 | [0.0002, 0.0102] | 3.5x | 9.2% | matched |
| basin | transport geography | river, power, 48 km | 0.0020 | [0.0003, 0.0097] | 3.2x | 9.6% | matched |
| basin | transport geography | river, power, 96 km | 0.0015 | [0.0002, 0.0087] | 4.1x | 7.2% | matched |
| basin | accumulation spans | window 8 +/- 0 | 0.0015 | [0.0002, 0.0116] | 4.0x | 10.0% | matched |
| basin | accumulation spans | window 8 +/- 2 | 0.0016 | [0.0002, 0.0090] | 3.9x | 8.4% | matched |
| basin | accumulation spans | window 8 +/- 4 | 0.0016 | [0.0002, 0.0095] | 4.0x | 7.6% | matched |
| basin | accumulation spans | window 8 +/- 6 | 0.0017 | [0.0003, 0.0098] | 3.7x | 7.6% | broken |
| basin | innovation boundary | boundary strength 0.0 | 0.0015 | [0.0002, 0.0116] | 4.0x | 10.0% | matched |
| basin | innovation boundary | boundary strength 0.1 | 0.0024 | [0.0003, 0.1824] | 2.7x | 25.6% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0034 | [0.0003, 0.1876] | 1.8x | 35.2% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0063 | [0.0004, 0.2090] | 1.0x | 50.0% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0227 | [0.0005, 0.2247] | 0.3x | 61.6% | matched |
| basin | innovation boundary | boundary strength 0.3 | 0.0641 | [0.0004, 0.2289] | 0.1x | 69.2% | matched |
| basin | innovation boundary | boundary strength 0.4 | 0.0852 | [0.0011, 0.2350] | 0.1x | 83.6% | broken |
| basin | innovation boundary | boundary strength 0.6 | 0.1669 | [0.0070, 0.2360] | 0.0x | 98.0% | broken |
| basin | innovation boundary | boundary strength 1.0 | 0.1984 | [0.0565, 0.2469] | 0.0x | 100.0% | broken |
| cmv | transport geography | river, exponential, 6 km | 0.0124 | [0.0041, 0.0385] | 2.7x | 6.0% | matched |
| cmv | transport geography | river, exponential, 12 km | 0.0127 | [0.0045, 0.0351] | 2.7x | 2.8% | matched |
| cmv | transport geography | river, exponential, 24 km | 0.0118 | [0.0037, 0.0356] | 2.9x | 3.6% | matched |
| cmv | transport geography | river, exponential, 48 km | 0.0112 | [0.0033, 0.0352] | 3.0x | 3.2% | matched |
| cmv | transport geography | river, exponential, 96 km | 0.0114 | [0.0039, 0.0377] | 3.0x | 3.6% | matched |
| cmv | transport geography | straight-line, exponential, 6 km | 0.0124 | [0.0041, 0.0385] | 2.7x | 6.0% | matched |
| cmv | transport geography | straight-line, exponential, 12 km | 0.0127 | [0.0045, 0.0351] | 2.7x | 2.8% | matched |
| cmv | transport geography | straight-line, exponential, 24 km | 0.0118 | [0.0037, 0.0356] | 2.9x | 3.6% | matched |
| cmv | transport geography | straight-line, exponential, 48 km | 0.0112 | [0.0033, 0.0352] | 3.0x | 3.2% | matched |
| cmv | transport geography | straight-line, exponential, 96 km | 0.0114 | [0.0039, 0.0377] | 3.0x | 3.6% | matched |
| cmv | transport geography | river, gaussian, 6 km | 0.0132 | [0.0043, 0.0361] | 2.6x | 3.6% | matched |
| cmv | transport geography | river, gaussian, 12 km | 0.0128 | [0.0039, 0.0336] | 2.6x | 2.8% | matched |
| cmv | transport geography | river, gaussian, 24 km | 0.0123 | [0.0036, 0.0422] | 2.8x | 5.2% | matched |
| cmv | transport geography | river, gaussian, 48 km | 0.0118 | [0.0033, 0.0347] | 2.8x | 2.8% | matched |
| cmv | transport geography | river, gaussian, 96 km | 0.0109 | [0.0036, 0.0318] | 3.1x | 1.6% | matched |
| cmv | transport geography | river, power, 6 km | 0.0122 | [0.0040, 0.0336] | 2.8x | 2.4% | matched |
| cmv | transport geography | river, power, 12 km | 0.0123 | [0.0046, 0.0316] | 2.7x | 2.4% | matched |
| cmv | transport geography | river, power, 24 km | 0.0124 | [0.0031, 0.0341] | 2.7x | 2.8% | matched |
| cmv | transport geography | river, power, 48 km | 0.0114 | [0.0039, 0.0296] | 3.0x | 0.4% | matched |
| cmv | transport geography | river, power, 96 km | 0.0111 | [0.0037, 0.0316] | 3.1x | 2.0% | matched |
| cmv | accumulation spans | window 8 +/- 0 | 0.0118 | [0.0037, 0.0356] | 2.9x | 3.6% | matched |
| cmv | accumulation spans | window 8 +/- 2 | 0.0116 | [0.0041, 0.0348] | 2.9x | 2.8% | matched |
| cmv | accumulation spans | window 8 +/- 4 | 0.0113 | [0.0036, 0.0355] | 3.0x | 2.8% | matched |
| cmv | accumulation spans | window 8 +/- 6 | 0.0112 | [0.0039, 0.0335] | 3.0x | 2.0% | matched |
| cmv | innovation boundary | boundary strength 0.0 | 0.0118 | [0.0037, 0.0356] | 2.9x | 3.6% | matched |
| cmv | innovation boundary | boundary strength 0.1 | 0.0186 | [0.0050, 0.4261] | 1.8x | 32.4% | matched |
| cmv | innovation boundary | boundary strength 0.15 | 0.0304 | [0.0068, 0.4340] | 1.1x | 48.4% | broken |
| cmv | innovation boundary | boundary strength 0.2 | 0.2501 | [0.0072, 0.4677] | 0.1x | 58.8% | broken |
| cmv | innovation boundary | boundary strength 0.25 | 0.2779 | [0.0075, 0.4627] | 0.1x | 65.2% | broken |
| cmv | innovation boundary | boundary strength 0.3 | 0.3161 | [0.0107, 0.4714] | 0.1x | 76.8% | broken |
| cmv | innovation boundary | boundary strength 0.4 | 0.3481 | [0.0125, 0.4805] | 0.1x | 88.4% | broken |
| cmv | innovation boundary | boundary strength 0.6 | 0.4025 | [0.1883, 0.5134] | 0.1x | 99.6% | broken |
| cmv | innovation boundary | boundary strength 1.0 | 0.4250 | [0.3321, 0.5163] | 0.1x | 100.0% | broken |

## Reading


**basin**

- **transport geography**: best shortfall 3.1x at river, power, 12 km; the best cell that keeps diversity matched is river, power, 12 km at 3.1x.
- **accumulation spans**: best shortfall 3.7x at window 8 +/- 6; the best cell that keeps diversity matched is window 8 +/- 2 at 3.9x.
- **innovation boundary**: best shortfall 0.0x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.3 at 0.1x.

**cmv**

- **transport geography**: best shortfall 2.6x at river, gaussian, 6 km; the best cell that keeps diversity matched is river, gaussian, 6 km at 2.6x.
- **accumulation spans**: best shortfall 2.9x at window 8 +/- 0; the best cell that keeps diversity matched is window 8 +/- 0 at 2.9x.
- **innovation boundary**: best shortfall 0.1x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.1 at 1.8x.
