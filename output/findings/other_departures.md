# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 250 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0062; calibrated cell N 2000, innovation 0.002, mixing 0.005.
- **cmv**: observed $F_{ST}$ 0.0337; calibrated cell N 120, innovation 0.004, mixing 0.02.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | transport geography | river, exponential, 6 km | 0.0032 | [0.0004, 0.0220] | 1.9x | 26.4% | broken |
| basin | transport geography | river, exponential, 12 km | 0.0034 | [0.0005, 0.0170] | 1.9x | 24.8% | broken |
| basin | transport geography | river, exponential, 24 km | 0.0025 | [0.0003, 0.0125] | 2.5x | 19.6% | matched |
| basin | transport geography | river, exponential, 48 km | 0.0019 | [0.0002, 0.0093] | 3.3x | 7.2% | matched |
| basin | transport geography | river, exponential, 96 km | 0.0015 | [0.0002, 0.0086] | 4.2x | 6.4% | matched |
| basin | transport geography | straight-line, exponential, 6 km | 0.0030 | [0.0003, 0.0174] | 2.1x | 20.4% | matched |
| basin | transport geography | straight-line, exponential, 12 km | 0.0020 | [0.0002, 0.0118] | 3.1x | 12.0% | matched |
| basin | transport geography | straight-line, exponential, 24 km | 0.0015 | [0.0002, 0.0088] | 4.1x | 8.4% | matched |
| basin | transport geography | straight-line, exponential, 48 km | 0.0013 | [0.0002, 0.0083] | 4.8x | 7.6% | matched |
| basin | transport geography | straight-line, exponential, 96 km | 0.0012 | [0.0001, 0.0076] | 5.1x | 5.2% | matched |
| basin | transport geography | river, gaussian, 6 km | 0.0041 | [0.0003, 0.0200] | 1.5x | 32.0% | broken |
| basin | transport geography | river, gaussian, 12 km | 0.0034 | [0.0004, 0.0193] | 1.9x | 31.6% | broken |
| basin | transport geography | river, gaussian, 24 km | 0.0033 | [0.0003, 0.0220] | 1.9x | 24.4% | broken |
| basin | transport geography | river, gaussian, 48 km | 0.0030 | [0.0003, 0.0115] | 2.1x | 20.8% | matched |
| basin | transport geography | river, gaussian, 96 km | 0.0017 | [0.0001, 0.0115] | 3.6x | 6.8% | matched |
| basin | transport geography | river, power, 6 km | 0.0027 | [0.0003, 0.0154] | 2.3x | 16.4% | matched |
| basin | transport geography | river, power, 12 km | 0.0021 | [0.0002, 0.0133] | 2.9x | 15.2% | matched |
| basin | transport geography | river, power, 24 km | 0.0020 | [0.0003, 0.0089] | 3.1x | 12.8% | matched |
| basin | transport geography | river, power, 48 km | 0.0017 | [0.0002, 0.0084] | 3.6x | 8.8% | matched |
| basin | transport geography | river, power, 96 km | 0.0014 | [0.0001, 0.0070] | 4.3x | 4.8% | matched |
| basin | accumulation spans | window 8 +/- 0 | 0.0025 | [0.0003, 0.0125] | 2.5x | 19.6% | matched |
| basin | accumulation spans | window 8 +/- 2 | 0.0023 | [0.0003, 0.0139] | 2.7x | 19.6% | matched |
| basin | accumulation spans | window 8 +/- 4 | 0.0026 | [0.0002, 0.0135] | 2.4x | 16.8% | matched |
| basin | accumulation spans | window 8 +/- 6 | 0.0024 | [0.0002, 0.0139] | 2.6x | 16.0% | matched |
| basin | innovation boundary | boundary strength 0.0 | 0.0025 | [0.0003, 0.0125] | 2.5x | 19.6% | matched |
| basin | innovation boundary | boundary strength 0.1 | 0.0033 | [0.0003, 0.1340] | 1.9x | 34.0% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0037 | [0.0004, 0.1309] | 1.7x | 35.6% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0054 | [0.0005, 0.1563] | 1.1x | 47.2% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0225 | [0.0004, 0.1580] | 0.3x | 60.4% | matched |
| basin | innovation boundary | boundary strength 0.3 | 0.0491 | [0.0006, 0.1653] | 0.1x | 68.0% | matched |
| basin | innovation boundary | boundary strength 0.4 | 0.0645 | [0.0012, 0.1717] | 0.1x | 80.0% | broken |
| basin | innovation boundary | boundary strength 0.6 | 0.1209 | [0.0053, 0.1835] | 0.1x | 96.4% | broken |
| basin | innovation boundary | boundary strength 1.0 | 0.1340 | [0.0437, 0.1836] | 0.0x | 99.6% | broken |
| cmv | transport geography | river, exponential, 6 km | 0.0260 | [0.0079, 0.0726] | 1.3x | 34.4% | broken |
| cmv | transport geography | river, exponential, 12 km | 0.0209 | [0.0057, 0.0742] | 1.6x | 18.0% | broken |
| cmv | transport geography | river, exponential, 24 km | 0.0128 | [0.0045, 0.0324] | 2.6x | 2.4% | matched |
| cmv | transport geography | river, exponential, 48 km | 0.0114 | [0.0036, 0.0321] | 3.0x | 2.4% | matched |
| cmv | transport geography | river, exponential, 96 km | 0.0106 | [0.0033, 0.0285] | 3.2x | 0.4% | matched |
| cmv | transport geography | straight-line, exponential, 6 km | 0.0260 | [0.0079, 0.0726] | 1.3x | 34.4% | broken |
| cmv | transport geography | straight-line, exponential, 12 km | 0.0209 | [0.0057, 0.0742] | 1.6x | 18.0% | broken |
| cmv | transport geography | straight-line, exponential, 24 km | 0.0128 | [0.0045, 0.0324] | 2.6x | 2.4% | matched |
| cmv | transport geography | straight-line, exponential, 48 km | 0.0114 | [0.0036, 0.0321] | 3.0x | 2.4% | matched |
| cmv | transport geography | straight-line, exponential, 96 km | 0.0106 | [0.0033, 0.0285] | 3.2x | 0.4% | matched |
| cmv | transport geography | river, gaussian, 6 km | 0.0290 | [0.0075, 0.0910] | 1.2x | 44.0% | broken |
| cmv | transport geography | river, gaussian, 12 km | 0.0265 | [0.0087, 0.0887] | 1.3x | 34.8% | broken |
| cmv | transport geography | river, gaussian, 24 km | 0.0186 | [0.0054, 0.0534] | 1.8x | 17.2% | broken |
| cmv | transport geography | river, gaussian, 48 km | 0.0120 | [0.0041, 0.0376] | 2.8x | 3.2% | matched |
| cmv | transport geography | river, gaussian, 96 km | 0.0101 | [0.0039, 0.0272] | 3.3x | 0.8% | matched |
| cmv | transport geography | river, power, 6 km | 0.0155 | [0.0048, 0.0445] | 2.2x | 9.6% | broken |
| cmv | transport geography | river, power, 12 km | 0.0129 | [0.0043, 0.0372] | 2.6x | 4.4% | matched |
| cmv | transport geography | river, power, 24 km | 0.0114 | [0.0037, 0.0322] | 3.0x | 2.0% | matched |
| cmv | transport geography | river, power, 48 km | 0.0111 | [0.0031, 0.0306] | 3.0x | 2.4% | matched |
| cmv | transport geography | river, power, 96 km | 0.0107 | [0.0032, 0.0247] | 3.2x | 0.0% | matched |
| cmv | accumulation spans | window 8 +/- 0 | 0.0128 | [0.0045, 0.0324] | 2.6x | 2.4% | matched |
| cmv | accumulation spans | window 8 +/- 2 | 0.0118 | [0.0041, 0.0324] | 2.9x | 2.4% | matched |
| cmv | accumulation spans | window 8 +/- 4 | 0.0121 | [0.0042, 0.0355] | 2.8x | 3.6% | matched |
| cmv | accumulation spans | window 8 +/- 6 | 0.0121 | [0.0039, 0.0375] | 2.8x | 3.6% | matched |
| cmv | innovation boundary | boundary strength 0.0 | 0.0128 | [0.0045, 0.0324] | 2.6x | 2.4% | matched |
| cmv | innovation boundary | boundary strength 0.1 | 0.0164 | [0.0041, 0.1080] | 2.1x | 27.2% | matched |
| cmv | innovation boundary | boundary strength 0.15 | 0.0239 | [0.0049, 0.1270] | 1.4x | 40.8% | broken |
| cmv | innovation boundary | boundary strength 0.2 | 0.0338 | [0.0066, 0.1270] | 1.0x | 50.0% | broken |
| cmv | innovation boundary | boundary strength 0.25 | 0.0412 | [0.0063, 0.1363] | 0.8x | 56.8% | broken |
| cmv | innovation boundary | boundary strength 0.3 | 0.0549 | [0.0061, 0.1444] | 0.6x | 67.2% | broken |
| cmv | innovation boundary | boundary strength 0.4 | 0.0687 | [0.0084, 0.1391] | 0.5x | 81.6% | broken |
| cmv | innovation boundary | boundary strength 0.6 | 0.0820 | [0.0227, 0.1438] | 0.4x | 92.8% | broken |
| cmv | innovation boundary | boundary strength 1.0 | 0.0860 | [0.0481, 0.1346] | 0.4x | 99.6% | broken |

## Reading


**basin**

- **transport geography**: best shortfall 1.5x at river, gaussian, 6 km; the best cell that keeps diversity matched is straight-line, exponential, 6 km at 2.1x.
- **accumulation spans**: best shortfall 2.4x at window 8 +/- 4; the best cell that keeps diversity matched is window 8 +/- 4 at 2.4x.
- **innovation boundary**: best shortfall 0.0x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.3 at 0.1x.

**cmv**

- **transport geography**: best shortfall 1.2x at river, gaussian, 6 km; the best cell that keeps diversity matched is river, power, 12 km at 2.6x.
- **accumulation spans**: best shortfall 2.6x at window 8 +/- 0; the best cell that keeps diversity matched is window 8 +/- 0 at 2.6x.
- **innovation boundary**: best shortfall 0.4x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.1 at 2.1x.
