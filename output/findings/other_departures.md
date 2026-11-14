# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 250 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0062; calibrated cell N 2000, innovation 0.001, mixing 0.02.
- **cmv**: observed $F_{ST}$ 0.0337; calibrated cell N 120, innovation 0.004, mixing 0.02.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | transport geography | river, exponential, 6 km | 0.0054 | [0.0006, 0.0330] | 1.2x | 44.0% | broken |
| basin | transport geography | river, exponential, 12 km | 0.0047 | [0.0004, 0.0244] | 1.3x | 39.6% | broken |
| basin | transport geography | river, exponential, 24 km | 0.0018 | [0.0002, 0.0114] | 3.5x | 12.8% | matched |
| basin | transport geography | river, exponential, 48 km | 0.0008 | [0.0001, 0.0051] | 7.7x | 0.4% | matched |
| basin | transport geography | river, exponential, 96 km | 0.0007 | [0.0001, 0.0040] | 9.5x | 0.8% | matched |
| basin | transport geography | straight-line, exponential, 6 km | 0.0037 | [0.0004, 0.0195] | 1.7x | 26.8% | matched |
| basin | transport geography | straight-line, exponential, 12 km | 0.0010 | [0.0001, 0.0054] | 5.9x | 2.0% | matched |
| basin | transport geography | straight-line, exponential, 24 km | 0.0006 | [0.0001, 0.0034] | 10.5x | 0.0% | matched |
| basin | transport geography | straight-line, exponential, 48 km | 0.0005 | [0.0001, 0.0035] | 11.5x | 0.0% | matched |
| basin | transport geography | straight-line, exponential, 96 km | 0.0005 | [0.0001, 0.0033] | 11.5x | 0.0% | matched |
| basin | transport geography | river, gaussian, 6 km | 0.0099 | [0.0011, 0.0547] | 0.6x | 71.2% | broken |
| basin | transport geography | river, gaussian, 12 km | 0.0063 | [0.0004, 0.0371] | 1.0x | 50.0% | broken |
| basin | transport geography | river, gaussian, 24 km | 0.0044 | [0.0005, 0.0323] | 1.4x | 36.8% | broken |
| basin | transport geography | river, gaussian, 48 km | 0.0024 | [0.0002, 0.0128] | 2.6x | 17.6% | matched |
| basin | transport geography | river, gaussian, 96 km | 0.0007 | [0.0001, 0.0036] | 9.6x | 0.4% | matched |
| basin | transport geography | river, power, 6 km | 0.0020 | [0.0002, 0.0097] | 3.2x | 8.0% | matched |
| basin | transport geography | river, power, 12 km | 0.0010 | [0.0001, 0.0055] | 6.0x | 0.8% | matched |
| basin | transport geography | river, power, 24 km | 0.0008 | [0.0001, 0.0035] | 8.1x | 0.4% | matched |
| basin | transport geography | river, power, 48 km | 0.0006 | [0.0001, 0.0040] | 9.9x | 0.4% | matched |
| basin | transport geography | river, power, 96 km | 0.0005 | [0.0001, 0.0030] | 12.6x | 0.0% | matched |
| basin | accumulation spans | window 8 +/- 0 | 0.0018 | [0.0002, 0.0114] | 3.5x | 12.8% | matched |
| basin | accumulation spans | window 8 +/- 2 | 0.0019 | [0.0002, 0.0114] | 3.3x | 13.6% | matched |
| basin | accumulation spans | window 8 +/- 4 | 0.0018 | [0.0002, 0.0119] | 3.5x | 12.8% | matched |
| basin | accumulation spans | window 8 +/- 6 | 0.0018 | [0.0002, 0.0113] | 3.4x | 13.2% | matched |
| basin | innovation boundary | boundary strength 0.0 | 0.0018 | [0.0002, 0.0114] | 3.5x | 12.8% | matched |
| basin | innovation boundary | boundary strength 0.1 | 0.0028 | [0.0003, 0.0424] | 2.2x | 29.2% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0034 | [0.0003, 0.0475] | 1.8x | 34.0% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0053 | [0.0003, 0.0562] | 1.2x | 47.2% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0097 | [0.0004, 0.0609] | 0.6x | 57.6% | matched |
| basin | innovation boundary | boundary strength 0.3 | 0.0154 | [0.0006, 0.0621] | 0.4x | 64.8% | matched |
| basin | innovation boundary | boundary strength 0.4 | 0.0219 | [0.0008, 0.0648] | 0.3x | 79.6% | broken |
| basin | innovation boundary | boundary strength 0.6 | 0.0366 | [0.0021, 0.0730] | 0.2x | 92.8% | broken |
| basin | innovation boundary | boundary strength 1.0 | 0.0414 | [0.0110, 0.0721] | 0.2x | 100.0% | broken |
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

- **transport geography**: best shortfall 0.6x at river, gaussian, 6 km; the best cell that keeps diversity matched is straight-line, exponential, 6 km at 1.7x.
- **accumulation spans**: best shortfall 3.3x at window 8 +/- 2; the best cell that keeps diversity matched is window 8 +/- 2 at 3.3x.
- **innovation boundary**: best shortfall 0.2x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.3 at 0.4x.

**cmv**

- **transport geography**: best shortfall 1.2x at river, gaussian, 6 km; the best cell that keeps diversity matched is river, power, 12 km at 2.6x.
- **accumulation spans**: best shortfall 2.6x at window 8 +/- 0; the best cell that keeps diversity matched is window 8 +/- 0 at 2.6x.
- **innovation boundary**: best shortfall 0.4x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.1 at 2.1x.
