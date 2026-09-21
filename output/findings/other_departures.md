# The three remaining departures from the drift null

Produced by `analyses/65_other_departures.py`, 250 realizations per cell, each region at its own calibrated cell.

The St. Francis basin is the test; the southeast-Missouri (`cmv`) set is a comparison, not a replication, because those deposits are on the whole earlier and were collected under a different sampling regime.

- **basin**: observed $F_{ST}$ 0.0207; calibrated cell N 10000, innovation 0.0002, mixing 0.005.
- **cmv**: observed $F_{ST}$ 0.0337; calibrated cell N 120, innovation 0.012, mixing 0.005.

A cell closes the gap only if it reaches its region's observed value AND stays inside the calibration's 10 percent diversity tolerance.

| region | departure | setting | median $F_{ST}$ | 95% | shortfall | reaches obs | diversity |
|---|---|---|---:|---|---:|---:|:---:|
| basin | transport geography | river, exponential, 6 km | 0.0070 | [0.0013, 0.0287] | 2.9x | 8.4% | broken |
| basin | transport geography | river, exponential, 12 km | 0.0051 | [0.0012, 0.0252] | 4.0x | 6.4% | matched |
| basin | transport geography | river, exponential, 24 km | 0.0029 | [0.0009, 0.0135] | 7.0x | 1.2% | matched |
| basin | transport geography | river, exponential, 48 km | 0.0013 | [0.0003, 0.0046] | 16.5x | 0.0% | matched |
| basin | transport geography | river, exponential, 96 km | 0.0009 | [0.0003, 0.0036] | 22.0x | 0.0% | matched |
| basin | transport geography | straight-line, exponential, 6 km | 0.0048 | [0.0011, 0.0215] | 4.3x | 2.8% | matched |
| basin | transport geography | straight-line, exponential, 12 km | 0.0018 | [0.0005, 0.0069] | 11.3x | 0.0% | matched |
| basin | transport geography | straight-line, exponential, 24 km | 0.0010 | [0.0003, 0.0033] | 20.2x | 0.0% | matched |
| basin | transport geography | straight-line, exponential, 48 km | 0.0008 | [0.0002, 0.0028] | 25.8x | 0.0% | matched |
| basin | transport geography | straight-line, exponential, 96 km | 0.0009 | [0.0002, 0.0030] | 23.1x | 0.0% | matched |
| basin | transport geography | river, gaussian, 6 km | 0.0098 | [0.0026, 0.0354] | 2.1x | 16.8% | broken |
| basin | transport geography | river, gaussian, 12 km | 0.0073 | [0.0017, 0.0262] | 2.8x | 8.4% | broken |
| basin | transport geography | river, gaussian, 24 km | 0.0061 | [0.0015, 0.0219] | 3.4x | 4.8% | broken |
| basin | transport geography | river, gaussian, 48 km | 0.0030 | [0.0006, 0.0124] | 6.9x | 0.8% | matched |
| basin | transport geography | river, gaussian, 96 km | 0.0012 | [0.0003, 0.0048] | 17.3x | 0.0% | matched |
| basin | transport geography | river, power, 6 km | 0.0020 | [0.0003, 0.0105] | 10.2x | 0.4% | matched |
| basin | transport geography | river, power, 12 km | 0.0016 | [0.0003, 0.0064] | 12.7x | 0.0% | matched |
| basin | transport geography | river, power, 24 km | 0.0012 | [0.0003, 0.0042] | 16.6x | 0.0% | matched |
| basin | transport geography | river, power, 48 km | 0.0010 | [0.0003, 0.0034] | 19.9x | 0.0% | matched |
| basin | transport geography | river, power, 96 km | 0.0009 | [0.0003, 0.0032] | 22.1x | 0.0% | matched |
| basin | accumulation spans | window 8 +/- 0 | 0.0029 | [0.0009, 0.0135] | 7.0x | 1.2% | matched |
| basin | accumulation spans | window 8 +/- 2 | 0.0029 | [0.0008, 0.0133] | 7.2x | 1.2% | matched |
| basin | accumulation spans | window 8 +/- 4 | 0.0029 | [0.0007, 0.0143] | 7.1x | 1.2% | matched |
| basin | accumulation spans | window 8 +/- 6 | 0.0029 | [0.0007, 0.0131] | 7.2x | 1.2% | matched |
| basin | innovation boundary | boundary strength 0.0 | 0.0029 | [0.0009, 0.0135] | 7.0x | 1.2% | matched |
| basin | innovation boundary | boundary strength 0.1 | 0.0038 | [0.0008, 0.0535] | 5.4x | 8.0% | matched |
| basin | innovation boundary | boundary strength 0.15 | 0.0043 | [0.0009, 0.0701] | 4.9x | 14.0% | matched |
| basin | innovation boundary | boundary strength 0.2 | 0.0072 | [0.0008, 0.0818] | 2.9x | 23.6% | matched |
| basin | innovation boundary | boundary strength 0.25 | 0.0110 | [0.0013, 0.0963] | 1.9x | 32.8% | broken |
| basin | innovation boundary | boundary strength 0.3 | 0.0133 | [0.0013, 0.1019] | 1.6x | 38.0% | broken |
| basin | innovation boundary | boundary strength 0.4 | 0.0268 | [0.0019, 0.1014] | 0.8x | 56.0% | broken |
| basin | innovation boundary | boundary strength 0.6 | 0.0470 | [0.0048, 0.1051] | 0.4x | 81.2% | broken |
| basin | innovation boundary | boundary strength 1.0 | 0.0579 | [0.0197, 0.0938] | 0.4x | 96.8% | broken |
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

- **transport geography**: best shortfall 2.1x at river, gaussian, 6 km; the best cell that keeps diversity matched is river, exponential, 12 km at 4.0x.
- **accumulation spans**: best shortfall 7.0x at window 8 +/- 0; the best cell that keeps diversity matched is window 8 +/- 0 at 7.0x.
- **innovation boundary**: best shortfall 0.4x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.2 at 2.9x.

**cmv**

- **transport geography**: best shortfall 2.6x at river, gaussian, 6 km; the best cell that keeps diversity matched is river, gaussian, 6 km at 2.6x.
- **accumulation spans**: best shortfall 2.9x at window 8 +/- 0; the best cell that keeps diversity matched is window 8 +/- 0 at 2.9x.
- **innovation boundary**: best shortfall 0.1x at boundary strength 1.0; the best cell that keeps diversity matched is boundary strength 0.1 at 1.8x.
