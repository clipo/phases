# Does connectivity-scaled participation close the shortfall?

Basin phase set, 28 assemblages, k = 5, 200 realizations each, calibrated pooled model (innovation 0.002, mixing 0.005, 2000 learners).
Per-node mixing is proportional to pre-normalisation connectivity, rescaled
so the mean stays at 0.005. Range 0.0007 to 0.0078.

| model | observed F_ST | drift median | shortfall |
|---|---|---|---|
| uniform participation | 0.0307 | 0.0067 | **4.6x** |
| scaled by connectivity | 0.0307 | 0.0072 | **4.3x** |

## Per cluster

| cluster | n | raw connectivity | its mixing rate | excess, uniform | excess, scaled |
|---|---|---|---|---|---|
| 1.0 | 5.0 | 2.17 | 0.0048 | 33.0x | **32.7x** |
| 2.0 | 5.0 | 2.67 | 0.0059 | 8.1x | **6.6x** |
| 0.0 | 10.0 | 2.14 | 0.0047 | 3.9x | **4.1x** |
| 4.0 | 4.0 | 2.74 | 0.0060 | 2.3x | **1.7x** |
| 3.0 | 4.0 | 1.84 | 0.0040 | 1.9x | **1.9x** |

## Did the diversity match survive?

The relaxation counts only if the assemblages' own diversity stays matched;
a model that reaches the observed differentiation by abandoning the
calibration has fitted one quantity at the other's expense.

| summary | observed | uniform | scaled |
|---|---|---|---|
| hs | 0.5633 | 0.5697 | 0.5696 |
| rich | 6.2500 | 5.7500 | 5.7143 |
| ht | 0.5916 | 0.5872 | 0.5891 |
