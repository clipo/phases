# Does connectivity-scaled participation close the shortfall?

Basin phase set, 28 assemblages, k = 5, 200 realizations each, calibrated pooled model (innovation 0.024, mixing 0.01, 120 learners).
Per-node mixing is proportional to pre-normalisation connectivity, rescaled
so the mean stays at 0.01. Range 0.0015 to 0.0155.

| model | observed F_ST | drift median | shortfall |
|---|---|---|---|
| uniform participation | 0.0307 | 0.0064 | **4.8x** |
| scaled by connectivity | 0.0307 | 0.0062 | **5.0x** |

## Per cluster

| cluster | n | raw connectivity | its mixing rate | excess, uniform | excess, scaled |
|---|---|---|---|---|---|
| 1.0 | 5.0 | 2.17 | 0.0095 | 35.7x | **34.9x** |
| 2.0 | 5.0 | 2.66 | 0.0117 | 7.6x | **8.5x** |
| 0.0 | 10.0 | 2.14 | 0.0094 | 5.5x | **6.2x** |
| 4.0 | 4.0 | 2.72 | 0.0119 | 2.1x | **1.9x** |
| 3.0 | 4.0 | 1.85 | 0.0081 | 1.9x | **2.0x** |

## Did the diversity match survive?

The relaxation counts only if the assemblages' own diversity stays matched;
a model that reaches the observed differentiation by abandoning the
calibration has fitted one quantity at the other's expense.

| summary | observed | uniform | scaled |
|---|---|---|---|
| hs | 0.5633 | 0.5680 | 0.5677 |
| rich | 6.2500 | 5.6786 | 5.6786 |
| ht | 0.5916 | 0.5912 | 0.5912 |
