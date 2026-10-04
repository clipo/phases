# Does connectivity-scaled participation close the shortfall?

Basin phase set, 28 assemblages, k = 5, 200 realizations each, calibrated pooled model (innovation 0.001, mixing 0.02, 2000 learners).
Per-node mixing is proportional to pre-normalisation connectivity, rescaled
so the mean stays at 0.02. Range 0.0035 to 0.0298.

| model | observed F_ST | drift median | shortfall |
|---|---|---|---|
| uniform participation | 0.0307 | 0.0044 | **7.0x** |
| scaled by connectivity | 0.0307 | 0.0047 | **6.5x** |

## Per cluster

| cluster | n | raw connectivity | its mixing rate | excess, uniform | excess, scaled |
|---|---|---|---|---|---|
| 1.0 | 5.0 | 1.96 | 0.0169 | 41.2x | **34.9x** |
| 2.0 | 5.0 | 2.52 | 0.0217 | 8.1x | **8.3x** |
| 0.0 | 10.0 | 2.46 | 0.0212 | 4.7x | **4.7x** |
| 4.0 | 4.0 | 2.46 | 0.0212 | 4.2x | **3.7x** |
| 3.0 | 4.0 | 2.05 | 0.0176 | 3.5x | **3.2x** |

## Did the diversity match survive?

The relaxation counts only if the assemblages' own diversity stays matched;
a model that reaches the observed differentiation by abandoning the
calibration has fitted one quantity at the other's expense.

| summary | observed | uniform | scaled |
|---|---|---|---|
| hs | 0.5633 | 0.5790 | 0.5780 |
| rich | 6.2500 | 6.0000 | 5.9821 |
| ht | 0.5916 | 0.5900 | 0.5877 |
