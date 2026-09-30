# Does connectivity-scaled participation close the shortfall?

Basin phase set, 28 assemblages, k = 5, 200 realizations each, calibrated pooled model (innovation 0.001, mixing 0.02, 2000 learners).
Per-node mixing is proportional to pre-normalisation connectivity, rescaled
so the mean stays at 0.02. Range 0.0031 to 0.0295.

| model | observed F_ST | drift median | shortfall |
|---|---|---|---|
| uniform participation | 0.0206 | 0.0046 | **4.5x** |
| scaled by connectivity | 0.0206 | 0.0049 | **4.2x** |

## Per cluster

| cluster | n | raw connectivity | its mixing rate | excess, uniform | excess, scaled |
|---|---|---|---|---|---|
| 2.0 | 1.0 | 0.35 | 0.0031 | 15.7x | **4.4x** |
| 1.0 | 5.0 | 2.52 | 0.0228 | 8.5x | **9.8x** |
| 4.0 | 8.0 | 2.22 | 0.0201 | 7.1x | **5.5x** |
| 0.0 | 9.0 | 2.35 | 0.0212 | 4.0x | **4.0x** |
| 3.0 | 5.0 | 2.02 | 0.0183 | 3.8x | **4.8x** |

## Did the diversity match survive?

The relaxation counts only if the assemblages' own diversity stays matched;
a model that reaches the observed differentiation by abandoning the
calibration has fitted one quantity at the other's expense.

| summary | observed | uniform | scaled |
|---|---|---|---|
| hs | 0.5633 | 0.5788 | 0.5796 |
| rich | 6.2500 | 6.0000 | 5.9821 |
| ht | 0.5916 | 0.5879 | 0.5894 |
