# Does connectivity-scaled participation close the shortfall?

Basin phase set, 43 assemblages, k = 5, 150 realizations each, calibrated pooled model (innovation 0.0005, mixing 0.02, 2000 learners).
Per-node mixing is proportional to pre-normalisation connectivity, rescaled
so the mean stays at 0.02. Range 0.0025 to 0.0405.

| model | observed F_ST | drift median | shortfall |
|---|---|---|---|
| uniform participation | 0.0280 | 0.0039 | **7.2x** |
| scaled by connectivity | 0.0280 | 0.0041 | **6.9x** |

## Per cluster

| cluster | n | raw connectivity | its mixing rate | excess, uniform | excess, scaled |
|---|---|---|---|---|---|
| 1.0 | 4.0 | 1.78 | 0.0124 | 77.7x | **57.4x** |
| 2.0 | 3.0 | 0.81 | 0.0057 | 15.8x | **6.9x** |
| 3.0 | 11.0 | 2.97 | 0.0207 | 12.2x | **15.7x** |
| 0.0 | 10.0 | 2.18 | 0.0152 | 4.1x | **3.7x** |
| 4.0 | 15.0 | 3.96 | 0.0276 | 0.5x | **0.5x** |

## Did the diversity match survive?

The relaxation counts only if the assemblages' own diversity stays matched;
a model that reaches the observed differentiation by abandoning the
calibration has fitted one quantity at the other's expense.

| summary | observed | uniform | scaled |
|---|---|---|---|
| hs | 0.5082 | 0.5113 | 0.5077 |
| rich | 6.0930 | 5.5000 | 5.5000 |
| ht | 0.5241 | 0.5164 | 0.5226 |
