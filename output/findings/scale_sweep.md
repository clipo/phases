# Is the drift verdict scale-invariant?

Basin phase set, 28 assemblages, 300 drift realizations,
calibrated pooled model (innovation 0.001, mixing 0.02, 2000 learners). Each realization is scored at
every k, so differences across scales are the partition and not a different draw.

| k | silhouette | observed F_ST | drift median | drift 95% | observed above? | shortfall | share of runs reaching observed |
|---|---|---|---|---|---|---|---|
| 2 | 0.533 | **0.0062** | 0.0021 | [0.0002, 0.0108] | NO | 3.0x | 11.7% |
| 3 | 0.524 | **0.0156** | 0.0036 | [0.0006, 0.0132] | yes | 4.3x | 1.3% |
| 4 | 0.488 | **0.0186** | 0.0040 | [0.0011, 0.0135] | yes | 4.7x | 0.7% |
| 5 | 0.503 | **0.0307** | 0.0046 | [0.0015, 0.0148] | yes | 6.7x | 0.0% |
| 6 | 0.425 | **0.0310** | 0.0050 | [0.0017, 0.0157] | yes | 6.1x | 0.0% |
| 7 | 0.450 | **0.0328** | 0.0055 | [0.0019, 0.0160] | yes | 5.9x | 0.0% |
| 8 | 0.442 | **0.0384** | 0.0057 | [0.0022, 0.0165] | yes | 6.7x | 0.0% |
| 9 | 0.461 | **0.0461** | 0.0062 | [0.0024, 0.0164] | yes | 7.4x | 0.0% |
| 10 | 0.490 | **0.0508** | 0.0066 | [0.0027, 0.0172] | yes | 7.7x | 0.0% |
| 11 | 0.468 | **0.0464** | 0.0072 | [0.0030, 0.0174] | yes | 6.5x | 0.0% |
| 12 | 0.438 | **0.0521** | 0.0076 | [0.0032, 0.0184] | yes | 6.8x | 0.0% |

**Verdict: the comparison does NOT hold at every scale; see the table.**

The level is a function of the scale: the observed F_ST rises with k because
a finer cut of a continuous distribution puts more variance between groups.
The comparison is not, and that is the quantity the argument uses.
