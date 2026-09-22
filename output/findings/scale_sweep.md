# Is the drift verdict scale-invariant?

Basin phase set, 28 assemblages, 300 drift realizations,
calibrated pooled model (innovation 0.024, mixing 0.01, 120 learners). Each realization is scored at
every k, so differences across scales are the partition and not a different draw.

| k | silhouette | observed F_ST | drift median | drift 95% | observed above? | shortfall | share of runs reaching observed |
|---|---|---|---|---|---|---|---|
| 2 | 0.529 | **0.0062** | 0.0017 | [0.0002, 0.0098] | NO | 3.6x | 9.7% |
| 3 | 0.501 | **0.0156** | 0.0043 | [0.0009, 0.0134] | yes | 3.6x | 1.7% |
| 4 | 0.469 | **0.0186** | 0.0054 | [0.0016, 0.0157] | yes | 3.5x | 0.3% |
| 5 | 0.469 | **0.0307** | 0.0064 | [0.0021, 0.0185] | yes | 4.8x | 0.0% |
| 6 | 0.371 | **0.0307** | 0.0070 | [0.0024, 0.0187] | yes | 4.4x | 0.0% |
| 7 | 0.352 | **0.0308** | 0.0091 | [0.0031, 0.0215] | yes | 3.4x | 0.0% |
| 8 | 0.358 | **0.0330** | 0.0112 | [0.0042, 0.0241] | yes | 3.0x | 0.0% |
| 9 | 0.433 | **0.0470** | 0.0115 | [0.0052, 0.0260] | yes | 4.1x | 0.0% |
| 10 | 0.423 | **0.0435** | 0.0135 | [0.0062, 0.0279] | yes | 3.2x | 0.0% |
| 11 | 0.450 | **0.0446** | 0.0145 | [0.0072, 0.0302] | yes | 3.1x | 0.0% |
| 12 | 0.386 | **0.0480** | 0.0140 | [0.0067, 0.0283] | yes | 3.4x | 0.0% |

**Verdict: the comparison does NOT hold at every scale; see the table.**

The level is a function of the scale: the observed F_ST rises with k because
a finer cut of a continuous distribution puts more variance between groups.
The comparison is not, and that is the quantity the argument uses.
