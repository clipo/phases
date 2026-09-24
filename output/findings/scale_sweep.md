# Is the drift verdict scale-invariant?

Basin phase set, 28 assemblages, 300 drift realizations,
calibrated pooled model (innovation 0.001, mixing 0.02, 2000 learners). Each realization is scored at
every k, so differences across scales are the partition and not a different draw.

| k | silhouette | observed F_ST | drift median | drift 95% | observed above? | shortfall | share of runs reaching observed |
|---|---|---|---|---|---|---|---|
| 2 | 0.527 | **0.0062** | 0.0019 | [0.0002, 0.0124] | NO | 3.3x | 13.3% |
| 3 | 0.506 | **0.0156** | 0.0035 | [0.0007, 0.0135] | yes | 4.5x | 1.3% |
| 4 | 0.466 | **0.0186** | 0.0041 | [0.0010, 0.0146] | yes | 4.6x | 1.0% |
| 5 | 0.433 | **0.0206** | 0.0043 | [0.0013, 0.0148] | yes | 4.8x | 0.3% |
| 6 | 0.453 | **0.0327** | 0.0050 | [0.0017, 0.0152] | yes | 6.5x | 0.0% |
| 7 | 0.475 | **0.0345** | 0.0053 | [0.0019, 0.0161] | yes | 6.5x | 0.0% |
| 8 | 0.418 | **0.0346** | 0.0058 | [0.0024, 0.0165] | yes | 6.0x | 0.0% |
| 9 | 0.432 | **0.0426** | 0.0065 | [0.0027, 0.0183] | yes | 6.5x | 0.0% |
| 10 | 0.444 | **0.0441** | 0.0067 | [0.0028, 0.0185] | yes | 6.6x | 0.0% |
| 11 | 0.414 | **0.0488** | 0.0065 | [0.0028, 0.0180] | yes | 7.5x | 0.0% |
| 12 | 0.435 | **0.0529** | 0.0078 | [0.0031, 0.0193] | yes | 6.8x | 0.0% |

**Verdict: the comparison does NOT hold at every scale; see the table.**

The level is a function of the scale: the observed F_ST rises with k because
a finer cut of a continuous distribution puts more variance between groups.
The comparison is not, and that is the quantity the argument uses.
