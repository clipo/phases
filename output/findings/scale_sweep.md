# Is the drift verdict scale-invariant?

Basin phase set, 28 assemblages, 300 drift realizations,
calibrated pooled model (innovation 0.002, mixing 0.005, 2000 learners). Each realization is scored at
every k, so differences across scales are the partition and not a different draw.

| k | silhouette | observed F_ST | drift median | drift 95% | observed above? | shortfall | share of runs reaching observed |
|---|---|---|---|---|---|---|---|
| 2 | 0.529 | **0.0062** | 0.0025 | [0.0002, 0.0157] | NO | 2.5x | 17.0% |
| 3 | 0.502 | **0.0156** | 0.0047 | [0.0011, 0.0186] | NO | 3.3x | 5.0% |
| 4 | 0.469 | **0.0186** | 0.0057 | [0.0016, 0.0210] | NO | 3.3x | 3.7% |
| 5 | 0.470 | **0.0307** | 0.0067 | [0.0022, 0.0214] | yes | 4.6x | 1.0% |
| 6 | 0.383 | **0.0254** | 0.0068 | [0.0024, 0.0224] | yes | 3.7x | 1.3% |
| 7 | 0.353 | **0.0308** | 0.0082 | [0.0034, 0.0232] | yes | 3.7x | 1.3% |
| 8 | 0.415 | **0.0328** | 0.0103 | [0.0048, 0.0253] | yes | 3.2x | 1.3% |
| 9 | 0.433 | **0.0470** | 0.0108 | [0.0050, 0.0265] | yes | 4.3x | 0.3% |
| 10 | 0.421 | **0.0435** | 0.0124 | [0.0059, 0.0271] | yes | 3.5x | 0.7% |
| 11 | 0.449 | **0.0446** | 0.0129 | [0.0068, 0.0284] | yes | 3.5x | 0.7% |
| 12 | 0.383 | **0.0480** | 0.0125 | [0.0064, 0.0284] | yes | 3.9x | 0.0% |

**Verdict: the comparison does NOT hold at every scale; see the table.**

The level is a function of the scale: the observed F_ST rises with k because
a finer cut of a continuous distribution puts more variance between groups.
The comparison is not, and that is the quantity the argument uses.
