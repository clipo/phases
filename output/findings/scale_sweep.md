# Is the drift verdict scale-invariant?

Basin phase set, 43 assemblages, 400 drift realizations,
calibrated pooled model (innovation 0.0005, mixing 0.02, 2000 learners). Each realization is scored at
every k, so differences across scales are the partition and not a different draw.

| k | silhouette | observed F_ST | drift median | drift 95% | observed above? | shortfall | share of runs reaching observed |
|---|---|---|---|---|---|---|---|
| 2 | 0.455 | **0.0076** | 0.0026 | [0.0002, 0.0136] | NO | 2.9x | 13.2% |
| 3 | 0.417 | **0.0163** | 0.0035 | [0.0004, 0.0135] | yes | 4.6x | 1.8% |
| 4 | 0.502 | **0.0207** | 0.0039 | [0.0007, 0.0143] | yes | 5.3x | 0.5% |
| 5 | 0.566 | **0.0280** | 0.0041 | [0.0008, 0.0151] | yes | 6.8x | 0.0% |
| 6 | 0.514 | **0.0420** | 0.0050 | [0.0012, 0.0168] | yes | 8.4x | 0.0% |

**Verdict: the comparison does NOT hold at every scale; see the table.**

The level is a function of the scale: the observed F_ST rises with k because
a finer cut of a continuous distribution puts more variance between groups.
The comparison is not, and that is the quantity the argument uses.
