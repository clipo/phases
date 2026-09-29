# The relaxations scored at every scale

Produced by `analyses/92_relaxations_by_scale.py`: 300 runs per setting at the calibrated combination (N 2000, innovation 0.001, mixing 0.02), with analysis 64's and 65's seeds; local innovation areas are the two spatial clusters, as in 65. Each cell gives the median between-group F_ST, its 95 percent range, and the share of runs at or above the observed value. Diversity is matched when the medians of all three within-assemblage summaries lie within 10 percent of the observed values. No setting is recalibrated.

| account | setting | 2 clusters (obs 0.0062) | 3 clusters (obs 0.0156) | 4 clusters (obs 0.0186) | phases (obs 0.0110) | diversity |
|---|---|---|---|---|---|---|
| baseline | calibrated copying model | 0.0022 [0.0003, 0.0119]; 13% | 0.0040 [0.0009, 0.0158]; 3% | 0.0042 [0.0014, 0.0166]; 1% | 0.0037 [0.0008, 0.0156]; 7% | matched |
| unequal populations | CV 0.5 | 0.0025 [0.0003, 0.0142]; 15% | 0.0042 [0.0009, 0.0180]; 4% | 0.0048 [0.0013, 0.0188]; 3% | 0.0040 [0.0009, 0.0178]; 6% | matched |
| unequal populations | CV 1.0 | 0.0033 [0.0004, 0.0219]; 29% | 0.0064 [0.0013, 0.0248]; 13% | 0.0069 [0.0015, 0.0258]; 9% | 0.0065 [0.0013, 0.0242]; 24% | broken |
| unequal populations | CV 1.5 | 0.0046 [0.0005, 0.0392]; 42% | 0.0089 [0.0016, 0.0483]; 26% | 0.0100 [0.0029, 0.0488]; 21% | 0.0084 [0.0015, 0.0437]; 40% | broken |
| transport | straight-line, exponential, 6 km | 0.0034 [0.0004, 0.0213]; 32% | 0.0050 [0.0010, 0.0241]; 11% | 0.0055 [0.0014, 0.0243]; 7% | 0.0050 [0.0010, 0.0236]; 19% | matched |
| transport | river, exponential, 6 km | 0.0071 [0.0007, 0.0328]; 54% | 0.0112 [0.0021, 0.0412]; 31% | 0.0133 [0.0029, 0.0455]; 31% | 0.0110 [0.0019, 0.0381]; 50% | broken |
| transport | river, exponential, 12 km | 0.0048 [0.0005, 0.0255]; 37% | 0.0074 [0.0015, 0.0306]; 16% | 0.0086 [0.0022, 0.0322]; 12% | 0.0072 [0.0014, 0.0299]; 27% | broken |
| transport | river, gaussian, 24 km | 0.0046 [0.0005, 0.0300]; 40% | 0.0086 [0.0014, 0.0351]; 18% | 0.0098 [0.0022, 0.0361]; 15% | 0.0084 [0.0014, 0.0337]; 33% | broken |
| transport | river, gaussian, 48 km | 0.0023 [0.0002, 0.0131]; 17% | 0.0043 [0.0009, 0.0187]; 3% | 0.0046 [0.0011, 0.0190]; 3% | 0.0042 [0.0008, 0.0176]; 8% | matched |
| local innovation | strength 0.2 | 0.0047 [0.0004, 0.0530]; 44% | 0.0063 [0.0011, 0.0544]; 31% | 0.0070 [0.0015, 0.0555]; 28% | 0.0061 [0.0010, 0.0545]; 38% | matched |
| local innovation | strength 0.25 | 0.0086 [0.0004, 0.0539]; 54% | 0.0107 [0.0011, 0.0571]; 42% | 0.0109 [0.0016, 0.0578]; 40% | 0.0103 [0.0010, 0.0565]; 49% | matched |
| local innovation | strength 0.3 | 0.0111 [0.0005, 0.0532]; 61% | 0.0125 [0.0011, 0.0551]; 44% | 0.0132 [0.0018, 0.0562]; 42% | 0.0123 [0.0012, 0.0545]; 53% | matched |

## Reading

Read each row on three things at once: whether the median reaches the observation, whether the 95 percent range covers it, and whether diversity stays matched. The baseline row is the reference; whether its range covers a finer-scale value depends on the seeds (compare Table 1), so read the share of runs reaching it rather than the range edge. A relaxation improves on the baseline where it raises that share with diversity matched. Settings that break the diversity match here were not recalibrated.
