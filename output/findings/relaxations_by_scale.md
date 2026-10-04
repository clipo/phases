# The relaxations scored at every scale

Produced by `analyses/92_relaxations_by_scale.py`: 300 runs per setting at the calibrated combination (N 2000, innovation 0.001, mixing 0.02), with analysis 64's and 65's seeds; local innovation areas are the two spatial clusters, as in 65, unless the row names neighborhoods (the five- and eight-cluster k-means divisions of the site map) or single sites. Each cell gives the median between-group F_ST, its 95 percent range, and the share of runs at or above the observed value. Diversity is matched when the medians of all three within-assemblage summaries lie within 10 percent of the observed values. No setting is recalibrated.

| account | setting | 2 clusters (obs 0.0062) | 3 clusters (obs 0.0156) | 4 clusters (obs 0.0186) | phases (obs 0.0123) | diversity |
|---|---|---|---|---|---|---|
| baseline | calibrated copying model | 0.0017 [0.0002, 0.0113]; 12% | 0.0035 [0.0006, 0.0134]; 2% | 0.0040 [0.0010, 0.0138]; 1% | 0.0024 [0.0005, 0.0098]; 1% | matched |
| unequal populations | CV 0.5 | 0.0025 [0.0003, 0.0124]; 17% | 0.0044 [0.0009, 0.0146]; 2% | 0.0049 [0.0013, 0.0153]; 1% | 0.0028 [0.0005, 0.0122]; 2% | matched |
| unequal populations | CV 1.0 | 0.0040 [0.0003, 0.0278]; 37% | 0.0073 [0.0011, 0.0314]; 16% | 0.0078 [0.0016, 0.0314]; 14% | 0.0047 [0.0008, 0.0250]; 15% | broken |
| unequal populations | CV 1.5 | 0.0052 [0.0003, 0.0379]; 47% | 0.0099 [0.0014, 0.0453]; 29% | 0.0108 [0.0019, 0.0454]; 24% | 0.0067 [0.0010, 0.0325]; 24% | broken |
| transport | straight-line, exponential, 6 km | 0.0035 [0.0004, 0.0194]; 26% | 0.0049 [0.0011, 0.0208]; 7% | 0.0052 [0.0013, 0.0213]; 5% | 0.0034 [0.0007, 0.0191]; 7% | matched |
| transport | river, exponential, 6 km | 0.0056 [0.0006, 0.0325]; 45% | 0.0087 [0.0019, 0.0397]; 26% | 0.0114 [0.0033, 0.0427]; 24% | 0.0077 [0.0014, 0.0337]; 29% | broken |
| transport | river, exponential, 12 km | 0.0047 [0.0004, 0.0247]; 39% | 0.0080 [0.0011, 0.0278]; 16% | 0.0086 [0.0019, 0.0297]; 12% | 0.0055 [0.0008, 0.0249]; 18% | broken |
| transport | river, gaussian, 24 km | 0.0043 [0.0005, 0.0315]; 38% | 0.0075 [0.0013, 0.0357]; 20% | 0.0089 [0.0022, 0.0372]; 17% | 0.0053 [0.0012, 0.0309]; 19% | broken |
| transport | river, gaussian, 48 km | 0.0026 [0.0002, 0.0130]; 18% | 0.0042 [0.0007, 0.0163]; 4% | 0.0047 [0.0010, 0.0169]; 1% | 0.0029 [0.0005, 0.0125]; 3% | matched |
| local innovation | strength 0.2 | 0.0054 [0.0003, 0.0566]; 47% | 0.0071 [0.0011, 0.0588]; 32% | 0.0075 [0.0013, 0.0597]; 28% | 0.0054 [0.0006, 0.0529]; 33% | matched |
| local innovation | strength 0.25 | 0.0094 [0.0004, 0.0608]; 57% | 0.0110 [0.0011, 0.0634]; 44% | 0.0113 [0.0014, 0.0655]; 41% | 0.0087 [0.0006, 0.0569]; 45% | matched |
| local innovation | strength 0.3 | 0.0139 [0.0006, 0.0615]; 63% | 0.0157 [0.0013, 0.0649]; 50% | 0.0165 [0.0016, 0.0667]; 47% | 0.0137 [0.0006, 0.0589]; 53% | matched |
| local innovation, five neighborhoods | strength 0.2 | 0.0050 [0.0004, 0.0448]; 45% | 0.0078 [0.0011, 0.0496]; 22% | 0.0081 [0.0016, 0.0512]; 20% | 0.0058 [0.0007, 0.0435]; 25% | broken |
| local innovation, five neighborhoods | strength 0.25 | 0.0068 [0.0006, 0.0503]; 52% | 0.0099 [0.0014, 0.0535]; 30% | 0.0107 [0.0017, 0.0545]; 28% | 0.0073 [0.0009, 0.0475]; 33% | broken |
| local innovation, eight neighborhoods | strength 0.2 | 0.0045 [0.0004, 0.0260]; 37% | 0.0066 [0.0014, 0.0288]; 15% | 0.0073 [0.0018, 0.0298]; 13% | 0.0050 [0.0010, 0.0237]; 19% | broken |
| local innovation, eight neighborhoods | strength 0.25 | 0.0067 [0.0005, 0.0276]; 52% | 0.0087 [0.0015, 0.0315]; 21% | 0.0096 [0.0021, 0.0318]; 15% | 0.0069 [0.0011, 0.0281]; 23% | broken |
| local innovation, each site its own | strength 0.1 | 0.0027 [0.0003, 0.0136]; 17% | 0.0043 [0.0007, 0.0149]; 2% | 0.0048 [0.0011, 0.0156]; 0% | 0.0031 [0.0005, 0.0128]; 3% | matched |
| local innovation, each site its own | strength 0.2 | 0.0031 [0.0005, 0.0131]; 21% | 0.0049 [0.0014, 0.0154]; 2% | 0.0054 [0.0019, 0.0154]; 1% | 0.0034 [0.0008, 0.0127]; 3% | broken |

## Reading

Read each row on three things at once: whether the median reaches the observation, whether the 95 percent range covers it, and whether diversity stays matched. The baseline row is the reference; whether its range covers a finer-scale value depends on the seeds (compare Table 1), so read the share of runs reaching it rather than the range edge. A relaxation improves on the baseline where it raises that share with diversity matched. Settings that break the diversity match here were not recalibrated.

## The dose behind each local-innovation row

An area's profile is reordered only when two or more of its classes are drawn, each with probability equal to the strength, and the reordering is not the identity. Measured over the same 300 runs: the mean share of areas whose profile differs from the regional one at all, and the mean share of a profile's frequency that sits in a different class (half the summed absolute difference from the regional profile, averaged over areas, unchanged ones included). Two rows are the same dose only when these agree; a row that reaches the finer-scale values no more often than the baseline at a small dose has not tested the grain it names.

| account | setting | areas | areas whose profile differs | frequency moved, mean |
|---|---|---|---|---|
| local innovation | strength 0.2 | 2 | 46% | 9.7% |
| local innovation | strength 0.25 | 2 | 59% | 14.1% |
| local innovation | strength 0.3 | 2 | 70% | 17.2% |
| local innovation, five neighborhoods | strength 0.2 | 5 | 46% | 9.4% |
| local innovation, five neighborhoods | strength 0.25 | 5 | 60% | 14.0% |
| local innovation, eight neighborhoods | strength 0.2 | 8 | 46% | 9.6% |
| local innovation, eight neighborhoods | strength 0.25 | 8 | 60% | 14.2% |
| local innovation, each site its own | strength 0.1 | 28 | 15% | 2.9% |
| local innovation, each site its own | strength 0.2 | 28 | 44% | 9.6% |
