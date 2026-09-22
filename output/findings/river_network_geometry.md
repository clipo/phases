# River-network geometry for the basin assemblages

Produced by `analyses/61_river_network_geometry.py`. 28 assemblages, 378 pairs.

The along-waterway metric used by the copying kernel in `33_time_aware_emergence.py` and `34_emergence_robustness.py`.

| quantity | value |
|---|---|
| assemblages connecting on the network | 28 of 28 |
| river-graph largest component | 1939 nodes |
| unreachable pairs | 0 |
| mean access distance to a mapped channel | 2.75 km |
| max access distance | 8.1 km |
| detour ratio, mean of per-pair ratios | 3.68 |
| detour ratio, median of per-pair ratios | 2.91 |
| detour ratio, ratio of mean distances | 3.11 |

## Reading

The three detour summaries differ because the per-pair ratio is right-skewed: a pair that is close in a straight line but must route around a meander carries a very large ratio and pulls the mean above the median. Quoting an unnamed 'average' is therefore ambiguous at the 0.7 level here, which is why the supplement names the summary it uses.

None of these figures enters a result. The simulations consume the distance matrix itself, never a summary of it.


## Which distance does the pottery follow?

Chi-square distance between decorated-class profiles, 378 pairs.

| | rank correlation with ceramic distance |
|---|---|
| straight-line distance | 0.358 |
| river-network distance | 0.238 |
| river, straight-line held fixed | -0.026 |
| straight-line, river held fixed | 0.277 |

The two metrics correlate at 0.710. If the pottery followed the waterways, the river
metric would carry information the straight-line one lacks. Read the third row for that. A null
there does not show that waterways were unimportant: the network is MODERN hydrography, and the
St. Francis, Tyronza and Mississippi have all moved since these sites were occupied, so the metric
may be measuring the wrong rivers.

### Within phases and between them

| pairs | n | straight-line | river | river, straight-line held fixed | straight-line, river held fixed |
|---|---|---|---|---|---|
| within a phase | 119 | 0.400 | 0.277 | 0.094 | 0.313 |
| between phases | 259 | 0.318 | 0.102 | -0.109 | 0.320 |

### The phases as units

Pooled class profile per phase; distances are means over member pairs.

| phase pair | ceramic distance | straight-line km | river km | detour |
|---|---|---|---|---|
| Parkin - Walls | 0.361 | 54.2 | 192.8 | 3.56 |
| Kent - Walls | 0.424 | 30.9 | 93.6 | 3.03 |
| Kent - Parkin | 0.459 | 52.4 | 127.9 | 2.44 |

Across the 3 phase pairs ceramic distance rank-correlates -0.500 with straight-line distance and -0.500 with river distance. Ten pairs is indicative, no more.

The largest detours join phases on the St. Francis to phases on the Mississippi, which the network can
connect only through their confluence far to the south. Those distances describe the routing of the modern
network, and they are what the drift model's copying kernel uses.

## The grain of the record

Nearest-neighbour spacing of the 28 assemblages: median 4.3 km, quartiles 3.2 to 8.0 km. No spatial statistic here
can resolve structure much finer than that.