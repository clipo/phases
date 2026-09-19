# River-network geometry for the basin assemblages

Produced by `analyses/61_river_network_geometry.py`. 29 assemblages, 406 pairs.

The along-waterway metric used by the copying kernel in `33_time_aware_emergence.py` and `34_emergence_robustness.py`.

| quantity | value |
|---|---|
| assemblages connecting on the network | 29 of 29 |
| river-graph largest component | 2016 nodes |
| unreachable pairs | 0 |
| mean access distance to a mapped channel | 3.12 km |
| max access distance | 10.6 km |
| detour ratio, mean of per-pair ratios | 3.10 |
| detour ratio, median of per-pair ratios | 2.39 |
| detour ratio, ratio of mean distances | 2.67 |

## Reading

The three detour summaries differ because the per-pair ratio is right-skewed: a pair that is close in a straight line but must route around a meander carries a very large ratio and pulls the mean above the median. Quoting an unnamed 'average' is therefore ambiguous at the 0.7 level here, which is why the supplement names the summary it uses.

None of these figures enters a result. The simulations consume the distance matrix itself, never a summary of it.
