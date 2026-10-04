# Turnover along the reconstructed phase boundaries: bandwidth and recovery

Produced by `analyses/91_turnover_sensitivity.py`. Basin set, 28 assemblages; 160 x 160 grid (analysis 75 uses 220); 200 same-size alternative divisions around random centers, built as analysis 75 builds them (seed 75). The boundaries are this study's reconstruction of territories from the assemblages' Phillips (1970) phases, not his drawn lines. Percentile = share of alternatives whose median turnover along their boundaries is below the phase boundaries'; a boundary between communities would sit high.

## A. Bandwidth

| bandwidth (km) | phase boundaries, median turnover per km | alternatives, median | percentile of the phase boundaries |
|---:|---:|---:|---:|
| 10 | 0.01588 | 0.01900 | 6.0 |
| 12.5 | 0.01188 | 0.01492 | 0.0 |
| 15 | 0.01038 | 0.01168 | 7.5 |
| 20 | 0.00949 | 0.00977 | 40.0 |
| 25 | 0.00918 | 0.00926 | 41.0 |
| 30 | 0.00909 | 0.00911 | 46.0 |

## B. Recovery of a copying boundary at the phase lines

100 simulated records per cell from analysis 84's grid and seeds (N 2000, innovation 0.001, mixing 0.02, 24 km along the rivers), sampled at the real sherd counts, bandwidth 15 km. Entries are percentiles of the phase boundaries among the same alternatives.

| copying factor | local innovation | percentile, 5% / 25% / 50% / 75% / 95% over records | records at or below the observed percentile |
|---|---|---|---:|
| 1 | 0 | 4 / 12 / 21 / 41 / 76 | 0.12 |
| 0.1 | 0 | 6 / 22 / 41 / 61 / 86 | 0.07 |
| 0.03 | 0 | 12 / 26 / 48 / 69 / 84 | 0.03 |
| 0.03 | 0.2 | 14 / 45 / 76 / 84 / 88 | 0.02 |

Observed at this grid and bandwidth: median 0.01038 per km, percentile 7.5.

## Reading

If the percentile stays low across bandwidths, the result is not an artifact of 15 km. If the boundary cells in B do not lift the phase lines' percentile above the no-boundary cell, the comparison cannot see a copying boundary on this sampling geometry, and the turnover result is description rather than evidence against a boundary.
