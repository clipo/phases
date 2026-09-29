# Turnover along the reconstructed phase boundaries: bandwidth and recovery

Produced by `analyses/91_turnover_sensitivity.py`. Basin set, 28 assemblages; 160 x 160 grid (analysis 75 uses 220); 200 same-size alternative divisions around random centers, built as analysis 75 builds them (seed 75). The boundaries are this study's reconstruction of territories from Mainfort's (1996) assignments, not published lines. Percentile = share of alternatives whose median turnover along their boundaries is below the phase boundaries'; a boundary between communities would sit high.

## A. Bandwidth

| bandwidth (km) | phase boundaries, median turnover per km | alternatives, median | percentile of the phase boundaries |
|---:|---:|---:|---:|
| 10 | 0.01624 | 0.02130 | 3.0 |
| 12.5 | 0.01232 | 0.01641 | 0.0 |
| 15 | 0.01051 | 0.01240 | 0.0 |
| 20 | 0.00936 | 0.00986 | 23.5 |
| 25 | 0.00910 | 0.00922 | 35.5 |
| 30 | 0.00896 | 0.00901 | 43.0 |

## B. Recovery of a copying boundary at the phase lines

100 simulated records per cell from analysis 84's grid and seeds (N 2000, innovation 0.001, mixing 0.02, 24 km along the rivers), sampled at the real sherd counts, bandwidth 15 km. Entries are percentiles of the phase boundaries among the same alternatives.

| copying factor | local innovation | percentile, 5% / 25% / 50% / 75% / 95% over records | records at or below the observed percentile |
|---|---|---|---:|
| 1 | 0 | 1 / 8 / 20 / 48 / 88 | 0.05 |
| 0.1 | 0 | 3 / 12 / 33 / 59 / 85 | 0.02 |
| 0.03 | 0 | 3 / 26 / 54 / 77 / 94 | 0.02 |
| 0.03 | 0.2 | 10 / 54 / 73 / 85 / 94 | 0.01 |

Observed at this grid and bandwidth: median 0.01051 per km, percentile 0.0.

## Reading

If the percentile stays low across bandwidths, the result is not an artifact of 15 km. If the boundary cells in B do not lift the phase lines' percentile above the no-boundary cell, the comparison cannot see a copying boundary on this sampling geometry, and the turnover result is description rather than evidence against a boundary.
