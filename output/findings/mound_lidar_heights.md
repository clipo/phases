# Mound height from lidar

One procedure, applied the same way at each site; see the header of
`analyses/68_mound_lidar_heights.py` for its definition. The spread across
the grid is what the choice of edge and baseline costs, not a confidence
interval.

## 11-N-1 (11-N-1)

- Recorded in the 2006 compilation: **23.0 ft**
- Lidar, at the reference operating point (edge 0.5 m, baseline the median of the 10-40 m annulus): **9.0 ft (2.73 m)**, footprint 16881 m2
- Across the 18-cell sensitivity grid: 5.6 to 39.5 ft

| edge (m) | annulus | pctl | footprint (m2) | height (m) | height (ft) |
|---|---|---|---|---|---|
| 0.3 | 10-40 m | 25 | 106395 | 9.38 | 30.8 |
| 0.3 | 10-40 m | 50 | 106395 | 8.02 | 26.3 |
| 0.3 | 10-40 m | 75 | 106395 | 7.60 | 24.9 |
| 0.3 | 20-60 m | 25 | 106395 | 12.03 | 39.5 |
| 0.3 | 20-60 m | 50 | 106395 | 8.01 | 26.3 |
| 0.3 | 20-60 m | 75 | 106395 | 7.67 | 25.2 |
| 0.5 | 10-40 m | 25 | 16881 | 4.21 | 13.8 |
| 0.5 | 10-40 m | 50 | 16881 | 2.73 | 9.0 |
| 0.5 | 10-40 m | 75 | 16881 | 1.70 | 5.6 |
| 0.5 | 20-60 m | 25 | 16881 | 3.81 | 12.5 |
| 0.5 | 20-60 m | 50 | 16881 | 2.68 | 8.8 |
| 0.5 | 20-60 m | 75 | 16881 | 1.82 | 6.0 |
| 1.0 | 10-40 m | 25 | 7178 | 6.61 | 21.7 |
| 1.0 | 10-40 m | 50 | 7178 | 3.64 | 11.9 |
| 1.0 | 10-40 m | 75 | 7178 | 1.71 | 5.6 |
| 1.0 | 20-60 m | 25 | 7178 | 5.92 | 19.4 |
| 1.0 | 20-60 m | 50 | 7178 | 2.86 | 9.4 |
| 1.0 | 20-60 m | 75 | 7178 | 1.71 | 5.6 |

## 13-N-3 (13-N-3)

- Recorded in the 2006 compilation: **23.0 ft**
- Lidar, at the reference operating point (edge 0.5 m, baseline the median of the 10-40 m annulus): **21.3 ft (6.48 m)**, footprint 4320 m2
- Across the 18-cell sensitivity grid: 20.7 to 22.3 ft

| edge (m) | annulus | pctl | footprint (m2) | height (m) | height (ft) |
|---|---|---|---|---|---|
| 0.3 | 10-40 m | 25 | 4527 | 6.68 | 21.9 |
| 0.3 | 10-40 m | 50 | 4527 | 6.49 | 21.3 |
| 0.3 | 10-40 m | 75 | 4527 | 6.35 | 20.8 |
| 0.3 | 20-60 m | 25 | 4527 | 6.78 | 22.3 |
| 0.3 | 20-60 m | 50 | 4527 | 6.57 | 21.5 |
| 0.3 | 20-60 m | 75 | 4527 | 6.40 | 21.0 |
| 0.5 | 10-40 m | 25 | 4320 | 6.67 | 21.9 |
| 0.5 | 10-40 m | 50 | 4320 | 6.48 | 21.3 |
| 0.5 | 10-40 m | 75 | 4320 | 6.34 | 20.8 |
| 0.5 | 20-60 m | 25 | 4320 | 6.77 | 22.2 |
| 0.5 | 20-60 m | 50 | 4320 | 6.56 | 21.5 |
| 0.5 | 20-60 m | 75 | 4320 | 6.40 | 21.0 |
| 1.0 | 10-40 m | 25 | 3860 | 6.65 | 21.8 |
| 1.0 | 10-40 m | 50 | 3860 | 6.47 | 21.2 |
| 1.0 | 10-40 m | 75 | 3860 | 6.32 | 20.7 |
| 1.0 | 20-60 m | 25 | 3860 | 6.76 | 22.2 |
| 1.0 | 20-60 m | 50 | 3860 | 6.55 | 21.5 |
| 1.0 | 20-60 m | 75 | 3860 | 6.40 | 21.0 |
