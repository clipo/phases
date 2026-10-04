# Could the phase-line comparison see a copying boundary if there were one?

Produced by `analyses/89_partition_power.py`: 100 simulated records per cell, rebuilt with analysis 84's seeds at its calibrated cell (N 2000, innovation 0.001, mixing 0.02, 24 km along the rivers); 50 same-size alternative divisions of each kind. Each entry is the median over simulated records of the share of alternatives the tested line beats, with the share of records in which it beats at least 90 percent of them in parentheses.

| copying factor at the phase lines | local innovation | line tested | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |
|---|---|---|---|---|---|---|
| 1.0 | 0.0 | phases | 0.38 (6%) | 0.38 (10%) | 0.65 (19%) | 0.63 (18%) |
| 1.0 | 0.0 | Parkin vs rest | 0.58 (7%) | 0.62 (19%) | 0.28 (2%) | 0.00 (22%) |
| 0.1 | 0.0 | phases | 0.67 (20%) | 0.69 (25%) | 0.86 (44%) | 0.84 (42%) |
| 0.1 | 0.0 | Parkin vs rest | 0.70 (14%) | 0.62 (37%) | 0.56 (21%) | 0.62 (42%) |
| 0.03 | 0.0 | phases | 0.83 (37%) | 0.79 (37%) | 0.92 (68%) | 0.94 (53%) |
| 0.03 | 0.0 | Parkin vs rest | 0.72 (24%) | 0.62 (43%) | 0.69 (29%) | 1.00 (52%) |
| 1.0 | 0.2 | phases | 0.56 (15%) | 0.46 (13%) | 0.64 (10%) | 0.60 (13%) |
| 1.0 | 0.2 | Parkin vs rest | 0.70 (8%) | 0.62 (33%) | 0.56 (15%) | 0.62 (41%) |
| 0.03 | 0.2 | phases | 0.95 (66%) | 1.00 (63%) | 0.96 (83%) | 1.00 (72%) |
| 0.03 | 0.2 | Parkin vs rest | 0.88 (48%) | 1.00 (72%) | 0.83 (37%) | 1.00 (66%) |

The observed record (plug-in, these alternatives):

| line | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |
|---|---|---|---|---|
| phases | 0.64 | 0.58 | 0.76 | 0.84 |
| Parkin vs rest | 0.58 | 1.00 | 0.42 | 0.38 |

## The distributions behind the medians

For each cell, the 5th, 25th, 50th, 75th and 95th percentiles of the share of alternatives the tested line beats, over simulated records, and the fraction of simulated records whose share is at or below the observed record's. A fraction of 0 means none of the simulated records scored that low, which bounds the frequency below about 1/reps, not at zero.

### 2000 learners, innovation 0.001, mixing 0.02 (calibrated combination)

| line | statistic | observed | copying factor | local innovation | 5% / 25% / 50% / 75% / 95% | records at or below observed |
|---|---|---:|---|---|---|---:|
| phases | F_ST vs random-center | 0.64 | 1.0 | 0.0 | 0.12 / 0.24 / 0.38 / 0.58 / 0.90 | 0.78 |
| phases | F_ST vs random-center | 0.64 | 0.1 | 0.0 | 0.20 / 0.43 / 0.67 / 0.86 / 0.96 | 0.49 |
| phases | F_ST vs random-center | 0.64 | 0.03 | 0.0 | 0.32 / 0.57 / 0.83 / 0.94 / 0.98 | 0.29 |
| phases | F_ST vs random-center | 0.64 | 1.0 | 0.2 | 0.20 / 0.38 / 0.56 / 0.77 / 0.96 | 0.64 |
| phases | F_ST vs random-center | 0.64 | 0.03 | 0.2 | 0.50 / 0.85 / 0.95 / 0.98 / 0.98 | 0.14 |
| phases | F_ST vs compact | 0.58 | 1.0 | 0.0 | 0.00 / 0.22 / 0.38 / 0.65 / 0.92 | 0.66 |
| phases | F_ST vs compact | 0.58 | 0.1 | 0.0 | 0.19 / 0.41 / 0.69 / 0.85 / 1.00 | 0.41 |
| phases | F_ST vs compact | 0.58 | 0.03 | 0.0 | 0.21 / 0.52 / 0.79 / 0.98 / 1.00 | 0.28 |
| phases | F_ST vs compact | 0.58 | 1.0 | 0.2 | 0.00 / 0.28 / 0.46 / 0.74 / 1.00 | 0.64 |
| phases | F_ST vs compact | 0.58 | 0.03 | 0.2 | 0.48 / 0.80 / 1.00 / 1.00 / 1.00 | 0.11 |
| phases | boundary excess vs random-center | 0.76 | 1.0 | 0.0 | 0.04 / 0.35 / 0.65 / 0.84 / 0.96 | 0.65 |
| phases | boundary excess vs random-center | 0.76 | 0.1 | 0.0 | 0.30 / 0.66 / 0.86 / 0.94 / 0.98 | 0.37 |
| phases | boundary excess vs random-center | 0.76 | 0.03 | 0.0 | 0.54 / 0.88 / 0.92 / 0.98 / 0.98 | 0.16 |
| phases | boundary excess vs random-center | 0.76 | 1.0 | 0.2 | 0.24 / 0.43 / 0.64 / 0.82 / 0.92 | 0.69 |
| phases | boundary excess vs random-center | 0.76 | 0.03 | 0.2 | 0.59 / 0.92 / 0.96 / 0.98 / 0.98 | 0.08 |
| phases | boundary excess vs compact | 0.84 | 1.0 | 0.0 | 0.00 / 0.32 / 0.63 / 0.84 / 1.00 | 0.82 |
| phases | boundary excess vs compact | 0.84 | 0.1 | 0.0 | 0.22 / 0.57 / 0.84 / 0.98 / 1.00 | 0.57 |
| phases | boundary excess vs compact | 0.84 | 0.03 | 0.0 | 0.41 / 0.78 / 0.94 / 1.00 / 1.00 | 0.45 |
| phases | boundary excess vs compact | 0.84 | 1.0 | 0.2 | 0.06 / 0.34 / 0.60 / 0.76 / 0.98 | 0.87 |
| phases | boundary excess vs compact | 0.84 | 0.03 | 0.2 | 0.54 / 0.84 / 1.00 / 1.00 / 1.00 | 0.28 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 1.0 | 0.0 | 0.12 / 0.36 / 0.58 / 0.72 / 0.90 | 0.52 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.1 | 0.0 | 0.12 / 0.50 / 0.70 / 0.86 / 0.90 | 0.36 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.03 | 0.0 | 0.30 / 0.60 / 0.72 / 0.86 / 0.90 | 0.23 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 1.0 | 0.2 | 0.42 / 0.64 / 0.70 / 0.86 / 0.90 | 0.20 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.03 | 0.2 | 0.46 / 0.68 / 0.88 / 0.90 / 0.90 | 0.16 |
| Parkin vs rest | F_ST vs compact | 1.00 | 1.0 | 0.0 | 0.00 / 0.38 / 0.62 / 0.62 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.1 | 0.0 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.03 | 0.0 | 0.00 / 0.62 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 1.0 | 0.2 | 0.00 / 0.62 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.03 | 0.2 | 0.38 / 0.62 / 1.00 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 1.0 | 0.0 | 0.00 / 0.07 / 0.28 / 0.58 / 0.86 | 0.67 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.1 | 0.0 | 0.00 / 0.28 / 0.56 / 0.83 / 0.90 | 0.37 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.03 | 0.0 | 0.00 / 0.36 / 0.69 / 0.90 / 0.90 | 0.30 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 1.0 | 0.2 | 0.00 / 0.26 / 0.56 / 0.79 / 0.90 | 0.40 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.03 | 0.2 | 0.04 / 0.39 / 0.83 / 0.90 / 0.90 | 0.27 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 1.0 | 0.0 | 0.00 / 0.00 / 0.00 / 0.62 / 1.00 | 0.68 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.1 | 0.0 | 0.00 / 0.00 / 0.62 / 1.00 / 1.00 | 0.44 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.03 | 0.0 | 0.00 / 0.38 / 1.00 / 1.00 / 1.00 | 0.39 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 1.0 | 0.2 | 0.00 / 0.00 / 0.62 / 1.00 / 1.00 | 0.42 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.03 | 0.2 | 0.00 / 0.38 / 1.00 / 1.00 / 1.00 | 0.27 |

### 10000 learners, innovation 0.0005, mixing 0.002 (sensitivity combination)

| line | statistic | observed | copying factor | local innovation | 5% / 25% / 50% / 75% / 95% | records at or below observed |
|---|---|---:|---|---|---|---:|
| phases | F_ST vs random-center | 0.64 | 1.0 | 0.0 | 0.10 / 0.26 / 0.46 / 0.65 / 0.90 | 0.75 |
| phases | F_ST vs random-center | 0.64 | 0.1 | 0.0 | 0.12 / 0.38 / 0.68 / 0.86 / 0.98 | 0.48 |
| phases | F_ST vs random-center | 0.64 | 0.03 | 0.0 | 0.16 / 0.43 / 0.75 / 0.92 / 0.98 | 0.43 |
| phases | F_ST vs random-center | 0.64 | 1.0 | 0.2 | 0.14 / 0.42 / 0.67 / 0.90 / 0.98 | 0.46 |
| phases | F_ST vs random-center | 0.64 | 0.03 | 0.2 | 0.22 / 0.61 / 0.92 / 0.98 / 0.98 | 0.27 |
| phases | F_ST vs compact | 0.58 | 1.0 | 0.0 | 0.02 / 0.22 / 0.48 / 0.72 / 0.94 | 0.64 |
| phases | F_ST vs compact | 0.58 | 0.1 | 0.0 | 0.02 / 0.34 / 0.70 / 0.92 / 1.00 | 0.41 |
| phases | F_ST vs compact | 0.58 | 0.03 | 0.0 | 0.10 / 0.42 / 0.77 / 1.00 / 1.00 | 0.34 |
| phases | F_ST vs compact | 0.58 | 1.0 | 0.2 | 0.02 / 0.30 / 0.69 / 0.84 / 1.00 | 0.44 |
| phases | F_ST vs compact | 0.58 | 0.03 | 0.2 | 0.22 / 0.71 / 1.00 / 1.00 / 1.00 | 0.22 |
| phases | boundary excess vs random-center | 0.76 | 1.0 | 0.0 | 0.08 / 0.35 / 0.64 / 0.86 / 0.96 | 0.65 |
| phases | boundary excess vs random-center | 0.76 | 0.1 | 0.0 | 0.24 / 0.55 / 0.82 / 0.94 / 0.98 | 0.36 |
| phases | boundary excess vs random-center | 0.76 | 0.03 | 0.0 | 0.14 / 0.62 / 0.83 / 0.92 / 0.98 | 0.42 |
| phases | boundary excess vs random-center | 0.76 | 1.0 | 0.2 | 0.14 / 0.47 / 0.77 / 0.90 / 0.98 | 0.50 |
| phases | boundary excess vs random-center | 0.76 | 0.03 | 0.2 | 0.40 / 0.86 / 0.94 / 0.98 / 0.98 | 0.16 |
| phases | boundary excess vs compact | 0.84 | 1.0 | 0.0 | 0.02 / 0.30 / 0.61 / 0.84 / 1.00 | 0.78 |
| phases | boundary excess vs compact | 0.84 | 0.1 | 0.0 | 0.08 / 0.50 / 0.78 / 1.00 / 1.00 | 0.56 |
| phases | boundary excess vs compact | 0.84 | 0.03 | 0.0 | 0.22 / 0.62 / 0.84 / 0.95 / 1.00 | 0.56 |
| phases | boundary excess vs compact | 0.84 | 1.0 | 0.2 | 0.06 / 0.41 / 0.76 / 0.84 / 1.00 | 0.76 |
| phases | boundary excess vs compact | 0.84 | 0.03 | 0.2 | 0.44 / 0.78 / 0.98 / 1.00 / 1.00 | 0.45 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 1.0 | 0.0 | 0.14 / 0.40 / 0.58 / 0.72 / 0.86 | 0.52 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.1 | 0.0 | 0.22 / 0.42 / 0.66 / 0.84 / 0.90 | 0.42 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.03 | 0.0 | 0.28 / 0.47 / 0.68 / 0.86 / 0.90 | 0.39 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 1.0 | 0.2 | 0.26 / 0.59 / 0.70 / 0.86 / 0.90 | 0.25 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.03 | 0.2 | 0.18 / 0.62 / 0.84 / 0.90 / 0.90 | 0.22 |
| Parkin vs rest | F_ST vs compact | 1.00 | 1.0 | 0.0 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.1 | 0.0 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.03 | 0.0 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 1.0 | 0.2 | 0.00 / 0.62 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.03 | 0.2 | 0.00 / 0.62 / 1.00 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 1.0 | 0.0 | 0.00 / 0.09 / 0.42 / 0.72 / 0.86 | 0.50 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.1 | 0.0 | 0.00 / 0.18 / 0.48 / 0.84 / 0.90 | 0.44 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.03 | 0.0 | 0.00 / 0.15 / 0.48 / 0.84 / 0.90 | 0.46 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 1.0 | 0.2 | 0.00 / 0.20 / 0.55 / 0.86 / 0.90 | 0.40 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.03 | 0.2 | 0.00 / 0.18 / 0.71 / 0.90 / 0.90 | 0.38 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 1.0 | 0.0 | 0.00 / 0.00 / 0.38 / 1.00 / 1.00 | 0.59 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.1 | 0.0 | 0.00 / 0.00 / 0.38 / 1.00 / 1.00 | 0.52 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.03 | 0.0 | 0.00 / 0.00 / 0.38 / 1.00 / 1.00 | 0.57 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 1.0 | 0.2 | 0.00 / 0.00 / 0.62 / 1.00 / 1.00 | 0.43 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.03 | 0.2 | 0.00 / 0.00 / 1.00 / 1.00 / 1.00 | 0.40 |

### 2000 learners, innovation 0.002, mixing 0.005 (sensitivity combination)

| line | statistic | observed | copying factor | local innovation | 5% / 25% / 50% / 75% / 95% | records at or below observed |
|---|---|---:|---|---|---|---:|
| phases | F_ST vs random-center | 0.64 | 1.0 | 0.0 | 0.04 / 0.20 / 0.36 / 0.61 / 0.88 | 0.77 |
| phases | F_ST vs random-center | 0.64 | 0.1 | 0.0 | 0.10 / 0.34 / 0.60 / 0.78 / 0.94 | 0.57 |
| phases | F_ST vs random-center | 0.64 | 0.03 | 0.0 | 0.12 / 0.32 / 0.56 / 0.80 / 0.94 | 0.59 |
| phases | F_ST vs random-center | 0.64 | 1.0 | 0.2 | 0.06 / 0.43 / 0.72 / 0.92 / 0.98 | 0.42 |
| phases | F_ST vs random-center | 0.64 | 0.03 | 0.2 | 0.18 / 0.59 / 0.88 / 0.98 / 0.98 | 0.28 |
| phases | F_ST vs compact | 0.58 | 1.0 | 0.0 | 0.00 / 0.17 / 0.37 / 0.64 / 0.92 | 0.69 |
| phases | F_ST vs compact | 0.58 | 0.1 | 0.0 | 0.00 / 0.28 / 0.56 / 0.80 / 1.00 | 0.53 |
| phases | F_ST vs compact | 0.58 | 0.03 | 0.0 | 0.02 / 0.30 / 0.58 / 0.84 / 1.00 | 0.52 |
| phases | F_ST vs compact | 0.58 | 1.0 | 0.2 | 0.06 / 0.40 / 0.70 / 0.92 / 1.00 | 0.36 |
| phases | F_ST vs compact | 0.58 | 0.03 | 0.2 | 0.10 / 0.56 / 0.89 / 1.00 / 1.00 | 0.27 |
| phases | boundary excess vs random-center | 0.76 | 1.0 | 0.0 | 0.04 / 0.23 / 0.52 / 0.74 / 0.94 | 0.80 |
| phases | boundary excess vs random-center | 0.76 | 0.1 | 0.0 | 0.16 / 0.49 / 0.76 / 0.90 / 0.96 | 0.51 |
| phases | boundary excess vs random-center | 0.76 | 0.03 | 0.0 | 0.06 / 0.40 / 0.64 / 0.84 / 0.96 | 0.62 |
| phases | boundary excess vs random-center | 0.76 | 1.0 | 0.2 | 0.24 / 0.57 / 0.82 / 0.93 / 0.98 | 0.46 |
| phases | boundary excess vs random-center | 0.76 | 0.03 | 0.2 | 0.44 / 0.80 / 0.94 / 0.98 / 0.98 | 0.23 |
| phases | boundary excess vs compact | 0.84 | 1.0 | 0.0 | 0.00 / 0.22 / 0.47 / 0.76 / 1.00 | 0.90 |
| phases | boundary excess vs compact | 0.84 | 0.1 | 0.0 | 0.08 / 0.42 / 0.73 / 0.86 / 1.00 | 0.75 |
| phases | boundary excess vs compact | 0.84 | 0.03 | 0.0 | 0.00 / 0.28 / 0.66 / 0.84 / 1.00 | 0.79 |
| phases | boundary excess vs compact | 0.84 | 1.0 | 0.2 | 0.16 / 0.58 / 0.77 / 0.98 / 1.00 | 0.67 |
| phases | boundary excess vs compact | 0.84 | 0.03 | 0.2 | 0.34 / 0.75 / 0.92 / 1.00 / 1.00 | 0.48 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 1.0 | 0.0 | 0.06 / 0.32 / 0.52 / 0.69 / 0.86 | 0.61 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.1 | 0.0 | 0.04 / 0.34 / 0.62 / 0.81 / 0.90 | 0.46 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.03 | 0.0 | 0.08 / 0.38 / 0.54 / 0.74 / 0.90 | 0.54 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 1.0 | 0.2 | 0.08 / 0.54 / 0.71 / 0.86 / 0.90 | 0.27 |
| Parkin vs rest | F_ST vs random-center | 0.58 | 0.03 | 0.2 | 0.26 / 0.55 / 0.74 / 0.86 / 0.90 | 0.29 |
| Parkin vs rest | F_ST vs compact | 1.00 | 1.0 | 0.0 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.1 | 0.0 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.03 | 0.0 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 1.0 | 0.2 | 0.00 / 0.38 / 0.62 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | F_ST vs compact | 1.00 | 0.03 | 0.2 | 0.00 / 0.56 / 1.00 / 1.00 / 1.00 | 1.00 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 1.0 | 0.0 | 0.00 / 0.10 / 0.36 / 0.68 / 0.88 | 0.60 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.1 | 0.0 | 0.00 / 0.17 / 0.48 / 0.82 / 0.90 | 0.44 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.03 | 0.0 | 0.00 / 0.12 / 0.38 / 0.77 / 0.90 | 0.52 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 1.0 | 0.2 | 0.00 / 0.28 / 0.67 / 0.86 / 0.90 | 0.34 |
| Parkin vs rest | boundary excess vs random-center | 0.42 | 0.03 | 0.2 | 0.00 / 0.14 / 0.68 / 0.88 / 0.90 | 0.44 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 1.0 | 0.0 | 0.00 / 0.00 / 0.19 / 1.00 / 1.00 | 0.58 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.1 | 0.0 | 0.00 / 0.00 / 0.38 / 1.00 / 1.00 | 0.52 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.03 | 0.0 | 0.00 / 0.00 / 0.62 / 1.00 / 1.00 | 0.49 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 1.0 | 0.2 | 0.00 / 0.00 / 1.00 / 1.00 / 1.00 | 0.37 |
| Parkin vs rest | boundary excess vs compact | 0.38 | 0.03 | 0.2 | 0.00 / 0.00 / 0.81 / 1.00 / 1.00 | 0.45 |

## Reading

Compare the rows with a copying boundary (factor below 1) against the row without. If the tested line beats most alternatives only when copying is restricted at it, the comparison responds, and the observed record's failure to beat them weighs against a restriction there by as much as the at-or-below fractions below differ between cells, which is moderate at best. If the rows look alike, the comparison cannot see a copying boundary and the argument from it should be dropped. Read the medians together with the distributions section: a median that moves shows the comparison responds, and only the at-or-below fractions say how strongly the observed score weighs against a boundary.

Figure written to figures/figS12_partition_power.png and its siblings; every simulated record's shares are in `output/findings/partition_power_records.csv`.
