# Could the phase-line comparison see a copying boundary if there were one?

Produced by `analyses/89_partition_power.py`: 100 simulated records per cell, rebuilt with analysis 84's seeds at its calibrated cell (N 2000, innovation 0.001, mixing 0.02, 24 km along the rivers); 50 same-size alternative divisions of each kind. Each entry is the median over simulated records of the share of alternatives the tested line beats, with the share of records in which it beats at least 90 percent of them in parentheses.

| copying factor at the phase lines | local innovation | line tested | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |
|---|---|---|---|---|---|---|
| 1.0 | 0.0 | phases | 0.62 (9%) | 0.41 (0%) | 0.56 (12%) | 0.27 (0%) |
| 1.0 | 0.0 | Parkin vs rest | 0.64 (0%) | 0.62 (0%) | 0.40 (0%) | 0.00 (0%) |
| 0.1 | 0.0 | phases | 0.80 (27%) | 0.58 (0%) | 0.76 (28%) | 0.58 (0%) |
| 0.1 | 0.0 | Parkin vs rest | 0.72 (0%) | 0.62 (0%) | 0.58 (0%) | 0.62 (0%) |
| 0.03 | 0.0 | phases | 0.90 (52%) | 0.74 (0%) | 0.90 (53%) | 0.78 (0%) |
| 0.03 | 0.0 | Parkin vs rest | 0.72 (0%) | 0.62 (0%) | 0.72 (0%) | 0.62 (0%) |
| 1.0 | 0.2 | phases | 0.76 (27%) | 0.58 (0%) | 0.72 (24%) | 0.64 (0%) |
| 1.0 | 0.2 | Parkin vs rest | 0.72 (0%) | 0.62 (0%) | 0.64 (0%) | 0.62 (0%) |
| 0.03 | 0.2 | phases | 0.94 (74%) | 0.80 (0%) | 0.94 (82%) | 0.83 (0%) |
| 0.03 | 0.2 | Parkin vs rest | 0.74 (0%) | 0.62 (0%) | 0.73 (0%) | 0.62 (0%) |

The observed record (plug-in, these alternatives):

| line | F_ST vs random-center | F_ST vs compact | boundary excess vs random-center | boundary excess vs compact |
|---|---|---|---|---|
| phases | 0.30 | 0.04 | 0.72 | 0.68 |
| Parkin vs rest | 0.28 | 0.00 | 0.62 | 0.00 |

## The distributions behind the medians

For each cell, the 5th, 25th, 50th, 75th and 95th percentiles of the share of alternatives the tested line beats, over simulated records, and the fraction of simulated records whose share is at or below the observed record's. A fraction of 0 means none of the simulated records scored that low, which bounds the frequency below about 1/reps, not at zero.

### 2000 learners, innovation 0.001, mixing 0.02 (calibrated combination)

| line | statistic | observed | copying factor | local innovation | 5% / 25% / 50% / 75% / 95% | records at or below observed |
|---|---|---:|---|---|---|---:|
| phases | F_ST vs random-center | 0.30 | 1.0 | 0.0 | 0.32 / 0.46 / 0.62 / 0.82 / 0.92 | 0.05 |
| phases | F_ST vs random-center | 0.30 | 0.1 | 0.0 | 0.38 / 0.64 / 0.80 / 0.90 / 0.96 | 0.02 |
| phases | F_ST vs random-center | 0.30 | 0.03 | 0.0 | 0.44 / 0.76 / 0.90 / 0.93 / 0.96 | 0.01 |
| phases | F_ST vs random-center | 0.30 | 1.0 | 0.2 | 0.34 / 0.62 / 0.76 / 0.90 / 0.96 | 0.05 |
| phases | F_ST vs random-center | 0.30 | 0.03 | 0.2 | 0.66 / 0.88 / 0.94 / 0.96 / 0.96 | 0.01 |
| phases | F_ST vs compact | 0.04 | 1.0 | 0.0 | 0.02 / 0.14 / 0.41 / 0.58 / 0.80 | 0.13 |
| phases | F_ST vs compact | 0.04 | 0.1 | 0.0 | 0.00 / 0.28 / 0.58 / 0.72 / 0.84 | 0.14 |
| phases | F_ST vs compact | 0.04 | 0.03 | 0.0 | 0.14 / 0.57 / 0.74 / 0.84 / 0.84 | 0.04 |
| phases | F_ST vs compact | 0.04 | 1.0 | 0.2 | 0.00 / 0.32 / 0.58 / 0.70 / 0.84 | 0.09 |
| phases | F_ST vs compact | 0.04 | 0.03 | 0.2 | 0.42 / 0.68 / 0.80 / 0.84 / 0.84 | 0.01 |
| phases | boundary excess vs random-center | 0.72 | 1.0 | 0.0 | 0.02 / 0.21 / 0.56 / 0.80 / 0.92 | 0.71 |
| phases | boundary excess vs random-center | 0.72 | 0.1 | 0.0 | 0.16 / 0.43 / 0.76 / 0.90 / 0.96 | 0.45 |
| phases | boundary excess vs random-center | 0.72 | 0.03 | 0.0 | 0.40 / 0.81 / 0.90 / 0.94 / 0.96 | 0.18 |
| phases | boundary excess vs random-center | 0.72 | 1.0 | 0.2 | 0.10 / 0.46 / 0.72 / 0.88 / 0.96 | 0.52 |
| phases | boundary excess vs random-center | 0.72 | 0.03 | 0.2 | 0.69 / 0.92 / 0.94 / 0.96 / 0.96 | 0.06 |
| phases | boundary excess vs compact | 0.68 | 1.0 | 0.0 | 0.00 / 0.06 / 0.27 / 0.66 / 0.84 | 0.77 |
| phases | boundary excess vs compact | 0.68 | 0.1 | 0.0 | 0.04 / 0.20 / 0.58 / 0.80 / 0.84 | 0.59 |
| phases | boundary excess vs compact | 0.68 | 0.03 | 0.0 | 0.20 / 0.66 / 0.78 / 0.84 / 0.84 | 0.32 |
| phases | boundary excess vs compact | 0.68 | 1.0 | 0.2 | 0.00 / 0.20 / 0.64 / 0.80 / 0.84 | 0.60 |
| phases | boundary excess vs compact | 0.68 | 0.03 | 0.2 | 0.46 / 0.70 / 0.83 / 0.84 / 0.84 | 0.23 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 1.0 | 0.0 | 0.12 / 0.40 / 0.64 / 0.72 / 0.80 | 0.16 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.1 | 0.0 | 0.08 / 0.56 / 0.72 / 0.78 / 0.80 | 0.14 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.03 | 0.0 | 0.14 / 0.49 / 0.72 / 0.80 / 0.80 | 0.10 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 1.0 | 0.2 | 0.10 / 0.61 / 0.72 / 0.80 / 0.80 | 0.10 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.03 | 0.2 | 0.18 / 0.53 / 0.74 / 0.80 / 0.80 | 0.08 |
| Parkin vs rest | F_ST vs compact | 0.00 | 1.0 | 0.0 | 0.00 / 0.46 / 0.62 / 0.62 / 0.62 | 0.25 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.1 | 0.0 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.17 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.03 | 0.0 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.20 |
| Parkin vs rest | F_ST vs compact | 0.00 | 1.0 | 0.2 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.14 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.03 | 0.2 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.19 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 1.0 | 0.0 | 0.00 / 0.11 / 0.40 / 0.72 / 0.80 | 0.67 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.1 | 0.0 | 0.00 / 0.18 / 0.58 / 0.79 / 0.80 | 0.54 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.03 | 0.0 | 0.00 / 0.45 / 0.72 / 0.80 / 0.80 | 0.34 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 1.0 | 0.2 | 0.00 / 0.11 / 0.64 / 0.80 / 0.80 | 0.50 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.03 | 0.2 | 0.04 / 0.50 / 0.73 / 0.80 / 0.80 | 0.33 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 1.0 | 0.0 | 0.00 / 0.00 / 0.00 / 0.62 / 0.62 | 0.54 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.1 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.37 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.03 | 0.0 | 0.00 / 0.46 / 0.62 / 0.62 / 0.62 | 0.25 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 1.0 | 0.2 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.40 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.03 | 0.2 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.18 |

### 10000 learners, innovation 0.0005, mixing 0.002 (sensitivity combination)

| line | statistic | observed | copying factor | local innovation | 5% / 25% / 50% / 75% / 95% | records at or below observed |
|---|---|---:|---|---|---|---:|
| phases | F_ST vs random-center | 0.30 | 1.0 | 0.0 | 0.26 / 0.46 / 0.63 / 0.78 / 0.90 | 0.07 |
| phases | F_ST vs random-center | 0.30 | 0.1 | 0.0 | 0.28 / 0.56 / 0.69 / 0.84 / 0.92 | 0.09 |
| phases | F_ST vs random-center | 0.30 | 0.03 | 0.0 | 0.38 / 0.68 / 0.82 / 0.88 / 0.96 | 0.03 |
| phases | F_ST vs random-center | 0.30 | 1.0 | 0.2 | 0.40 / 0.66 / 0.83 / 0.92 / 0.96 | 0.01 |
| phases | F_ST vs random-center | 0.30 | 0.03 | 0.2 | 0.54 / 0.76 / 0.92 / 0.96 / 0.96 | 0.02 |
| phases | F_ST vs compact | 0.04 | 1.0 | 0.0 | 0.06 / 0.28 / 0.48 / 0.67 / 0.84 | 0.05 |
| phases | F_ST vs compact | 0.04 | 0.1 | 0.0 | 0.14 / 0.35 / 0.59 / 0.70 / 0.84 | 0.04 |
| phases | F_ST vs compact | 0.04 | 0.03 | 0.0 | 0.12 / 0.46 / 0.66 / 0.80 / 0.84 | 0.03 |
| phases | F_ST vs compact | 0.04 | 1.0 | 0.2 | 0.14 / 0.44 / 0.62 / 0.80 / 0.84 | 0.04 |
| phases | F_ST vs compact | 0.04 | 0.03 | 0.2 | 0.30 / 0.68 / 0.80 / 0.84 / 0.84 | 0.00 |
| phases | boundary excess vs random-center | 0.72 | 1.0 | 0.0 | 0.02 / 0.22 / 0.49 / 0.76 / 0.94 | 0.68 |
| phases | boundary excess vs random-center | 0.72 | 0.1 | 0.0 | 0.10 / 0.47 / 0.74 / 0.90 / 0.96 | 0.48 |
| phases | boundary excess vs random-center | 0.72 | 0.03 | 0.0 | 0.10 / 0.48 / 0.79 / 0.92 / 0.96 | 0.43 |
| phases | boundary excess vs random-center | 0.72 | 1.0 | 0.2 | 0.08 / 0.48 / 0.77 / 0.92 / 0.96 | 0.45 |
| phases | boundary excess vs random-center | 0.72 | 0.03 | 0.2 | 0.28 / 0.78 / 0.92 / 0.96 / 0.96 | 0.18 |
| phases | boundary excess vs compact | 0.68 | 1.0 | 0.0 | 0.00 / 0.16 / 0.41 / 0.67 / 0.84 | 0.78 |
| phases | boundary excess vs compact | 0.68 | 0.1 | 0.0 | 0.06 / 0.32 / 0.69 / 0.80 / 0.84 | 0.50 |
| phases | boundary excess vs compact | 0.68 | 0.03 | 0.0 | 0.06 / 0.28 / 0.67 / 0.82 / 0.84 | 0.51 |
| phases | boundary excess vs compact | 0.68 | 1.0 | 0.2 | 0.04 / 0.20 / 0.63 / 0.80 / 0.84 | 0.58 |
| phases | boundary excess vs compact | 0.68 | 0.03 | 0.2 | 0.20 / 0.70 / 0.82 / 0.84 / 0.84 | 0.23 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 1.0 | 0.0 | 0.16 / 0.40 / 0.63 / 0.72 / 0.80 | 0.14 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.1 | 0.0 | 0.12 / 0.44 / 0.64 / 0.72 / 0.80 | 0.17 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.03 | 0.0 | 0.14 / 0.35 / 0.62 / 0.72 / 0.80 | 0.19 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 1.0 | 0.2 | 0.26 / 0.53 / 0.66 / 0.80 / 0.80 | 0.07 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.03 | 0.2 | 0.32 / 0.53 / 0.72 / 0.80 / 0.80 | 0.05 |
| Parkin vs rest | F_ST vs compact | 0.00 | 1.0 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.29 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.1 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.29 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.03 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.29 |
| Parkin vs rest | F_ST vs compact | 0.00 | 1.0 | 0.2 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.17 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.03 | 0.2 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.18 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 1.0 | 0.0 | 0.00 / 0.00 / 0.39 / 0.70 / 0.80 | 0.71 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.1 | 0.0 | 0.00 / 0.20 / 0.56 / 0.80 / 0.80 | 0.59 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.03 | 0.0 | 0.00 / 0.18 / 0.54 / 0.76 / 0.80 | 0.58 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 1.0 | 0.2 | 0.00 / 0.14 / 0.48 / 0.80 / 0.80 | 0.64 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.03 | 0.2 | 0.00 / 0.41 / 0.69 / 0.80 / 0.80 | 0.38 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 1.0 | 0.0 | 0.00 / 0.00 / 0.00 / 0.62 / 0.62 | 0.56 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.1 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.44 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.03 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.46 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 1.0 | 0.2 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.47 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.03 | 0.2 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.21 |

### 2000 learners, innovation 0.002, mixing 0.005 (sensitivity combination)

| line | statistic | observed | copying factor | local innovation | 5% / 25% / 50% / 75% / 95% | records at or below observed |
|---|---|---:|---|---|---|---:|
| phases | F_ST vs random-center | 0.30 | 1.0 | 0.0 | 0.18 / 0.45 / 0.60 / 0.76 / 0.88 | 0.12 |
| phases | F_ST vs random-center | 0.30 | 0.1 | 0.0 | 0.28 / 0.52 / 0.66 / 0.82 / 0.94 | 0.06 |
| phases | F_ST vs random-center | 0.30 | 0.03 | 0.0 | 0.26 / 0.47 / 0.67 / 0.84 / 0.94 | 0.09 |
| phases | F_ST vs random-center | 0.30 | 1.0 | 0.2 | 0.36 / 0.72 / 0.84 / 0.92 / 0.96 | 0.04 |
| phases | F_ST vs random-center | 0.30 | 0.03 | 0.2 | 0.36 / 0.70 / 0.87 / 0.96 / 0.96 | 0.03 |
| phases | F_ST vs compact | 0.04 | 1.0 | 0.0 | 0.00 / 0.23 / 0.42 / 0.62 / 0.80 | 0.10 |
| phases | F_ST vs compact | 0.04 | 0.1 | 0.0 | 0.06 / 0.28 / 0.49 / 0.70 / 0.84 | 0.05 |
| phases | F_ST vs compact | 0.04 | 0.03 | 0.0 | 0.06 / 0.22 / 0.50 / 0.71 / 0.84 | 0.05 |
| phases | F_ST vs compact | 0.04 | 1.0 | 0.2 | 0.14 / 0.48 / 0.66 / 0.81 / 0.84 | 0.04 |
| phases | F_ST vs compact | 0.04 | 0.03 | 0.2 | 0.06 / 0.53 / 0.71 / 0.84 / 0.84 | 0.05 |
| phases | boundary excess vs random-center | 0.72 | 1.0 | 0.0 | 0.10 / 0.22 / 0.50 / 0.81 / 0.92 | 0.67 |
| phases | boundary excess vs random-center | 0.72 | 0.1 | 0.0 | 0.02 / 0.26 / 0.65 / 0.86 / 0.94 | 0.55 |
| phases | boundary excess vs random-center | 0.72 | 0.03 | 0.0 | 0.02 / 0.21 / 0.63 / 0.88 / 0.96 | 0.56 |
| phases | boundary excess vs random-center | 0.72 | 1.0 | 0.2 | 0.04 / 0.56 / 0.81 / 0.92 / 0.96 | 0.38 |
| phases | boundary excess vs random-center | 0.72 | 0.03 | 0.2 | 0.20 / 0.65 / 0.92 / 0.96 / 0.96 | 0.30 |
| phases | boundary excess vs compact | 0.68 | 1.0 | 0.0 | 0.02 / 0.18 / 0.46 / 0.68 / 0.84 | 0.77 |
| phases | boundary excess vs compact | 0.68 | 0.1 | 0.0 | 0.00 / 0.17 / 0.53 / 0.78 / 0.84 | 0.67 |
| phases | boundary excess vs compact | 0.68 | 0.03 | 0.0 | 0.00 / 0.14 / 0.30 / 0.72 / 0.84 | 0.64 |
| phases | boundary excess vs compact | 0.68 | 1.0 | 0.2 | 0.04 / 0.29 / 0.67 / 0.82 / 0.84 | 0.53 |
| phases | boundary excess vs compact | 0.68 | 0.03 | 0.2 | 0.06 / 0.46 / 0.80 / 0.84 / 0.84 | 0.37 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 1.0 | 0.0 | 0.08 / 0.26 / 0.56 / 0.72 / 0.80 | 0.29 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.1 | 0.0 | 0.08 / 0.30 / 0.56 / 0.72 / 0.80 | 0.24 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.03 | 0.0 | 0.14 / 0.32 / 0.55 / 0.72 / 0.80 | 0.23 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 1.0 | 0.2 | 0.15 / 0.49 / 0.66 / 0.79 / 0.80 | 0.09 |
| Parkin vs rest | F_ST vs random-center | 0.28 | 0.03 | 0.2 | 0.14 / 0.40 / 0.66 / 0.78 / 0.80 | 0.09 |
| Parkin vs rest | F_ST vs compact | 0.00 | 1.0 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.41 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.1 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.34 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.03 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.38 |
| Parkin vs rest | F_ST vs compact | 0.00 | 1.0 | 0.2 | 0.00 / 0.62 / 0.62 / 0.62 / 0.62 | 0.21 |
| Parkin vs rest | F_ST vs compact | 0.00 | 0.03 | 0.2 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.27 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 1.0 | 0.0 | 0.00 / 0.08 / 0.38 / 0.72 / 0.80 | 0.67 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.1 | 0.0 | 0.00 / 0.01 / 0.51 / 0.76 / 0.80 | 0.61 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.03 | 0.0 | 0.00 / 0.00 / 0.41 / 0.68 / 0.80 | 0.70 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 1.0 | 0.2 | 0.00 / 0.22 / 0.57 / 0.80 / 0.80 | 0.57 |
| Parkin vs rest | boundary excess vs random-center | 0.62 | 0.03 | 0.2 | 0.00 / 0.42 / 0.66 / 0.80 / 0.80 | 0.46 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 1.0 | 0.0 | 0.00 / 0.00 / 0.00 / 0.62 / 0.62 | 0.57 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.1 | 0.0 | 0.00 / 0.00 / 0.00 / 0.62 / 0.62 | 0.51 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.03 | 0.0 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.49 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 1.0 | 0.2 | 0.00 / 0.00 / 0.62 / 0.62 / 0.62 | 0.39 |
| Parkin vs rest | boundary excess vs compact | 0.00 | 0.03 | 0.2 | 0.00 / 0.46 / 0.62 / 0.62 / 0.62 | 0.25 |

## Reading

Compare the rows with a copying boundary (factor below 1) against the row without. If the tested line beats most alternatives only when copying is restricted at it, the comparison has power, and the observed record's failure to beat them is evidence against a restriction there. If the rows look alike, the comparison cannot see a copying boundary and the argument from it should be dropped. Read the medians together with the distributions section: a median that moves shows the comparison responds, and only the at-or-below fractions say how strongly the observed score weighs against a boundary.

Figure written to figures/figS12_partition_power.png and its siblings; every simulated record's shares are in `output/findings/partition_power_records.csv`.
