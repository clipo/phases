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

## Reading

Compare the rows with a copying boundary (factor below 1) against the row without. If the tested line beats most alternatives only when copying is restricted at it, the comparison has power, and the observed record's failure to beat them is evidence against a restriction there. If the rows look alike, the comparison cannot see a copying boundary and the argument from it should be dropped.
