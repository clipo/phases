# Drift extension: joint tabulations (exploratory)

Produced by `.venv-np254/bin/python analyses/per_type/drift_extension_joint.py` from `output/findings/drift_extension_cells.csv`. Exploratory: the grid was run on 2026-11-14 without the Codex build audit, on one seed set of 50 records per cell, and these tabulations were chosen after the run was seen. They are not the registered reading rules, whose statements are in `output/findings/drift_extension.md`. Every count is a count of grid cells at which each named summary of the record lies inside its own 10th to 90th percentile band (marginal coverage). Cells that share a window share simulated histories, so the counts are not counts of independent confirmations.

## Primary point

- Level: all three common-class F values inside their bands at 5 of 504 evaluable cells (504 tested), all at F_W in {0.7, 0.9}. By span limit: Ld 0, Ld/4 2, 1 3. By window: 8LR 0, 2LR 0, 32LR 5.
- Rare classes in those cells: number of the seven with record occupancy inside its band, as cells per count: {0: 5} (maximum 0). In every one of those cells at least 7 of the seven have the record's occupancy above the band (the record has the class at more sites).
- Converse: 266 cells have at least 5 of the 7 rare classes inside; the largest number of common-class F values inside among them is 0. Among cells with at least one common-class F inside, the largest number of rare classes inside is 3. Cells with level and at least 5 rare classes inside: 0.
- Latitude trend of the painted types (A_lat, class 2) inside its band: dates independent of location (c0) 0 of 72 cells (record above the band at 72); half coupling (lat50N, lat50S) 13 of 144; strict ordering (lat100N, lat100S) 91 of 144.
- D_dist and the latitude trends of classes 0 and 2 all inside, by coupling: c0 0 of 72, lat50N 3 of 72, lat100N 43 of 72, lat50S 10 of 72, lat100S 31 of 72, stepL 11 of 72, stepS 5 of 72. Under strict ordering, by F_W: 0: 0, 0.01: 1, 0.03: 0, 0.1: 18, 0.3: 18, 0.5: 16, 0.7: 11, 0.9: 10 (of 18 each).
- Level and those three spatial readouts together: 0 cells.

Cells by number of common-class F values inside (rows) and number of rare classes inside (columns):

| F inside \ rare inside | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 4 | 16 | 76 | 88 | 45 | 133 |
| 1 | 8 | 2 | 7 | 3 | 0 | 0 | 0 | 0 |
| 2 | 55 | 50 | 11 | 1 | 0 | 0 | 0 | 0 |
| 3 | 5 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## Secondary point

- Level: all three common-class F values inside their bands at 45 of 496 evaluable cells (504 tested), all at F_W in {0.7, 0.9}. By span limit: Ld 17, Ld/4 14, 1 14. By window: 8LR 21, 2LR 0, 32LR 24.
- Rare classes in those cells: number of the seven with record occupancy inside its band, as cells per count: {0: 43, 1: 2} (maximum 1). In every one of those cells at least 6 of the seven have the record's occupancy above the band (the record has the class at more sites).
- Converse: 185 cells have at least 5 of the 7 rare classes inside; the largest number of common-class F values inside among them is 0. Among cells with at least one common-class F inside, the largest number of rare classes inside is 4. Cells with level and at least 5 rare classes inside: 0.
- Latitude trend of the painted types (A_lat, class 2) inside its band: dates independent of location (c0) 0 of 72 cells (record above the band at 70); half coupling (lat50N, lat50S) 1 of 144; strict ordering (lat100N, lat100S) 82 of 144.
- D_dist and the latitude trends of classes 0 and 2 all inside, by coupling: c0 0 of 72, lat50N 0 of 72, lat100N 38 of 72, lat50S 1 of 72, lat100S 38 of 72, stepL 5 of 72, stepS 0 of 72. Under strict ordering, by F_W: 0: 0, 0.01: 0, 0.03: 0, 0.1: 13, 0.3: 18, 0.5: 18, 0.7: 18, 0.9: 9 (of 18 each).
- Level and those three spatial readouts together: 13 cells, all under lat100N, lat100S; rare classes inside at most 0; status of the latitude trend of class 1 in them: {'below': 13}.

Cells by number of common-class F values inside (rows) and number of rare classes inside (columns):

| F inside \ rare inside | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 8 | 52 | 68 | 74 | 81 | 30 |
| 1 | 18 | 22 | 21 | 7 | 2 | 0 | 0 | 0 |
| 2 | 62 | 13 | 1 | 0 | 0 | 0 | 0 | 0 |
| 3 | 43 | 2 | 0 | 0 | 0 | 0 | 0 | 0 |
