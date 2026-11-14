# The two-cluster F_ST across every diversity-matched screen cell, without the F_ST selection

Produced by `analyses/probe_latent_groups/screen_fst_unselected.py` from `resemblance_screen.csv` (2,970 cells, six seeds each). Observed two-cluster F_ST 0.0062, compared at full precision. Descriptive; exploratory design validation, not evidence about the record.

## What the selection did

Of 2970 screen cells, 565 match the three diversity summaries within 10 percent and 498 of those have an accepted time-scale measurement. The stage-2 recheck took the 34 of these whose six-seed F_ST median lay within 50 percent of the observed value, and the operating points were chosen among them by the four-summary loss, whose fourth term is the squared relative error of that median (D22). By six-seed median the primary point ranks 25 of 498 matched cells and the secondary 11 of 498.

## Where the observed value lies in the unselected matched cells

| quantity | value |
|---|---|
| matched cells with accepted time scales | 498 |
| six-seed F_ST median: 10th, 25th, 50th, 75th, 90th percentile over cells | 0.0003, 0.0005, 0.0010, 0.0016, 0.0026 |
| cells whose six-seed median is at or above the observed value | 8 |
| cells whose six-seed 90th percentile is below the observed value (observed above that percentile) | 436 (88 percent) |
| cells with no seed of six at or above the observed value | 400 (80 percent) |
| mean fraction of seeds at or above the observed value, over cells | 0.050 |

## By copying-network family

| family | matched cells | median of six-seed medians | cells with observed above the 90th percentile | cells with median at or above observed |
|---|---|---|---|---|
| exponential-12km | 49 | 0.0009 | 48 | 0 |
| exponential-24km | 53 | 0.0008 | 51 | 0 |
| exponential-6km | 44 | 0.0013 | 36 | 0 |
| gabriel | 45 | 0.0009 | 41 | 0 |
| gaussian-12km | 51 | 0.0016 | 37 | 3 |
| gaussian-24km | 48 | 0.0010 | 45 | 0 |
| gaussian-6km | 51 | 0.0017 | 33 | 5 |
| power-12km | 54 | 0.0007 | 50 | 0 |
| power-24km | 57 | 0.0008 | 53 | 0 |
| power-6km | 46 | 0.0009 | 42 | 0 |

## The ten matched cells with the largest six-seed F_ST median

| family | learners | innovation | mixing | F_ST median | 90th percentile | seeds at or above observed |
|---|---|---|---|---|---|---|
| gaussian-12km | 2000 | 0.001 | 0.05 | 0.0096 | 0.0235 | 0.67 |
| gaussian-6km | 2000 | 0.001 | 0.4 | 0.0094 | 0.0195 | 0.67 |
| gaussian-6km | 2000 | 0.001 | 0.2 | 0.0074 | 0.0411 | 0.50 |
| gaussian-6km | 10000 | 0.0005 | 0.02 | 0.0069 | 0.0136 | 0.50 |
| gaussian-12km | 120 | 0.008 | 0.4 | 0.0069 | 0.0146 | 0.50 |
| gaussian-6km | 10000 | 0.0005 | 0.2 | 0.0066 | 0.0163 | 0.50 |
| gaussian-12km | 2000 | 0.0005 | 0.1 | 0.0066 | 0.0118 | 0.50 |
| gaussian-6km | 10000 | 0.0005 | 0.4 | 0.0064 | 0.0212 | 0.50 |
| gaussian-12km | 2000 | 0.002 | 0.005 | 0.0061 | 0.0135 | 0.50 |
| gaussian-12km | 2000 | 0.0005 | 0.05 | 0.0057 | 0.0087 | 0.50 |

Reading: a comparison of the observed two-cluster F_ST against the full set of diversity-matched cells, rather than against the cells selected for their closeness to it, places the observed value above the six-seed 90th percentile of most cells. The 50-seed bands at the two operating points (17 of 30 cells containing the observed value) describe cells selected for relatively high, observed-like F_ST, although neither point has the highest. Six seeds per cell resolve the position coarsely, so the counts above are stated as counts of cells, not as a probability.
