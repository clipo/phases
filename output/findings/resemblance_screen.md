# Restricted resemblance screen under the probe's pipeline (pilot step 5)

Produced by `analyses/probe_latent_groups/resemblance_screen.py`; rules in `mls_emergence.probe.resemblance` (spec 3.3, R2-6). Exploratory design validation, not evidence about the record.

Grid: 10 copying-network families by 3 population sizes by 11 innovation rates by 9 mixing rates = 2970 cells, 6 seeds each, 17820 runs in 1254 s on 28 processes. Each cell at its own measured time scales (L_d from 2 to 2500, stored length 50 to 25000; 179 cells with an unaccepted measurement, flagged). Observed counts sha256 `41bd753fe1c7`; spatial clusters as in step 4 (2 clusters).

| summary | observed |
|---|---|
| hs | 0.5633 |
| rich | 6.2500 |
| ht | 0.5916 |
| fst_clusters | 0.0062 |

## Families: matched cells and the most spatially structured matched cell of each (D21)

| family | matched cells | best matched cell by F_ST median | L_d, S | mean hs | mean rich | mean ht | F_ST median | F_ST 90th pct | seeds at or above observed | loss |
|---|---|---|---|---|---|---|---|---|---|---|
| exponential-6km | 52 | n_ind 120, innovation 0.024, mixing 0.4 | 10, 360 | 0.5935 | 5.720 | 0.6063 | 0.0055 | 0.0180 | 0.50 | 0.0107 |
| exponential-12km | 56 | n_ind 2000, innovation 0.0005, mixing 0.02 | 103, 4040 | 0.5937 | 5.827 | 0.6070 | 0.0031 | 0.0095 | 0.33 | 0.0082 |
| exponential-24km | 62 | n_ind 120, innovation 0.012, mixing 0.1 | 12, 650 | 0.5850 | 5.655 | 0.6032 | 0.0052 | 0.0070 | 0.33 | 0.0109 |
| gaussian-6km | 56 | n_ind 10000, innovation 0.0002, mixing 0.05 | 2500, 25000 | 0.5636 | 5.869 | 0.5897 | 0.0135 | 0.0393 | 0.83 | 0.0037 |
| gaussian-12km | 57 | n_ind 2000, innovation 0.001, mixing 0.05 | 177, 3580 | 0.5809 | 5.679 | 0.5911 | 0.0096 | 0.0235 | 0.67 | 0.0093 |
| gaussian-24km | 54 | n_ind 10000, innovation 0.0002, mixing 0.001 | 2500, 25000 | 0.5555 | 6.095 | 0.5730 | 0.0058 | 0.0144 | 0.33 | 0.0018 |
| power-6km | 54 | n_ind 2000, innovation 0.008, mixing 0.001 | 76, 980 | 0.5848 | 5.988 | 0.6029 | 0.0029 | 0.0063 | 0.17 | 0.0036 |
| power-12km | 61 | n_ind 2000, innovation 0.002, mixing 0.01 | 80, 2880 | 0.5594 | 5.714 | 0.5754 | 0.0030 | 0.0060 | 0.17 | 0.0081 |
| power-24km | 61 | n_ind 2000, innovation 0.001, mixing 0.005 | 106, 3780 | 0.5877 | 5.720 | 0.6168 | 0.0050 | 0.0109 | 0.50 | 0.0109 |
| gabriel | 52 | n_ind 120, innovation 0.024, mixing 0.4 | 9, 390 | 0.5740 | 5.720 | 0.5831 | 0.0043 | 0.0070 | 0.33 | 0.0078 |

## Stop rule on the inherited operating point (exponential 24 km; 2,000; 0.001; 0.02)

Means across 6 seeds: hs 0.5839, rich 5.8214, ht 0.5961; relative errors hs 0.037, rich 0.069, ht 0.008 at tolerance 0.1; two-cluster F_ST median 0.0012 against observed 0.0062 (no tolerance). **Verdict: PASS.**

## Matched cells, all families, ten most structured

565 of 2970 cells match all three diversity summaries within 10 percent.

| cell | mean hs | mean rich | mean ht | F_ST median | seeds at or above observed | loss |
|---|---|---|---|---|---|---|
| gaussian-6km: n_ind 10000, innovation 0.0002, mixing 0.05 | 0.5636 | 5.869 | 0.5897 | 0.0135 | 0.83 | 0.0037 |
| gaussian-6km: n_ind 10000, innovation 0.0002, mixing 0.02 | 0.5776 | 6.143 | 0.5971 | 0.0118 | 0.83 | 0.0010 |
| gaussian-12km: n_ind 2000, innovation 0.001, mixing 0.05 | 0.5809 | 5.679 | 0.5911 | 0.0096 | 0.67 | 0.0093 |
| gaussian-6km: n_ind 2000, innovation 0.001, mixing 0.4 | 0.5739 | 5.744 | 0.5934 | 0.0094 | 0.67 | 0.0069 |
| gaussian-6km: n_ind 10000, innovation 0.0002, mixing 0.4 | 0.6097 | 5.863 | 0.6176 | 0.0082 | 0.83 | 0.0125 |
| gaussian-6km: n_ind 10000, innovation 0.0002, mixing 0.1 | 0.5620 | 5.762 | 0.5761 | 0.0079 | 0.67 | 0.0068 |
| gaussian-6km: n_ind 2000, innovation 0.001, mixing 0.2 | 0.5655 | 5.655 | 0.5867 | 0.0074 | 0.50 | 0.0092 |
| gaussian-6km: n_ind 10000, innovation 0.0005, mixing 0.02 | 0.5789 | 6.190 | 0.5890 | 0.0069 | 0.50 | 0.0009 |
| gaussian-12km: n_ind 120, innovation 0.008, mixing 0.4 | 0.5688 | 5.732 | 0.5847 | 0.0069 | 0.50 | 0.0071 |
| gaussian-6km: n_ind 10000, innovation 0.0005, mixing 0.2 | 0.5806 | 6.804 | 0.5897 | 0.0066 | 0.50 | 0.0088 |

## S0-stress cell within the inherited family (the stress cell is re-chosen once the operating family is decided)

**exponential-24km: n_ind 10000, innovation 0.0002, mixing 0.2** (grid index 205): the matched cell with the smallest F_ST median (0.0001; ties by loss 0.0057, then grid order).

Full table: `output/probe_latent_groups/resemblance_screen.csv`.
