# Span variants at the operating point: where the between-cluster contrast goes (pilot step 5)

Produced by `analyses/probe_latent_groups/resemblance_variants.py`, 20 seeds per variant, same functions as the screen. Exploratory design validation, not evidence about the record.

Observed: hs 0.5633, rich 6.250, ht 0.5916, two-cluster F_ST 0.0062. Operating point {'n_ind': 2000, 'innovation': 0.001, 'mixing': 0.02}, stored length 3890, L_d 78, burn-in 5000.

| variant | mean hs | mean rich | mean ht | F_ST median [10%, 90%] |
|---|---|---|---|---|
| T1 (screen): lengths 1..L_d, ends over last L_d | 0.5848 | 6.321 | 0.5927 | 0.0006 [0.0002, 0.0021] |
| T0: same lengths, common end | 0.5821 | 6.302 | 0.5899 | 0.0007 [0.0002, 0.0021] |
| single step, no averaging | 0.5834 | 6.018 | 0.5916 | 0.0012 [0.0004, 0.0025] |
| eight snapshots 20 steps apart (old window), common end | 0.5890 | 6.609 | 0.5931 | 0.0004 [0.0002, 0.0019] |
| last step, expected counts (no multinomial noise) | 0.5828 | 6.264 | 0.5905 | 0.0012 [0.0005, 0.0025] |
