# Generator diagnostic: is 'only F_ST' a record property or a generator artifact?

Real config: n = 29, K = 10 types, k = 3 spatial clusters, 6 bins, rarefied to 50. 60 seeds. Entries are mean Spearman rho of the signature vs ordinal position (emergence should drive the channel's signature toward +1).

## Generator A: spatial divergence+conformity (F_ST channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.50 | -0.04 | -0.40 | +0.03 |
| 0.4 | +0.31 | -0.18 | +0.38 | -0.03 |
| 0.8 | -0.66 | -0.32 | +0.83 | -0.02 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+1.23). Per-signature response: neutral -1.15, seriability -0.28, fst +1.23, spatial -0.05.

## Generator B: single-pool conformity (neutral channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.51 | -0.07 | -0.40 | -0.05 |
| 0.4 | -0.77 | -0.02 | -0.32 | +0.08 |
| 0.8 | -0.58 | -0.07 | -0.41 | +0.08 |

- Largest response (rho at s=0.8 minus s=0): **spatial** (+0.13). Per-signature response: neutral -1.09, seriability -0.00, fst -0.01, spatial +0.13.

## Generator C: social non-spatial divergence (seriation channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.49 | -0.03 | -0.29 | -0.04 |
| 0.4 | +0.49 | -0.39 | -0.16 | -0.00 |
| 0.8 | +0.24 | -0.35 | +0.54 | -0.24 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+0.83). Per-signature response: neutral -0.25, seriability -0.32, fst +0.83, spatial -0.21.

## Rarefaction-depth sensitivity (generator A, s = 0.6)

| NRARE | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 50 | +0.12 | -0.33 | +0.70 | +0.06 |
| 100 | +0.15 | -0.29 | +0.79 | +0.05 |
| 150 | +0.21 | -0.27 | +0.86 | +0.15 |

Larger NRARE keeps more rare types; if neutral/seriation climb with depth, their weakness at NRARE=50 is partly a rarefaction-floor artifact, not a record property.

## Cluster-count sweep (generator A, s = 0.6)

| k | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 3 | +0.12 | -0.33 | +0.70 | +0.06 |
| 4 | +0.28 | -0.21 | +0.78 | -0.04 |
| 5 | +0.19 | -0.12 | +0.83 | +0.00 |

At k=3 the spatial boundary and F_ST summarize one partition; if the spatial signature climbs at k=4,5 its weakness is a low-cluster-count artifact.

## Verdict

Read the three generator tables together. If each signature shows its largest response to the generator whose channel it measures (A->F_ST, B->neutral, C->seriation), the four signatures are each empirically sufficient for a distinct emergence channel, and the four-signature criterion is vindicated in principle even though, for the spatial bounded-group emergence the Parkin hypothesis predicts, F_ST is the cleanest at this resolution. If B and C fail to light any signature, the signatures genuinely fail at this resolution regardless of channel, and the F_ST-led reading is a record property.