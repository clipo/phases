# Generator diagnostic: is 'only F_ST' a record property or a generator artifact?

Real config: n = 29, K = 10 types, k = 3 spatial clusters, 6 bins, rarefied to 50. 60 seeds. Entries are mean Spearman rho of the signature vs ordinal position (emergence should drive the channel's signature toward +1).

## Generator A: spatial divergence+conformity (F_ST channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.42 | +0.09 | -0.51 | +0.03 |
| 0.4 | +0.33 | +0.21 | +0.11 | +0.04 |
| 0.8 | -0.62 | +0.31 | +0.76 | +0.17 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+1.27). Per-signature response: neutral -1.04, seriability +0.22, fst +1.27, spatial +0.14.

## Generator B: single-pool conformity (neutral channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.41 | +0.04 | -0.46 | +0.16 |
| 0.4 | -0.73 | +0.19 | -0.43 | +0.12 |
| 0.8 | -0.67 | +0.32 | -0.49 | +0.07 |

- Largest response (rho at s=0.8 minus s=0): **seriability** (+0.28). Per-signature response: neutral -1.08, seriability +0.28, fst -0.03, spatial -0.09.

## Generator C: social non-spatial divergence (seriation channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.48 | +0.18 | -0.48 | +0.01 |
| 0.4 | +0.39 | +0.07 | -0.28 | +0.08 |
| 0.8 | +0.27 | +0.12 | +0.29 | -0.21 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+0.77). Per-signature response: neutral -0.21, seriability -0.06, fst +0.77, spatial -0.23.

## Rarefaction-depth sensitivity (generator A, s = 0.6)

| NRARE | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 50 | +0.08 | +0.14 | +0.59 | +0.03 |
| 100 | +0.19 | +0.19 | +0.70 | +0.13 |
| 150 | +0.26 | +0.24 | +0.78 | +0.16 |

Larger NRARE keeps more rare types; if neutral/seriation climb with depth, their weakness at NRARE=50 is partly a rarefaction-floor artifact, not a record property.

## Cluster-count sweep (generator A, s = 0.6)

| k | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 3 | +0.08 | +0.14 | +0.59 | +0.03 |
| 4 | +0.29 | +0.09 | +0.69 | -0.02 |
| 5 | +0.21 | +0.17 | +0.68 | +0.19 |

At k=3 the spatial boundary and F_ST summarize one partition; if the spatial signature climbs at k=4,5 its weakness is a low-cluster-count artifact.

## Verdict

Read the three generator tables together. If each signature shows its largest response to the generator whose channel it measures (A->F_ST, B->neutral, C->seriation), the four signatures are each empirically sufficient for a distinct emergence channel, and the four-signature criterion is vindicated in principle even though, for the spatial bounded-group emergence the Parkin hypothesis predicts, F_ST is the cleanest at this resolution. If B and C fail to light any signature, the signatures genuinely fail at this resolution regardless of channel, and the F_ST-led reading is a record property.