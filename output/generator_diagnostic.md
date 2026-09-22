# Generator diagnostic: is 'only F_ST' a record property or a generator artifact?

Real config: n = 28, K = 10 types, k = 2 spatial clusters, 6 bins, rarefied to 50. 60 seeds. Entries are mean Spearman rho of the signature vs ordinal position (emergence should drive the channel's signature toward +1).

## Generator A: spatial divergence+conformity (F_ST channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.05 | -0.06 | +0.01 | -0.02 |
| 0.4 | +0.00 | -0.17 | +0.54 | +0.03 |
| 0.8 | -0.69 | -0.15 | +0.96 | -0.01 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+0.95). Per-signature response: neutral -0.74, seriability -0.09, fst +0.95, spatial +0.01.

## Generator B: single-pool conformity (neutral channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.12 | -0.07 | +0.16 | -0.01 |
| 0.4 | -0.81 | -0.13 | +0.00 | +0.01 |
| 0.8 | -0.61 | -0.15 | -0.05 | -0.06 |

- Largest response (rho at s=0.8 minus s=0): **spatial** (-0.04). Per-signature response: neutral -0.73, seriability -0.08, fst -0.21, spatial -0.04.

## Generator C: social non-spatial divergence (seriation channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.15 | -0.14 | +0.07 | -0.06 |
| 0.4 | +0.03 | -0.35 | +0.36 | -0.04 |
| 0.8 | -0.14 | -0.35 | +0.77 | -0.11 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+0.69). Per-signature response: neutral -0.29, seriability -0.21, fst +0.69, spatial -0.05.

## Rarefaction-depth sensitivity (generator A, s = 0.6)

| NRARE | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 50 | -0.27 | -0.23 | +0.86 | +0.03 |
| 100 | -0.21 | -0.21 | +0.94 | -0.06 |
| 150 | -0.23 | -0.23 | +0.96 | -0.14 |

Larger NRARE keeps more rare types; if neutral/seriation climb with depth, their weakness at NRARE=50 is partly a rarefaction-floor artifact, not a record property.

## Cluster-count sweep (generator A, s = 0.6)

| k | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 3 | -0.25 | -0.25 | +0.86 | +0.11 |
| 4 | +0.02 | -0.10 | +0.85 | +0.12 |
| 5 | +0.08 | -0.14 | +0.91 | -0.04 |

At k=3 the spatial boundary and F_ST summarize one partition; if the spatial signature climbs at k=4,5 its weakness is a low-cluster-count artifact.

## Verdict

Read the three generator tables together. If each signature shows its largest response to the generator whose channel it measures (A->F_ST, B->neutral, C->seriation), the four signatures are each empirically sufficient for a distinct emergence channel, and the four-signature criterion is vindicated in principle even though, for the spatial bounded-group emergence the Parkin hypothesis predicts, F_ST is the cleanest at this resolution. If B and C fail to light any signature, the signatures genuinely fail at this resolution regardless of channel, and the F_ST-led reading is a record property.