# Generator diagnostic: is 'only F_ST' a record property or a generator artifact?

Real config: n = 43, K = 10 types, k = 5 spatial clusters, 6 bins, rarefied to 50. 60 seeds. Entries are mean Spearman rho of the signature vs ordinal position (emergence should drive the channel's signature toward +1).

## Generator A: spatial divergence+conformity (F_ST channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.63 | -0.16 | -0.67 | +0.08 |
| 0.4 | +0.45 | -0.67 | +0.62 | -0.01 |
| 0.8 | -0.50 | -0.65 | +0.98 | -0.05 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+1.65). Per-signature response: neutral -1.13, seriability -0.49, fst +1.65, spatial -0.13.

## Generator B: single-pool conformity (neutral channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.63 | -0.15 | -0.71 | +0.03 |
| 0.4 | -0.85 | -0.17 | -0.68 | -0.03 |
| 0.8 | -0.55 | -0.07 | -0.69 | +0.01 |

- Largest response (rho at s=0.8 minus s=0): **seriability** (+0.08). Per-signature response: neutral -1.18, seriability +0.08, fst +0.02, spatial -0.02.

## Generator C: social non-spatial divergence (seriation channel)

| s | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 0.0 | +0.56 | -0.23 | -0.64 | -0.00 |
| 0.4 | +0.48 | -0.67 | -0.39 | +0.01 |
| 0.8 | +0.35 | -0.70 | +0.44 | -0.02 |

- Largest response (rho at s=0.8 minus s=0): **fst** (+1.08). Per-signature response: neutral -0.21, seriability -0.47, fst +1.08, spatial -0.02.

## Rarefaction-depth sensitivity (generator A, s = 0.6)

| NRARE | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 50 | +0.14 | -0.71 | +0.93 | -0.07 |
| 100 | +0.18 | -0.82 | +0.97 | -0.09 |
| 150 | +0.17 | -0.84 | +0.98 | -0.07 |

Larger NRARE keeps more rare types; if neutral/seriation climb with depth, their weakness at NRARE=50 is partly a rarefaction-floor artifact, not a record property.

## Cluster-count sweep (generator A, s = 0.6)

| k | neutral | seriability | fst | spatial |
|---|---|---|---|---|
| 3 | -0.36 | -0.73 | +0.97 | +0.04 |
| 4 | -0.22 | -0.67 | +0.94 | -0.10 |
| 5 | +0.14 | -0.71 | +0.93 | -0.07 |

At k=3 the spatial boundary and F_ST summarize one partition; if the spatial signature climbs at k=4,5 its weakness is a low-cluster-count artifact.

## Verdict

Read the three generator tables together. If each signature shows its largest response to the generator whose channel it measures (A->F_ST, B->neutral, C->seriation), the four signatures are each empirically sufficient for a distinct emergence channel, and the four-signature criterion is vindicated in principle even though, for the spatial bounded-group emergence the Parkin hypothesis predicts, F_ST is the cleanest at this resolution. If B and C fail to light any signature, the signatures genuinely fail at this resolution regardless of channel, and the F_ST-led reading is a record property.