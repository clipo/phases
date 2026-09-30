# Are the PFG phases real, or spatial-drift artifacts? (n = 28 curated assemblages)

Test of whether the within>between decorated-type similarity structure
underwriting the Parkin phase requires bounded social interaction, or is
what neutral drift on a spatially structured network produces with no
groupness (Lipo et al. 2021).

## Part A. Distance-decay strength
- Mantel r (BR similarity vs geographic distance): -0.322 (p = 0.0004).
- Geographic distance alone explains r^2 = 0.10 of the pairwise
  ceramic similarity variance. Negative r = closer assemblages are more
  similar, the isolation-by-distance signature.

## Part B. Excess beyond isolation-by-distance
- Ceramic communities (greedy modularity on BR graph): 2 communities, modularity Q = 0.045; sizes {0: 16, 1: 12}.
- Raw within-minus-between BR similarity gap: +41.6 BR units.
- Distance-CONTROLLED within-minus-between gap (boundary excess): +38.6 BR units.
  Under pure isolation-by-distance this collapses to ~0; a genuine
  interaction boundary leaves a positive excess.
- Partial Mantel, ceramic distance vs community membership controlling for geography: r = +0.606 (p = 0.0002).
- Parkin falls in community 0.

## Part C. Generative spatial-drift null (Lipo et al. 2021)
Neutral drift on the real coordinate layout, distance-decayed
interaction (length scale 12 km, between-node copy rate 0.05,
innovation 0.01), no imposed boundaries. 95% null interval over
200 seeds vs observed:

| statistic | observed | null mean | null 95% | obs inside null? |
|---|---|---|---|---|
| distance-decay r | -0.322 | -0.440 | [-0.619, -0.213] | yes |
| modularity Q | +0.045 | +0.172 | [+0.087, +0.279] | NO |
| boundary excess (BR) | +38.606 | +36.508 | [+17.835, +60.353] | yes |
| cultural F_ST | +0.039 | +0.103 | [+0.050, +0.183] | NO |

Interpretation: where observed falls INSIDE the spatial-drift null,
that aspect of the phase structure is reproduced by neutral drift on
geography alone and needs no social boundary. Where it falls OUTSIDE,
the data carry structure beyond spatially structured drift.
