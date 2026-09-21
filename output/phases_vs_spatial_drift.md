# Are the PFG phases real, or spatial-drift artifacts? (n = 43 curated assemblages)

Test of whether the within>between decorated-type similarity structure
underwriting the Parkin phase requires bounded social interaction, or is
what neutral drift on a spatially structured network produces with no
groupness (Lipo et al. 2021).

## Part A. Distance-decay strength
- Mantel r (BR similarity vs geographic distance): -0.360 (p = 0.0002).
- Geographic distance alone explains r^2 = 0.13 of the pairwise
  ceramic similarity variance. Negative r = closer assemblages are more
  similar, the isolation-by-distance signature.

## Part B. Excess beyond isolation-by-distance
- Ceramic communities (greedy modularity on BR graph): 2 communities, modularity Q = 0.057; sizes {0: 26, 1: 17}.
- Raw within-minus-between BR similarity gap: +47.0 BR units.
- Distance-CONTROLLED within-minus-between gap (boundary excess): +48.0 BR units.
  Under pure isolation-by-distance this collapses to ~0; a genuine
  interaction boundary leaves a positive excess.
- Partial Mantel, ceramic distance vs community membership controlling for geography: r = +0.613 (p = 0.0002).
- Parkin falls in community 0.

## Part C. Generative spatial-drift null (Lipo et al. 2021)
Neutral drift on the real coordinate layout, distance-decayed
interaction (length scale 12 km, between-node copy rate 0.05,
innovation 0.01), no imposed boundaries. 95% null interval over
200 seeds vs observed:

| statistic | observed | null mean | null 95% | obs inside null? |
|---|---|---|---|---|
| distance-decay r | -0.360 | -0.330 | [-0.495, -0.181] | yes |
| modularity Q | +0.057 | +0.191 | [+0.121, +0.294] | NO |
| boundary excess (BR) | +48.009 | +42.544 | [+27.009, +64.922] | yes |
| cultural F_ST | +0.040 | +0.147 | [+0.080, +0.247] | NO |

Interpretation: where observed falls INSIDE the spatial-drift null,
that aspect of the phase structure is reproduced by neutral drift on
geography alone and needs no social boundary. Where it falls OUTSIDE,
the data carry structure beyond spatially structured drift.
