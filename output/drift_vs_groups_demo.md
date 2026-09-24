# Spatial drift vs bounded groups: which best explains the basin?

Observed (n = 28 curated assemblages), two generative
models on the same coordinate layout, 150 seeds each. A model 'brackets'
the observed statistic if the observed value falls in its 95% envelope.

| statistic | observed | spatial-drift 95% | bounded-groups 95% | consistent with |
|---|---|---|---|---|
| distance-decay r | -0.322 | [-0.450, -0.046] | [-0.737, -0.632] | drift |
| modularity Q | +0.045 | [+0.034, +0.102] | [+0.289, +0.451] | drift |
| boundary excess (BR) | +38.606 | [+10.672, +41.533] | [+62.299, +112.176] | drift |
| cultural $F_{ST}$ | +0.039 | [+0.022, +0.121] | [+0.087, +0.186] | drift |

Reading: where the observed value sits inside the spatial-drift envelope
but outside the bounded-groups envelope, neutral drift on geography is
the better explanation. The bounded-groups model, by construction,
produces sharper community structure (higher Q, higher F_ST, larger
distance-controlled boundary excess) than the basin actually shows.
