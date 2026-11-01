# Spatial drift vs bounded groups: which best explains the basin?

Observed (n = 28 curated assemblages), two generative
models on the same coordinate layout, 150 seeds each. A model 'brackets'
the observed statistic if the observed value falls in its 95% envelope.

| statistic | observed | spatial-drift 95% | bounded-groups 95% | consistent with |
|---|---|---|---|---|
| distance-decay r | -0.340 | [-0.495, -0.093] | [-0.778, -0.654] | drift |
| modularity Q | +0.045 | [+0.038, +0.109] | [+0.259, +0.457] | drift |
| boundary excess (BR) | +34.500 | [+15.243, +36.454] | [+56.929, +108.088] | drift |
| cultural $F_{ST}$ | +0.039 | [+0.028, +0.104] | [+0.095, +0.200] | drift |

Reading: where the observed value sits inside the spatial-drift envelope
but outside the bounded-groups envelope, neutral drift on geography is
the better explanation. The bounded-groups model, by construction,
produces sharper community structure (higher Q, higher F_ST, larger
distance-controlled boundary excess) than the basin actually shows.
