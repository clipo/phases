# Spatial drift vs bounded groups: which best explains the basin?

Observed (n = 43 curated assemblages), two generative
models on the same coordinate layout, 150 seeds each. A model 'brackets'
the observed statistic if the observed value falls in its 95% envelope.

| statistic | observed | spatial-drift 95% | bounded-groups 95% | consistent with |
|---|---|---|---|---|
| distance-decay r | -0.360 | [-0.384, +0.011] | [-0.595, -0.447] | drift |
| modularity Q | +0.057 | [+0.045, +0.097] | [+0.212, +0.361] | drift |
| boundary excess (BR) | +48.009 | [+14.979, +46.272] | [+48.650, +93.811] | neither |
| cultural $F_{ST}$ | +0.040 | [+0.033, +0.139] | [+0.095, +0.308] | drift |

Reading: where the observed value sits inside the spatial-drift envelope
but outside the bounded-groups envelope, neutral drift on geography is
the better explanation. The bounded-groups model, by construction,
produces sharper community structure (higher Q, higher F_ST, larger
distance-controlled boundary excess) than the basin actually shows.
