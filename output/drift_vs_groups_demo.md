# Spatial drift vs bounded groups: which best explains the basin?

Observed (n = 28 curated assemblages), two generative
models on the same coordinate layout, 150 seeds each. A model 'brackets'
the observed statistic if the observed value falls in its 95% envelope.

| statistic | observed | spatial-drift 95% | bounded-groups 95% | consistent with |
|---|---|---|---|---|
| distance-decay r | -0.321 | [-0.472, -0.044] | [-0.738, -0.625] | drift |
| modularity Q | +0.045 | [+0.034, +0.107] | [+0.299, +0.454] | drift |
| boundary excess (BR) | +38.636 | [+14.257, +40.291] | [+63.177, +113.806] | drift |
| cultural $F_{ST}$ | +0.039 | [+0.026, +0.123] | [+0.088, +0.181] | drift |

Reading: where the observed value sits inside the spatial-drift envelope
but outside the bounded-groups envelope, neutral drift on geography is
the better explanation. The bounded-groups model, by construction,
produces sharper community structure (higher Q, higher F_ST, larger
distance-controlled boundary excess) than the basin actually shows.
