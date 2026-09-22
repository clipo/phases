# Spatial drift vs bounded groups: which best explains the basin?

Observed (n = 28 curated assemblages), two generative
models on the same coordinate layout, 150 seeds each. A model 'brackets'
the observed statistic if the observed value falls in its 95% envelope.

| statistic | observed | spatial-drift 95% | bounded-groups 95% | consistent with |
|---|---|---|---|---|
| distance-decay r | -0.322 | [-0.486, -0.067] | [-0.741, -0.626] | drift |
| modularity Q | +0.045 | [+0.034, +0.106] | [+0.265, +0.458] | drift |
| boundary excess (BR) | +38.636 | [+12.897, +39.874] | [+59.277, +113.417] | drift |
| cultural $F_{ST}$ | +0.039 | [+0.028, +0.112] | [+0.082, +0.202] | drift |

Reading: where the observed value sits inside the spatial-drift envelope
but outside the bounded-groups envelope, neutral drift on geography is
the better explanation. The bounded-groups model, by construction,
produces sharper community structure (higher Q, higher F_ST, larger
distance-controlled boundary excess) than the basin actually shows.
