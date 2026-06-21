# Spatial drift vs bounded groups: which best explains the basin?

Observed (n = 29 curated assemblages), two generative
models on the same coordinate layout, 150 seeds each. A model 'brackets'
the observed statistic if the observed value falls in its 95% envelope.

| statistic | observed | spatial-drift 95% | bounded-groups 95% | consistent with |
|---|---|---|---|---|
| distance-decay r | -0.301 | [-0.598, -0.149] | [-0.713, -0.595] | drift |
| modularity Q | +0.057 | [+0.040, +0.146] | [+0.321, +0.487] | drift |
| boundary excess (BR) | +43.438 | [+9.461, +46.247] | [+67.182, +115.506] | drift |
| cultural F_ST | +0.043 | [+0.027, +0.121] | [+0.107, +0.180] | drift |

Reading: where the observed value sits inside the spatial-drift envelope
but outside the bounded-groups envelope, neutral drift on geography is
the better explanation. The bounded-groups model, by construction,
produces sharper community structure (higher Q, higher F_ST, larger
distance-controlled boundary excess) than the basin actually shows.
