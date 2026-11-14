# Is the phases' ABC posterior a property of the phase lines or of the map?

Produced by `analyses/90_placebo_divisions.py`, which runs `analyses/84_phases_as_groups.py --placebo <seed>` for each seed: the phases are replaced by a division with their group sizes around random centers, and the whole grid, posterior and all, is rerun on it (300 runs per cell, 5 percent acceptance, prior 0.75 on both).

| groups | adjusted Rand with the phases | P(copying boundary at the lines) | P(local innovation) |
|---|---|---|---|
| the published phases | 1.000 | 0.75 | 0.68 |
| placebo, seed 11 | 0.264 | 0.75 | 0.70 |
| placebo, seed 12 | 0.171 | 0.70 | 0.66 |
| placebo, seed 13 | 0.346 | 0.85 | 0.80 |

If the placebo rows match or exceed the phases, the lean toward a copying boundary is what any division of this map produces, not evidence about the phase lines.
