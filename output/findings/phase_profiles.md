# How different are the phases, in sherds?

Produced by `analyses/99_phase_profiles.py`. 28 assemblages, phases Kent 11 (5,785 decorated sherds), Parkin 12 (7,398 decorated sherds), Walls 5 (918 decorated sherds). Model: the calibrated copying model with no groups, 300 runs with analysis 84's seeds (its no-boundary, regional-pool setting), sampled at the real sherd counts.

## Each phase's decorated sherds, percent by class

| class | Kent | Parkin | Walls | all |
|---|---:|---:|---:|---:|
| Parkin Punctated | 47.9 | 59.5 | 57.4 | 54.6 |
| Barton/Kent/MPI | 34.7 | 32.0 | 13.9 | 31.9 |
| Painted | 12.5 | 5.5 | 10.3 | 8.7 |
| Ranch Incised | 1.3 | 1.1 | 6.4 | 1.5 |
| Fortune Noded | 0.7 | 1.5 | 0.2 | 1.1 |
| Walls Engraved | 1.3 | 0.0 | 6.6 | 1.0 |
| Rhodes Incised | 0.5 | 0.1 | 3.5 | 0.5 |
| Wallace Incised | 0.7 | 0.1 | 0.0 | 0.3 |
| Vernon Paul Applique | 0.3 | 0.2 | 0.2 | 0.2 |
| Hull Engraved | 0.1 | 0.0 | 1.3 | 0.1 |

## The record against the model

| measure | observed | model median (95 percent range) | runs reaching observed |
|---|---:|---|---:|
| sherds that would have to change class to make two phases alike, mean over pairs (percent) | 19.4 | 7.4 (4.0 to 13.2) | 0 of 300 (0%) |
| range across phases in the share of the commonest class, Parkin Punctated (percentage points) | 11.7 | 6.7 (1.3 to 16.9) | 36 of 300 (12%) |
| between-phase cultural F_ST | 0.0123 | 0.0024 (0.0006 to 0.0116) | 7 of 300 (2%) |

## Pair by pair

The share of sherds that would have to change class to make the two phases' profiles alike.

| pair | observed (percent) | model median (95 percent range) | runs reaching observed |
|---|---:|---|---:|
| Kent and Parkin | 12.5 | 6.2 (2.5 to 14.0) | 17 of 300 (6%) |
| Kent and Walls | 24.2 | 6.5 (2.5 to 12.9) | 0 of 300 (0%) |
| Parkin and Walls | 21.6 | 9.0 (3.8 to 19.2) | 3 of 300 (1%) |

Straight-line distance between the nearest assemblages of two phases: Kent and Parkin, 12 km; Kent and Walls, 16 km; Parkin and Walls, 25 km.

## Each phase against the rest, beside same-size divisions

Analysis 86's comparison (2,000 posterior draws of the counts, 300 alternative divisions of each kind), with one phase against the other two taken together, at the seed analysis 86 uses for the Parkin phase, so the Parkin row is that analysis's own result.

| division | assemblages | F_ST, posterior median (95% CrI) | median F_ST, random-center / compact divisions | P(separates better than random-center / compact) | boundary excess | P(larger than random-center / compact) |
|---|---|---|---|---|---|---|
| Kent against the rest | 11 of 28 | 0.0077 (0.0058 to 0.0100) | 0.0064 / 0.0076 | 0.62 / 0.49 | +6.5 | 0.46 / 0.28 |
| Parkin against the rest | 12 of 28 | 0.0066 (0.0052 to 0.0083) | 0.0062 / 0.0057 | 0.52 / 0.94 | +5.9 | 0.43 / 0.49 |
| Walls against the rest | 5 of 28 | 0.0047 (0.0037 to 0.0057) | 0.0048 / 0.0048 | 0.43 / 0.47 | +9.1 | 0.54 / 0.54 |

## Reading

The three measures say the same thing in different units. In a run, the commonest class is whichever class is commonest in that run, so the second row compares like with like. The F_ST row repeats Table 1's comparison and is here so the three can be read side by side; it should equal analysis 84's no-boundary row.

The last table asks a different question from the first two. A phase can differ from its neighbors by more than the copying model produces and still be no more distinct than any other group of neighboring assemblages of its size; a probability near one half in that table says exactly that.
