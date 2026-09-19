# Does the Bayesian F_ST fit run on the basin? Measured.

Produced by `analyses/47_basin_scope_check.py`. Seed 0, 2000 draws x 4 chains, F ~ Uniform(0,1), clusters by the analysis-07 silhouette rule (k in 2..6, kmeans seed 7).

## The two fits

| | A: as analysis 43 runs today | B: restricted to the basin |
|---|---|---|
| assemblages | **55** | **29** |
| spatial clusters | 5 | 3 |
| BN F_ST median | **0.0665** | **0.0731** |
| BN 95% HDI | [0.0341, 0.1140] | [0.0220, 0.1579] |
| Gini-Simpson median | 0.0327 | 0.0179 |
| Gini-Simpson 95% HDI | [0.0301, 0.0352] | [0.0158, 0.0202] |
| plug-in F_ST | 0.0330 | 0.0179 |
| max R-hat | 1.0001 | 1.0004 |
| min ESS | 4672 | 4299 |
| divergences | 0 | 0 |

## What is in A but not in the basin

26 assemblages: `40LA007`, `40TP026`, `Beck`, `Bishop`, `Cheatham`, `Chuccalissa`, `Fullen`, `Graves_Lake`, `Hatchie`, `Irby`, `Jeter`, `Jones_Bayou`, `Lake_Cormorant`, `Mound_Place`, `Norfolk`, `Parchman`, `Porter`, `Pouncey`, `Rast`, `Richardsons_Landing`, `Salomon`, `Wall`, `Walls`, `Wilder`, `Woodlyn`, `Young`.

Silhouette scores on the basin coordinates, k = 2..6: k=2: 0.5528, k=3: 0.6008, k=4: 0.5051, k=5: 0.4935, k=6: 0.4774. k=3 wins clearly, which is the "three spatial clusters" the manuscript describes. The k=5 in column A is the silhouette optimum for the wider curated set, not for the basin.

## Reading

`prepare_inputs()` returns the whole curated set with coordinates. `43_bayesian_fst.basin_group_counts()` groups it by `inp.cluster_of` without restricting to `data/processed/basin_members_curated.txt`, so despite its name and despite the report headed "Observed St. Francis basin", column A is a regional quantity. The excluded 26 include `Walls`, `Wall`, `Chuccalissa` and `Parchman`, which are the Mississippi-proximal sites that `analyses/16_basin_membership.py` names in its own docstring as the reason the earlier latitude cut was replaced by the hydrological rule.

The Gini-Simpson readout, the quantity that matches the manuscript's estimator, moves from 0.0327 to 0.0179, a 45 percent reduction. The BN parameter moves less but its interval widens as it should on 29 assemblages rather than 55. Both fits are numerically healthy, so this is a scope defect and not a sampling one.

## Scope of the consequence

- **No manuscript number changes.** Figures 4 and 5 reach the data through `make_figures._load_curated()`, which does apply the basin membership. Only analyses 43 and 44 use the unrestricted path.
- `output/bayesian_fst.md` and the `STATUS.md` entry quoting "basin BN F = 0.067" and "Gini-Simpson readout 0.0327" attribute column A's numbers to the basin. Those are column B's questions with column A's answers.
- `analyses/44_bayesian_fst_validation.py` calibrates coverage, SBC and small-sample behavior at column A's geometry (5 clusters over 55) while describing it as the basin fit. The validation is not wrong, but it validates a design the paper does not use.

## Disposition

Restrict `basin_group_counts()` to the canonical membership and re-select k on the basin coordinates, then regenerate 43 and 44. Add a test that asserts the fitted set equals `data/processed/basin_members_curated.txt`, since a name alone did not prevent this. Scheduled as the first task of item 1.

## Why this matters beyond the numbers

This is the same defect class the 2026-06-10 alignment run found and fixed, where five figures were computed on the whole curated set while the text reported the basin. It recurred in code written a month later. The Verification Regime preamble in `CLAUDE.md` cites that episode as this project's own precedent for gates being necessary and never sufficient; the recurrence is the evidence that the lesson needed the test it is now getting.
