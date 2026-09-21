# Does the Bayesian F_ST fit run on the basin? Measured.

Produced by `analyses/47_basin_scope_check.py`. Seed 0, 2000 draws x 4 chains, F ~ Uniform(0,1), clusters by the analysis-07 silhouette rule (k in 2..6, kmeans seed 7).

## The two fits

| | A: as analysis 43 runs today | B: restricted to the basin |
|---|---|---|
| assemblages | **55** | **43** |
| spatial clusters | 5 | 5 |
| BN F_ST median | **0.0754** | **0.0659** |
| BN 95% HDI | [0.0392, 0.1305] | [0.0339, 0.1139] |
| Gini-Simpson median | 0.0388 | 0.0277 |
| Gini-Simpson 95% HDI | [0.0360, 0.0414] | [0.0254, 0.0300] |
| plug-in F_ST | 0.0391 | 0.0280 |
| max R-hat | 1.0001 | 1.0010 |
| min ESS | 5939 | 4375 |
| divergences | 0 | 0 |

## What is in A but not in the basin

12 assemblages: `40LA007`, `40TP026`, `Bishop`, `Fullen`, `Graves_Lake`, `Hatchie`, `Jeter`, `Jones_Bayou`, `Porter`, `Rast`, `Richardsons_Landing`, `Wilder`.

Silhouette scores on the basin coordinates, k = 2..6: k=2: 0.4555, k=3: 0.4173, k=4: 0.5022, k=5: 0.5665, k=6: 0.5138. k=5 wins clearly, which is the "three spatial clusters" the manuscript describes. The k=5 in column A is the silhouette optimum for the wider curated set, not for the basin.

## Reading

`prepare_inputs()` returns the whole curated set with coordinates. `43_bayesian_fst.basin_group_counts()` groups it by `inp.cluster_of` without restricting to `data/processed/basin_members_curated.txt`, so despite its name and despite the report headed "Observed St. Francis basin", column A is a regional quantity. The excluded 26 include `Walls`, `Wall`, `Chuccalissa` and `Parchman`, which are the Mississippi-proximal sites that `analyses/16_basin_membership.py` names in its own docstring as the reason the earlier latitude cut was replaced by the hydrological rule.

The Gini-Simpson readout, the quantity that matches the manuscript's estimator, moves from 0.0388 to 0.0277, a 29 percent reduction. The BN parameter moves less but its interval widens as it should on 29 assemblages rather than 55. Both fits are numerically healthy, so this is a scope defect and not a sampling one.

## Scope of the consequence

- **No manuscript number changes.** Figures 4 and 5 reach the data through `make_figures._load_curated()`, which does apply the basin membership. Only analyses 43 and 44 use the unrestricted path.
- `output/bayesian_fst.md` and the `STATUS.md` entry quoting "basin BN F = 0.067" and "Gini-Simpson readout 0.0327" attribute column A's numbers to the basin. Those are column B's questions with column A's answers.
- `analyses/44_bayesian_fst_validation.py` calibrates coverage, SBC and small-sample behavior at column A's geometry (5 clusters over 55) while describing it as the basin fit. The validation is not wrong, but it validates a design the paper does not use.

## Disposition

Restrict `basin_group_counts()` to the canonical membership and re-select k on the basin coordinates, then regenerate 43 and 44. Add a test that asserts the fitted set equals `data/processed/basin_members_curated.txt`, since a name alone did not prevent this. Scheduled as the first task of item 1.

## Why this matters beyond the numbers

This is the same defect class the 2026-06-10 alignment run found and fixed, where five figures were computed on the whole curated set while the text reported the basin. It recurred in code written a month later. The Verification Regime preamble in `CLAUDE.md` cites that episode as this project's own precedent for gates being necessary and never sufficient; the recurrence is the evidence that the lesson needed the test it is now getting.
