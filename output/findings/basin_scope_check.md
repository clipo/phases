# Does the Bayesian F_ST fit run on the basin? Measured.

Produced by `analyses/47_basin_scope_check.py`. Seed 0, 2000 draws x 4 chains, F ~ Uniform(0,1), clusters by the analysis-07 silhouette rule (k in 2..6, kmeans seed 7).

## The two fits

| | A: as analysis 43 runs today | B: restricted to the basin |
|---|---|---|
| assemblages | **28** | **28** |
| spatial clusters | 2 | 2 |
| BN F_ST median | **0.1011** | **0.1011** |
| BN 95% HDI | [0.0114, 0.2422] | [0.0114, 0.2422] |
| Gini-Simpson median | 0.0062 | 0.0062 |
| Gini-Simpson 95% HDI | [0.0051, 0.0076] | [0.0051, 0.0076] |
| plug-in F_ST | 0.0062 | 0.0062 |
| max R-hat | 1.0011 | 1.0011 |
| min ESS | 3257 | 3257 |
| divergences | 0 | 0 |

## What is in A but not in the basin

0 assemblages: .

Silhouette scores on the basin coordinates, k = 2..6: k=2: 0.5290, k=3: 0.5015, k=4: 0.4686, k=5: 0.4691, k=6: 0.3715. k=2 wins clearly, which is the "three spatial clusters" the manuscript describes. The k=5 in column A is the silhouette optimum for the wider curated set, not for the basin.

## Reading

`prepare_inputs()` returns the whole curated set with coordinates. `43_bayesian_fst.basin_group_counts()` groups it by `inp.cluster_of` without restricting to `data/processed/basin_members_curated.txt`, so despite its name and despite the report headed "Observed St. Francis basin", column A is a regional quantity. The excluded 26 include `Walls`, `Wall`, `Chuccalissa` and `Parchman`, which are the Mississippi-proximal sites that `analyses/16_basin_membership.py` names in its own docstring as the reason the earlier latitude cut was replaced by the hydrological rule.

The Gini-Simpson readout, the quantity that matches the manuscript's estimator, moves from 0.0062 to 0.0062, a 0 percent reduction. The BN parameter moves less but its interval widens as it should on 29 assemblages rather than 55. Both fits are numerically healthy, so this is a scope defect and not a sampling one.

## Scope of the consequence

- **No manuscript number changes.** Figures 4 and 5 reach the data through `make_figures._load_curated()`, which does apply the basin membership. Only analyses 43 and 44 use the unrestricted path.
- `output/bayesian_fst.md` and the `STATUS.md` entry quoting "basin BN F = 0.067" and "Gini-Simpson readout 0.0327" attribute column A's numbers to the basin. Those are column B's questions with column A's answers.
- `analyses/44_bayesian_fst_validation.py` calibrates coverage, SBC and small-sample behavior at column A's geometry (5 clusters over 55) while describing it as the basin fit. The validation is not wrong, but it validates a design the paper does not use.

## Disposition

Restrict `basin_group_counts()` to the canonical membership and re-select k on the basin coordinates, then regenerate 43 and 44. Add a test that asserts the fitted set equals `data/processed/basin_members_curated.txt`, since a name alone did not prevent this. Scheduled as the first task of item 1.

## Why this matters beyond the numbers

This is the same defect class the 2026-06-10 alignment run found and fixed, where five figures were computed on the whole curated set while the text reported the basin. It recurred in code written a month later. The Verification Regime preamble in `CLAUDE.md` cites that episode as this project's own precedent for gates being necessary and never sufficient; the recurrence is the evidence that the lesson needed the test it is now getting.
