# The frequency increment test, made hierarchical

Produced by `analyses/52_bayesian_fit_increments.py` (full). 9 testable decorated classes over 6 seriation bins.

Replaces nine independent one-sample t-tests and a rejection count with one hierarchical model (`docs/FREQUENTIST_INVENTORY.md` item F).

## Population-level departure from drift

- mu_pop posterior median **+0.001**, 95% credible interval [-0.059, +0.064].
- P(mu_pop > 0) = 0.518.

Under neutral drift the rescaled increments have mean zero, so this interval is the direct statement about whether the assemblage as a whole departs from drift.

## Per-class posteriors, with the frequentist test alongside

| class | n increments | mean Y | t-test p | **posterior mu** | **95% CI** | P(mu > 0) |
|---|---|---|---|---|---|---|
| Parkin_Punctated | 5 | -0.15 | 0.000 | **-0.151** | [-0.167, -0.001] | 0.023 |
| Barton/Kent/MPI | 5 | +0.07 | 0.052 | **+0.059** | [-0.028, +0.127] | 0.918 |
| Painted | 5 | +0.15 | 0.015 | **+0.101** | [-0.042, +0.202] | 0.912 |
| Fortune_Noded | 5 | -0.00 | 0.885 | **-0.002** | [-0.045, +0.043] | 0.455 |
| Ranch_Incised | 5 | -0.01 | 0.225 | **-0.013** | [-0.045, +0.023] | 0.171 |
| Walls_Engraved | 5 | +0.02 | 0.450 | **+0.015** | [-0.051, +0.076] | 0.729 |
| Rhodes_Incised | 5 | +0.01 | 0.591 | **+0.009** | [-0.044, +0.058] | 0.673 |
| Vernon_Paul_Applique | 4 | +0.00 | 0.902 | **+0.002** | [-0.060, +0.063] | 0.539 |
| Hull_Engraved | 4 | +0.01 | 0.753 | **+0.006** | [-0.058, +0.067] | 0.604 |

## Diagnostics (rule 16)

- empirical fit: R-hat 1.0183, min ESS 487, divergences 14/8000 (0.17%), E-BFMI 0.433
- recovery fit: R-hat 1.0044, min ESS 1485, divergences 2/8000, E-BFMI 0.738

## Rule 20(c): recovery against a prior that leans toward neutrality

`mu_pop ~ Normal(0, 1)` centres on exactly the conclusion the paper draws, so the prior is sympathetic and has to be shown not to be doing the work. One class was simulated with a true mean of +1.5, a value the prior disfavours and partial pooling actively shrinks; the other classes were simulated neutral.

Recovered: **+0.381** with 95% interval [-0.382, +2.407]. Coverage of the true +1.5: yes. Shrinkage of the median: **3.9x**.

**Coverage holds only because the interval is very wide. The median of a class simulated at +1.5 is dragged to +0.381, a shrinkage of 3.9x, and the interval spans zero. This design has little power to resolve a per-class departure: a genuinely drifting class would be reported as unresolved. The population-level mu_pop remains reportable, but per-class posteriors must not be read as evidence that individual classes are neutral.**

**A second caution about the empirical interval.** mu_pop on the real data is [-0.032, +0.052], far tighter than the recovery interval. That is not because the real data are more informative. With every class near zero, tau collapses toward zero and the classes are pooled almost completely, which makes mu_pop precise; in the recovery data one genuinely departing class inflates tau and widens everything. The tight empirical interval is therefore partly a consequence of the answer, and should be quoted with that stated.

## Why this also fixes a multiplicity problem

"Seven of nine testable classes are consistent with unbiased drift" counts rejections across nine uncorrected tests, where about 0.45 false rejections are expected at alpha = 0.05 even if every class is neutral. Partial pooling replaces the count with per-class posteriors that are already shrunk toward the group, so no correction is applied because none is needed.
