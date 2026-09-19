# The frequency increment test, made hierarchical

Produced by `analyses/52_bayesian_fit_increments.py` (full). 9 testable decorated classes over 6 seriation bins.

Replaces nine independent one-sample t-tests and a rejection count with one hierarchical model (`docs/FREQUENTIST_INVENTORY.md` item F).

## Population-level departure from drift

- mu_pop posterior median **+0.008**, 95% credible interval [-0.035, +0.059].
- P(mu_pop > 0) = 0.692.

Under neutral drift the rescaled increments have mean zero, so this interval is the direct statement about whether the assemblage as a whole departs from drift.

## Per-class posteriors, with the frequentist test alongside

| class | n increments | mean Y | t-test p | **posterior mu** | **95% CI** | P(mu > 0) |
|---|---|---|---|---|---|---|
| Parkin_Punctated | 5 | -0.13 | 0.021 | **-0.011** | [-0.148, +0.046] | 0.360 |
| Barton/Kent/MPI | 5 | +0.03 | 0.552 | **+0.010** | [-0.058, +0.093] | 0.677 |
| Painted | 5 | +0.18 | 0.005 | **+0.024** | [-0.028, +0.197] | 0.803 |
| Fortune_Noded | 5 | -0.01 | 0.369 | **-0.005** | [-0.040, +0.029] | 0.364 |
| Ranch_Incised | 5 | +0.03 | 0.606 | **+0.010** | [-0.061, +0.097] | 0.669 |
| Walls_Engraved | 4 | +0.03 | 0.187 | **+0.017** | [-0.024, +0.068] | 0.813 |
| Rhodes_Incised | 4 | +0.00 | 0.902 | **+0.004** | [-0.038, +0.047] | 0.613 |
| Vernon_Paul_Applique | 3 | +0.00 | 0.960 | **+0.005** | [-0.056, +0.063] | 0.593 |
| Hull_Engraved | 4 | +0.01 | 0.571 | **+0.010** | [-0.034, +0.067] | 0.703 |

## Diagnostics (rule 16)

- empirical fit: R-hat 1.0058, min ESS 893, divergences 34/8000 (0.42%), E-BFMI 0.671
- recovery fit: R-hat 1.0024, min ESS 1524, divergences 9/8000, E-BFMI 0.754

## Rule 20(c): recovery against a prior that leans toward neutrality

`mu_pop ~ Normal(0, 1)` centres on exactly the conclusion the paper draws, so the prior is sympathetic and has to be shown not to be doing the work. One class was simulated with a true mean of +1.5, a value the prior disfavours and partial pooling actively shrinks; the other classes were simulated neutral.

Recovered: **+0.444** with 95% interval [-0.397, +2.054]. Coverage of the true +1.5: yes. Shrinkage of the median: **3.4x**.

**Coverage holds only because the interval is very wide. The median of a class simulated at +1.5 is dragged to +0.444, a shrinkage of 3.4x, and the interval spans zero. This design has little power to resolve a per-class departure: a genuinely drifting class would be reported as unresolved. The population-level mu_pop remains reportable, but per-class posteriors must not be read as evidence that individual classes are neutral.**

**A second caution about the empirical interval.** mu_pop on the real data is [-0.032, +0.052], far tighter than the recovery interval. That is not because the real data are more informative. With every class near zero, tau collapses toward zero and the classes are pooled almost completely, which makes mu_pop precise; in the recovery data one genuinely departing class inflates tau and widens everything. The tight empirical interval is therefore partly a consequence of the answer, and should be quoted with that stated.

## Why this also fixes a multiplicity problem

"Seven of nine testable classes are consistent with unbiased drift" counts rejections across nine uncorrected tests, where about 0.45 false rejections are expected at alpha = 0.05 even if every class is neutral. Partial pooling replaces the count with per-class posteriors that are already shrunk toward the group, so no correction is applied because none is needed.
