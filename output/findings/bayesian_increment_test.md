# The frequency increment test, made hierarchical

Produced by `analyses/52_bayesian_fit_increments.py` (full). 10 testable decorated classes over 6 seriation bins.

Replaces nine independent one-sample t-tests and a rejection count with one hierarchical model (`docs/FREQUENTIST_INVENTORY.md` item F).

## Population-level departure from drift

- mu_pop posterior median **+0.009**, 95% credible interval [-0.041, +0.083].
- P(mu_pop > 0) = 0.640.

Under neutral drift the rescaled increments have mean zero, so this interval is the direct statement about whether the assemblage as a whole departs from drift.

## Per-class posteriors, with the frequentist test alongside

| class | n increments | mean Y | t-test p | **posterior mu** | **95% CI** | P(mu > 0) |
|---|---|---|---|---|---|---|
| Parkin_Punctated | 5 | -0.12 | 0.113 | **-0.014** | [-0.158, +0.077] | 0.354 |
| Barton/Kent/MPI | 5 | +0.00 | 0.926 | **+0.005** | [-0.072, +0.096] | 0.568 |
| Painted | 5 | +0.21 | 0.002 | **+0.049** | [-0.036, +0.230] | 0.791 |
| Fortune_Noded | 4 | -0.02 | 0.302 | **-0.011** | [-0.063, +0.049] | 0.308 |
| Ranch_Incised | 5 | +0.04 | 0.652 | **+0.010** | [-0.080, +0.129] | 0.623 |
| Walls_Engraved | 4 | +0.08 | 0.592 | **+0.011** | [-0.094, +0.163] | 0.619 |
| Wallace_Incised | 3 | -0.01 | 0.795 | **+0.000** | [-0.078, +0.090] | 0.507 |
| Rhodes_Incised | 4 | +0.02 | 0.742 | **+0.008** | [-0.069, +0.104] | 0.610 |
| Vernon_Paul_Applique | 5 | +0.00 | 0.824 | **+0.004** | [-0.041, +0.053] | 0.592 |
| Hull_Engraved | 4 | +0.03 | 0.534 | **+0.012** | [-0.057, +0.111] | 0.663 |

## Diagnostics (rule 16)

- empirical fit: R-hat 1.0030, min ESS 1140, divergences 3/8000 (0.04%), E-BFMI 0.708
- recovery fit: R-hat 1.0026, min ESS 1806, divergences 1/8000, E-BFMI 0.790

## Rule 20(c): recovery against a prior that leans toward neutrality

`mu_pop ~ Normal(0, 1)` centres on exactly the conclusion the paper draws, so the prior is sympathetic and has to be shown not to be doing the work. One class was simulated with a true mean of +1.5, a value the prior disfavours and partial pooling actively shrinks; the other classes were simulated neutral.

Recovered: **+0.292** with 95% interval [-0.345, +2.414]. Coverage of the true +1.5: yes. Shrinkage of the median: **5.1x**.

**Coverage holds only because the interval is very wide. The median of a class simulated at +1.5 is dragged to +0.292, a shrinkage of 5.1x, and the interval spans zero. This design has little power to resolve a per-class departure: a genuinely drifting class would be reported as unresolved. The population-level mu_pop remains reportable, but per-class posteriors must not be read as evidence that individual classes are neutral.**

**A second caution about the empirical interval.** mu_pop on the real data is [-0.032, +0.052], far tighter than the recovery interval. That is not because the real data are more informative. With every class near zero, tau collapses toward zero and the classes are pooled almost completely, which makes mu_pop precise; in the recovery data one genuinely departing class inflates tau and widens everything. The tight empirical interval is therefore partly a consequence of the answer, and should be quoted with that stated.

## Why this also fixes a multiplicity problem

"Seven of nine testable classes are consistent with unbiased drift" counts rejections across nine uncorrected tests, where about 0.45 false rejections are expected at alpha = 0.05 even if every class is neutral. Partial pooling replaces the count with per-class posteriors that are already shrunk toward the group, so no correction is applied because none is needed.
