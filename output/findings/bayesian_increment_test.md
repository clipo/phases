# The frequency increment test, made hierarchical

Produced by `analyses/52_bayesian_fit_increments.py` (full). 10 testable decorated classes over 6 seriation bins.

Replaces nine independent one-sample t-tests and a rejection count with one hierarchical model (`docs/FREQUENTIST_INVENTORY.md` item F).

## Population-level departure from drift

- mu_pop posterior median **+0.008**, 95% credible interval [-0.040, +0.081].
- P(mu_pop > 0) = 0.638.

Under neutral drift the rescaled increments have mean zero, so this interval is the direct statement about whether the assemblage as a whole departs from drift.

## Per-class posteriors, with the frequentist test alongside

| class | n increments | mean Y | t-test p | **posterior mu** | **95% CI** | P(mu > 0) |
|---|---|---|---|---|---|---|
| Parkin_Punctated | 5 | -0.12 | 0.113 | **-0.014** | [-0.161, +0.073] | 0.349 |
| Barton/Kent/MPI | 5 | +0.00 | 0.926 | **+0.004** | [-0.074, +0.099] | 0.559 |
| Painted | 5 | +0.21 | 0.002 | **+0.045** | [-0.038, +0.228] | 0.784 |
| Fortune_Noded | 4 | -0.02 | 0.302 | **-0.011** | [-0.062, +0.048] | 0.313 |
| Ranch_Incised | 5 | +0.04 | 0.652 | **+0.010** | [-0.078, +0.136] | 0.613 |
| Walls_Engraved | 4 | +0.08 | 0.592 | **+0.010** | [-0.093, +0.166] | 0.623 |
| Wallace_Incised | 3 | -0.01 | 0.795 | **-0.000** | [-0.076, +0.094] | 0.494 |
| Rhodes_Incised | 4 | +0.02 | 0.742 | **+0.008** | [-0.067, +0.108] | 0.613 |
| Vernon_Paul_Applique | 5 | +0.00 | 0.824 | **+0.004** | [-0.038, +0.052] | 0.582 |
| Hull_Engraved | 4 | +0.03 | 0.534 | **+0.012** | [-0.064, +0.113] | 0.653 |

## Diagnostics (rule 16)

- empirical fit: R-hat 1.0039, min ESS 989, divergences 1/8000 (0.01%), E-BFMI 0.714
- recovery fit: R-hat 1.0017, min ESS 1714, divergences 1/8000, E-BFMI 0.775

## Rule 20(c): recovery against a prior that leans toward neutrality

`mu_pop ~ Normal(0, 1)` centres on exactly the conclusion the paper draws, so the prior is sympathetic and has to be shown not to be doing the work. One class was simulated with a true mean of +1.5, a value the prior disfavours and partial pooling actively shrinks; the other classes were simulated neutral.

Recovered: **+0.291** with 95% interval [-0.363, +2.430]. Coverage of the true +1.5: yes. Shrinkage of the median: **5.2x**.

**Coverage holds only because the interval is very wide. The median of a class simulated at +1.5 is dragged to +0.291, a shrinkage of 5.2x, and the interval spans zero. This design has little power to resolve a per-class departure: a genuinely drifting class would be reported as unresolved. The population-level mu_pop remains reportable, but per-class posteriors must not be read as evidence that individual classes are neutral.**

**A second caution about the empirical interval.** mu_pop on the real data is [-0.032, +0.052], far tighter than the recovery interval. That is not because the real data are more informative. With every class near zero, tau collapses toward zero and the classes are pooled almost completely, which makes mu_pop precise; in the recovery data one genuinely departing class inflates tau and widens everything. The tight empirical interval is therefore partly a consequence of the answer, and should be quoted with that stated.

## Why this also fixes a multiplicity problem

"Seven of nine testable classes are consistent with unbiased drift" counts rejections across nine uncorrected tests, where about 0.45 false rejections are expected at alpha = 0.05 even if every class is neutral. Partial pooling replaces the count with per-class posteriors that are already shrunk toward the group, so no correction is applied because none is needed.
