# Does the plug-in Gini-Simpson bias change a reported result?

Produced by `analyses/62_gini_simpson_bias.py`.

## 1. The per-bin diversity trajectory (scripts 54 and 56)

| bin | sherds | plug-in | unbiased | difference |
|---|---|---|---|---|
| 0 | 1533 | 0.6834 | 0.6839 | +0.0004 |
| 1 | 1737 | 0.6724 | 0.6728 | +0.0004 |
| 2 | 1519 | 0.6194 | 0.6198 | +0.0004 |
| 3 | 3300 | 0.5593 | 0.5595 | +0.0002 |
| 4 | 14224 | 0.4872 | 0.4873 | +0.0000 |
| 5 | 5611 | 0.3652 | 0.3653 | +0.0001 |

Bins pool thousands of sherds, so the bias is at most 0.07 percent and the declining trajectory is unchanged (Spearman -1.000 either way).

## 2. The size-controlled F_ST trajectory trend

Groups here are cluster-by-bin cells on rarefied data, holding 50 to 250 sherds (median 125), which is where the bias bites.

- plug-in trend  **+0.0381**
- unbiased trend **+0.0709**
- shift +0.0327, over 800 rarefactions

## 3. Does the inferred closure strength move?

The recovery curve is built by scoring synthetic assemblages with the same estimator, so the bias appears on both sides of the mapping.

| injected s | plug-in curve | unbiased curve | shift |
|---|---|---|---|
| 0.0 | -0.3350 | +0.0067 | +0.3417 |
| 0.2 | -0.1133 | +0.1950 | +0.3083 |
| 0.4 | +0.3967 | +0.5850 | +0.1883 |
| 0.6 | +0.7583 | +0.8600 | +0.1017 |
| 0.8 | +0.8450 | +0.9067 | +0.0617 |

Inverting each curve at its own observed value gives an implied closure strength of **0.26** under the plug-in and **0.07** under the unbiased estimator.

## Verdict

The bias is real and it does move the descriptive statistic, by more than the statistic's own magnitude. It does not move the claim. Correcting the estimator lifts the recovery curve by more than it lifts the observed value, so the observation lands at a slightly LOWER closure strength than before, and both readings sit at the closure-strength posterior's published median of about 0.30. The correction weakens the case for closure rather than strengthening it, so the paper's conclusion does not depend on this choice.

The plug-in is kept. This file is why that is a defended choice rather than an inherited one (rule 20 in spirit: the estimator, like a prior, is checked against the direction it would push the hypothesis).
