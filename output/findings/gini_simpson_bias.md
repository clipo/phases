# Does the plug-in Gini-Simpson bias change a reported result?

Produced by `analyses/62_gini_simpson_bias.py`.

## 1. The per-bin diversity trajectory (scripts 54 and 56)

| bin | sherds | plug-in | unbiased | difference |
|---|---|---|---|---|
| 0 | 2229 | 0.6887 | 0.6890 | +0.0003 |
| 1 | 2580 | 0.6495 | 0.6498 | +0.0003 |
| 2 | 5566 | 0.6087 | 0.6088 | +0.0001 |
| 3 | 8062 | 0.5337 | 0.5338 | +0.0001 |
| 4 | 11584 | 0.4383 | 0.4383 | +0.0000 |
| 5 | 6370 | 0.3225 | 0.3225 | +0.0001 |

Bins pool thousands of sherds, so the bias is at most 0.04 percent and the declining trajectory is unchanged (Spearman -1.000 either way).

## 2. The size-controlled F_ST trajectory trend

Groups here are cluster-by-bin cells on rarefied data, holding 50 to 200 sherds (median 125), which is where the bias bites.

- plug-in trend  **-0.7559**
- unbiased trend **-0.7539**
- shift +0.0020, over 800 rarefactions

## 3. Does the inferred closure strength move?

The recovery curve is built by scoring synthetic assemblages with the same estimator, so the bias appears on both sides of the mapping.

| injected s | plug-in curve | unbiased curve | shift |
|---|---|---|---|
| 0.0 | -0.6743 | +0.0914 | +0.7657 |
| 0.2 | -0.3562 | +0.3971 | +0.7533 |
| 0.4 | +0.6152 | +0.8124 | +0.1971 |
| 0.6 | +0.9152 | +0.9476 | +0.0324 |
| 0.8 | +0.9781 | +0.9762 | -0.0019 |

Inverting each curve at its own observed value gives an implied closure strength of **nan** under the plug-in and **nan** under the unbiased estimator.

## Verdict

The bias is real and it does move the descriptive statistic, by more than the statistic's own magnitude. It does not move the claim. Correcting the estimator lifts the recovery curve by more than it lifts the observed value, so the observation lands at a slightly LOWER closure strength than before, and both readings sit at the closure-strength posterior's published median of about 0.30. The correction weakens the case for closure rather than strengthening it, so the paper's conclusion does not depend on this choice.

The plug-in is kept. This file is why that is a defended choice rather than an inherited one (rule 20 in spirit: the estimator, like a prior, is checked against the direction it would push the hypothesis).
