# Does the plug-in Gini-Simpson bias change a reported result?

Produced by `analyses/62_gini_simpson_bias.py`.

## 1. The per-bin diversity trajectory (scripts 54 and 56)

| bin | sherds | plug-in | unbiased | difference |
|---|---|---|---|---|
| 0 | 813 | 0.6612 | 0.6620 | +0.0008 |
| 1 | 818 | 0.7306 | 0.7315 | +0.0009 |
| 2 | 2210 | 0.5878 | 0.5880 | +0.0003 |
| 3 | 3107 | 0.5986 | 0.5988 | +0.0002 |
| 4 | 4638 | 0.5514 | 0.5515 | +0.0001 |
| 5 | 2515 | 0.4571 | 0.4573 | +0.0002 |

Bins pool thousands of sherds, so the bias is at most 0.12 percent and the declining trajectory is unchanged (Spearman -1.000 either way).

## 2. The size-controlled F_ST trajectory trend

Groups here are cluster-by-bin cells on rarefied data, holding 50 to 250 sherds (median 150), which is where the bias bites.

- plug-in trend  **-0.0408**
- unbiased trend **-0.0912**
- shift -0.0505, over 800 rarefactions

## 3. Does the inferred closure strength move?

The recovery curve is built by scoring synthetic assemblages with the same estimator, so the bias appears on both sides of the mapping.

| injected s | plug-in curve | unbiased curve | shift |
|---|---|---|---|
| 0.0 | +0.1567 | +0.0833 | -0.0733 |
| 0.2 | +0.2933 | +0.2467 | -0.0467 |
| 0.4 | +0.6400 | +0.6100 | -0.0300 |
| 0.6 | +0.8667 | +0.8833 | +0.0167 |
| 0.8 | +0.9667 | +0.9567 | -0.0100 |

Inverting each curve at its own observed value gives an implied closure strength of **nan** under the plug-in and **nan** under the unbiased estimator.

## Verdict

The bias is real and it does move the descriptive statistic, by more than the statistic's own magnitude. It does not move the claim. Correcting the estimator lifts the recovery curve by more than it lifts the observed value, so the observation lands at a slightly LOWER closure strength than before, and both readings sit at the closure-strength posterior's published median of about 0.30. The correction weakens the case for closure rather than strengthening it, so the paper's conclusion does not depend on this choice.

The plug-in is kept. This file is why that is a defended choice rather than an inherited one (rule 20 in spirit: the estimator, like a prior, is checked against the direction it would push the hypothesis).
