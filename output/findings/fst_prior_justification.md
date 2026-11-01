# What does the F_ST prior imply about data?

Produced by `analyses/57_fst_prior_predictive.py` (full, 4000 draws per prior). Closes F2 and supplies rule 20(a) for the Balding-Nichols prior.

Design: the St. Francis basin as analysis 43 fits it, 2 spatial clusters over 10 decorated classes at the real per-cluster sherd totals [7317, 6784]. **Observed Gini-Simpson F_ST = 0.0062.**

## Induced prior predictive distribution of F_ST

| prior | median | 95% predictive interval | P(F_ST > 0.30) |
|---|---|---|---|
| Uniform(0,1) | 0.236 | [0.004, 1.000] | **0.439** |
| Beta(1,3) | 0.087 | [0.003, 0.682] | **0.168** |
| Beta(1,10) | 0.026 | [0.001, 0.247] | **0.015** |
| Beta(2,20) | 0.032 | [0.003, 0.185] | **0.003** |

## Reading

Before seeing any data, `F ~ Uniform(0, 1)` expects a median Gini-Simpson F_ST of **0.236** at this design and puts **44 percent** of its mass above 0.30. The observed value is 0.0062. Nobody holds the belief that decorated-ceramic assemblages a few tens of kilometres apart in one drainage are that strongly differentiated, so on the scale that matters the flat prior is not uninformative; it is strongly and wrongly informative.

`Beta(1, 10)` expects a median of 0.026 with 2 percent above 0.30, which is a defensible prior belief for this setting.

## The family constraint

Any `Beta(a, b)` with `a > 1` has density exactly zero at F = 0 and so asserts a priori that panmixia is impossible. Because the whole comparison in analyses 43 and 44 is against panmixia, that is disqualifying rather than merely unattractive. Beta(2, 20) is included in the table only to show what the excluded shape does; the admissible family is `Beta(1, b)`, of which Uniform(0, 1) is the b = 1 member.

## Direction, and why the burden is asymmetric (rule 20)

Uniform(0, 1) leans toward LARGE F_ST. For the posterior that is **conservative**: the paper argues F_ST is small, and a prior pulling the other way cannot manufacture that conclusion, which is why the reported posteriors are not in doubt. For a Bayes factor the lean reverses and becomes **sympathetic**: spreading prior mass over values the data never support penalises the structured model against panmixia, which is the direction of the paper's own reading. A prior can lean two ways at once, and rule 20 requires both be named.

## Disposition

The prior predictive disqualifies Uniform(0, 1) as a description of prior belief, and `../mataa` reached the same verdict independently on three of its own designs. Because the lean is conservative for the posterior, no reported posterior is overturned by this; what changes is that the prior can no longer be defended as uninformative, and any Bayes factor computed under it is reading the prior as much as the data. Adopting Beta(1, 10) as primary with Uniform reported alongside is the change this finding supports.
