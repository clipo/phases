# A posterior on closure strength

Produced by `analyses/53_closure_strength_posterior.py` (full: 120 seeds per grid point, 400 empirical rarefactions).

Replaces the calibrated detector, its five percent false-positive rate, its power curve and its detection threshold with one quantity (`docs/FREQUENTIST_INVENTORY.md` item A). No test is performed.

Observed statistic: the rarefied F_ST trajectory rank correlation, T = **-0.7549**.

## Posterior over injected closure strength s

| time-averaging | median | 95% credible interval | P(s >= 0.3) | P(s >= 0.5) |
|---|---|---|---|---|
| none (w = 1) | **0.069** | [0.000, 0.252] | 0.010 | 0.000 |
| record's (w = 3) | **0.096** | [0.000, 0.314] | 0.057 | 0.000 |

## Rule 20(c): recovery of a closure the data do not contain

Data simulated at a true s = 0.7, well above anything the basin shows, returns a posterior median of **0.843** with interval [0.477, 0.993], which covers the truth.

A strong closure is recovered rather than dragged toward the middle, so a low posterior on the real data is a statement about the data.

## Reading

The old apparatus said closure of strength 0.5 or more would be detected with power above 0.8, and read the empirical value as a nominal strength near 0.35. The posterior says the same thing without a null: strengths at or above 0.5 carry posterior mass 0.000 without time-averaging and 0.000 with it, so moderate and strong closure are excluded by the data rather than by a threshold, while weak closure remains credible. That is the bounded non-detection the paper argues for, stated as a parameter.

## Declared limitation

p(T | s) is approximated by a normal at each grid point, with mean and standard deviation taken across seeds and interpolated between them. T is a bounded rank correlation, so the approximation is reasonable in the interior and degrades near the ends. This is a simulation-based posterior, not an exact one, and is reported as such.
