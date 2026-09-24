# The neutrality envelope as a posterior predictive check

Produced by `analyses/56_neutrality_ppc.py` (full, 1500 predictive draws). Replaces the fixed-Ne sweep in `analyses/18_kandler_shennan_neutrality.py` (`docs/FREQUENTIST_INVENTORY.md` item G).

Effective population size and innovation rate are drawn from the POOLED ABC-SMC posterior (4000 draws over 5 seeds), not swept over fixed values and not taken from a single run.

## Result

Observed per-bin Gini-Simpson diversity: 0.661, 0.731, 0.588, 0.599, 0.551, 0.457.

| posterior | bins inside the 95% predictive interval | Bayesian predictive p |
|---|---|---|
| **pooled** | **6/6** | **0.885** |
| seed 101 alone | 6/6 | 0.634 |
| seed 202 alone | 6/6 | 0.885 |
| seed 303 alone | 6/6 | 0.907 |
| seed 404 alone | 6/6 | 0.915 |
| seed 505 alone | 6/6 | 0.802 |

Across single-seed posteriors the predictive p ranges 0.634 to 0.915 and the bins-inside count ranges 6 to 6 of 6. That spread is the dependency `analyses/55_abc_stability.py` measured (N median 98 to 263 across seeds) propagating into the check, and it is the reason the pooled row is the one to report.

## Reading

The neutral model reproduces the observed diversity trajectory: the predictive p sits away from both tails and the observed trajectory lies inside the predictive interval in 6 of 6 bins.

**What this does not say.** A passing predictive check says the neutral model is adequate for this trajectory, not that transmission was neutral; other processes produce the same trajectory. Kandler and Shennan's own caveat also applies, that with an ordinal axis and an uncertain effective size the predictive interval is wide, and drawing Ne from a posterior whose median moves by a factor of nearly three across seeds widens it further rather than narrowing it. This is concordant evidence, and weak, which is what the supplement already says.

## What replaces what

| before | after |
|---|---|
| sweep over fixed Ne | Ne drawn from the ABC posterior |
| 95% band read as an acceptance region | 95% posterior predictive interval |
| "the test cannot reject neutral transmission" | a Bayesian predictive p with its across-seed range |
