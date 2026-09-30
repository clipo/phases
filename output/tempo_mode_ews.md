# Tempo-and-mode and early-warning-signal tests

Basin curated set (n = 28), 10 decorated types, oriented CA axis.

## (1) Tempo and mode (Hunt 2006/2008; OU-vs-BM model selection)

Per-type trajectories tested: 5. Models fit on the 6-bin level vector with sampling variances; comparison by AICc.

**Decorated-type trajectories, best-model tally:**

| model | types selecting it | mean Akaike weight |
|---|---|---|
| BM | 2/5 | 0.55 |
| GRW | 0/5 | 0.03 |
| Stasis | 3/5 | 0.42 |
| OU | 0/5 | 0.00 |

**Gini-Simpson diversity trajectory:**

| model | AICc | Akaike weight |
|---|---|---|
| BM | -10.6 | 0.94 |
| GRW | -1.9 | 0.01 |
| Stasis | -4.5 | 0.05 |
| OU | 27.9 | 0.00 |

- Diversity trajectory best model: **BM** (weight 0.94), directional GRW weight 0.01.
- Directional (GRW) is selected for 0 of 5 type trajectories.

**Reading.** We lead with the Gini-Simpson diversity trajectory, which the correspondence-analysis ordering does not arrange by construction. It selects the unbiased random walk decisively (BM weight 0.94), the neutral-drift expectation, with negligible support for the directional, stasis, or OU alternatives. The per-type tally favors the directional GRW for some types, but that comparison is weak and partly tautological: the CA axis is built to maximize monotonic type-frequency turnover, so individual type trajectories trend along it almost by construction, and only a few types are abundant enough across the sequence to fit. A regime transition would register as a coherent directional shift in the aggregate trait, which it does not. The robust reading is drift, not a sustained directional transition.

## (2) Early-warning signals (Scheffer et al. 2009; Dakos et al. 2012)

Trait: Gini-Simpson diversity of the 28 assemblages ordered by CA position. Gaussian detrending (bandwidth 3), rolling window 14 of 28.

- Rolling-variance trend: Kendall tau = +0.09 (permutation p = 0.450).
- Lag-1 autocorrelation trend: Kendall tau = -0.75 (permutation p = 0.959).
- Sensitivity across window 0.4-0.6 and bandwidth 2-4: variance tau in [-0.22, +0.58], AR(1) tau in [-0.77, -0.57].

**Reading.** The indicators are not jointly rising beyond the permutation null (no early-warning signal). Critical slowing down toward a bifurcation inflates the variance and the lag-1 autocorrelation together, and a reliable signal requires both. Here the rolling variance shows no rising trend (tau +0.09, p 0.45) and is unstable in sign across the sensitivity sweep, while the autocorrelation rises (tau -0.75, p 0.96). The lone autocorrelation trend is the less diagnostic of the two and is partly expected from the CA ordering itself, which places compositionally similar assemblages adjacent and so structures short-range autocorrelation. With the two indicators not rising together, and with only 29 ordinal points, there is no robust early-warning signal of an approaching transition. The result is concordant with the static criterion and the tempo test.

## Verdict

Neither dynamic probe recovers a transition the static criterion missed. The aggregate diversity trajectory is best described by neutral drift (BM) rather than a directional shift, and the early-warning indicators do not rise together toward contact. Both are concordant with the no-convergence reading, now tested in the time domain rather than only on static states.