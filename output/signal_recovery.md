# Record-matched, size-controlled signal recovery

Real configuration: 43 assemblages, 10 decorated types, k = 5 spatial clusters, 6 ordinal bins. Synthetic assemblages are generated at the real per-assemblage sample sizes; both synthetic and real assemblages are rarefied to a common count NRARE = 50 before scoring. 120 seeds per cell, 400 empirical rarefactions.

The s = 0 baseline is a fixed shared Zipf profile with multinomial sampling, not a spatially evolving drift null; the injection imposes increasing divergence and profile concentration, and every recovery statement below refers only to that alternative. The seriation statistic is the number of maximal deterministic IDSS groups per assemblage within a bin (continuity 0.10), not an ordered-row proxy.

## The sample-size confound

- Assemblage sample size trends with seriation position at Spearman rho = +0.33, so a size-sensitive signature can rise or fall along the axis through sampling alone.
- Raw (uncontrolled) empirical F_ST trend = -0.71 (the Figure 5 value); after rarefaction it is -0.756 (Monte Carlo 95% interval [-0.758, -0.755] over 4000 draws; a single rarefaction has SD 0.05, so this statistic is not stable in sign at a few hundred draws). The raw rise is a sampling artifact.

## Which signatures recover the injected emergence (rarefied, no averaging)

| s | neutral departure | IDSS groups per assemblage | cultural $F_{ST}$ | spatial boundary |
|---|---|---|---|---|
| 0.0 | +0.59 | -0.08 | -0.71 | +0.06 |
| 0.1 | +0.53 | -0.17 | -0.58 | -0.02 |
| 0.2 | +0.49 | -0.26 | -0.38 | +0.01 |
| 0.3 | +0.48 | -0.51 | +0.18 | +0.02 |
| 0.4 | +0.38 | -0.68 | +0.64 | +0.02 |
| 0.5 | +0.40 | -0.68 | +0.82 | +0.06 |
| 0.6 | +0.03 | -0.70 | +0.91 | -0.00 |
| 0.7 | -0.33 | -0.74 | +0.95 | -0.00 |
| 0.8 | -0.56 | -0.67 | +0.98 | -0.05 |
| 0.9 | -0.71 | -0.62 | +0.99 | -0.02 |
| 1.0 | -0.70 | -0.62 | +0.99 | +0.00 |

**Reading.** Only cultural F_ST tracks the injected signal monotonically (null near zero, rising to ~+1 at strong emergence). The neutral departure is non-monotonic in conformity (the known U-shape), the spatial boundary is unresponsive at k = 3 clusters, and the seriation group count does not rise in the hypothesized direction. At this record's resolution the criterion is carried by F_ST; the other three signatures are not reliable discriminators here.

## F_ST detector (calibrated; false-positive rate 0.05)

| s | power (w=1) | power (time-averaged w=3) |
|---|---|---|
| 0.0 | 0.07 | 0.03 |
| 0.1 | 0.12 | 0.10 |
| 0.2 | 0.35 | 0.22 |
| 0.3 | 0.82 | 0.45 |
| 0.4 | 1.00 | 0.87 |
| 0.5 | 1.00 | 1.00 |
| 0.6 | 1.00 | 1.00 |
| 0.7 | 1.00 | 1.00 |
| 0.8 | 1.00 | 1.00 |
| 0.9 | 1.00 | 1.00 |
| 1.0 | 1.00 | 1.00 |

- Detection threshold s* = 0.30 (no averaging) / 0.40 (time-averaged): the weakest emergence recovered at power >= 80%. Time-averaging penalty: +0.10.
- Null (s=0) F_ST trend mean -0.71; detection threshold -0.26.
- Thresholds use 400 separate null calibration draws per window; the s = 0 test row estimates the achieved false-positive rate independently.

## Empirical placement (size-controlled)

- Rarefied empirical F_ST trend = -0.75 [-0.83, -0.60].
- On the recovery curve this corresponds to a nominal injected strength s ~ 0.00, below the resolution limit s* = 0.30 that this record can reliably detect.
- The data show no resolvable emergence: the faint trend is not distinguishable from the no-emergence null at this resolution.
- The rarefaction percentiles condition on the observed sites, partition and ordering; they are not a confidence interval over all archaeological uncertainties. The signed trend reverses with axis orientation.

## Verdict

At the record's own resolution and sample sizes, and with the size confound removed by rarefaction, the discrimination is carried by cultural F_ST. The apparent raw F_ST rise (-0.71) is an artifact of the sample-size-versus-position trend (rho +0.33) and vanishes under size control, leaving a trend of -0.756 [-0.758, -0.755]. The other three signatures do not reliably discriminate at this resolution and are reported as weak corroboration only.

**How the strength of any closure is reported.** The detection threshold, power curve and five percent false-positive rate that this section used to report are withdrawn under the project's rule 18 (no frequentist inference supports a substantive claim). The simulation grid above is a simulated likelihood, so the reportable quantity is a POSTERIOR over the injected closure strength, computed by analyses/53_closure_strength_posterior.py and written to output/findings/closure_strength_posterior.md. It is what Figure 4's right panel now draws. The power numbers this script still computes internally are retained only to build that grid and must not be quoted.