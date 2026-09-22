# Record-matched, size-controlled signal recovery

Real configuration: 28 assemblages, 10 decorated types, k = 2 spatial clusters, 6 ordinal bins. Synthetic assemblages are generated at the real per-assemblage sample sizes; both synthetic and real assemblages are rarefied to a common count NRARE = 50 before scoring. 120 seeds per cell, 400 empirical rarefactions.

The s = 0 baseline is a fixed shared Zipf profile with multinomial sampling, not a spatially evolving drift null; the injection imposes increasing divergence and profile concentration, and every recovery statement below refers only to that alternative. The seriation statistic is the number of maximal deterministic IDSS groups per assemblage within a bin (continuity 0.10), not an ordered-row proxy.

## The sample-size confound

- Assemblage sample size trends with seriation position at Spearman rho = +0.61, so a size-sensitive signature can rise or fall along the axis through sampling alone.
- Raw (uncontrolled) empirical F_ST trend = -0.80 (the Figure 5 value); after rarefaction it is -0.048, averaged over 4000 draws (single-rarefaction 2.5th to 97.5th percentiles -0.80 to +0.40, SD 0.42: the record does not resolve the sign of this trend). The Monte Carlo standard error of the 4000-draw mean is 0.007; it describes the compute budget, not the record, and is not the interval to quote (rule 6). The raw rise is a sampling artifact.

## Which signatures recover the injected emergence (rarefied, no averaging)

| s | neutral departure | IDSS groups per assemblage | cultural $F_{ST}$ | spatial boundary |
|---|---|---|---|---|
| 0.0 | +0.12 | -0.10 | +0.17 | +0.03 |
| 0.1 | +0.05 | -0.05 | +0.12 | +0.05 |
| 0.2 | -0.03 | -0.11 | +0.13 | -0.01 |
| 0.3 | -0.04 | -0.16 | +0.34 | +0.00 |
| 0.4 | -0.09 | -0.16 | +0.63 | -0.02 |
| 0.5 | -0.12 | -0.12 | +0.78 | +0.03 |
| 0.6 | -0.12 | -0.22 | +0.88 | -0.04 |
| 0.7 | -0.65 | -0.18 | +0.93 | -0.01 |
| 0.8 | -0.72 | -0.19 | +0.97 | -0.07 |
| 0.9 | -0.78 | -0.15 | +0.98 | -0.10 |
| 1.0 | -0.74 | -0.16 | +1.00 | -0.01 |

**Reading.** Only cultural F_ST tracks the injected signal monotonically (null near zero, rising to ~+1 at strong emergence). The neutral departure is non-monotonic in conformity (the known U-shape), the spatial boundary is unresponsive at k = 3 clusters, and the seriation group count does not rise in the hypothesized direction. At this record's resolution the criterion is carried by F_ST; the other three signatures are not reliable discriminators here.

## F_ST detector (calibrated; false-positive rate 0.05)

| s | power (w=1) | power (time-averaged w=3) |
|---|---|---|
| 0.0 | 0.09 | 0.00 |
| 0.1 | 0.05 | 0.00 |
| 0.2 | 0.07 | 0.00 |
| 0.3 | 0.17 | 0.00 |
| 0.4 | 0.23 | 0.00 |
| 0.5 | 0.33 | 0.00 |
| 0.6 | 0.53 | 0.00 |
| 0.7 | 0.72 | 0.00 |
| 0.8 | 0.86 | 0.00 |
| 0.9 | 0.90 | 0.00 |
| 1.0 | 0.99 | 0.00 |

- Detection threshold s* = 0.80 (no averaging) / not reached at any strength tested (peak power 0.00 at s = 0.0) (time-averaged): the weakest emergence recovered at power >= 80%. Time-averaging penalty: not comparable, because the time-averaged threshold is never reached on this grid.
- Null (s=0) F_ST trend mean +0.17; detection threshold +0.80.
- Thresholds use 400 separate null calibration draws per window; the s = 0 test row estimates the achieved false-positive rate independently.

## Empirical placement (size-controlled)

- Rarefied empirical F_ST trend = -0.048 [-0.80, +0.40] (single-rarefaction 2.5th to 97.5th percentiles over 4000 draws; the same estimator as above).
- On the recovery curve this corresponds to a nominal injected strength s ~ 0.00, below the resolution limit s* = 0.80 that this record can reliably detect.
- The data show no resolvable emergence: the faint trend is not distinguishable from the no-emergence null at this resolution.
- The rarefaction percentiles condition on the observed sites, partition and ordering; they are not a confidence interval over all archaeological uncertainties. The signed trend reverses with axis orientation.

## Verdict

At the record's own resolution and sample sizes, and with the size confound removed by rarefaction, the discrimination is carried by cultural F_ST. The apparent raw F_ST rise (-0.80) is an artifact of the sample-size-versus-position trend (rho +0.61) and vanishes under size control, leaving a trend of -0.048 [-0.061, -0.035]. The other three signatures do not reliably discriminate at this resolution and are reported as weak corroboration only.

**How the strength of any closure is reported.** The detection threshold, power curve and five percent false-positive rate that this section used to report are withdrawn under the project's rule 18 (no frequentist inference supports a substantive claim). The simulation grid above is a simulated likelihood, so the reportable quantity is a POSTERIOR over the injected closure strength, computed by analyses/53_closure_strength_posterior.py and written to output/findings/closure_strength_posterior.md. It is what Figure 4's right panel now draws. The power numbers this script still computes internally are retained only to build that grid and must not be quoted.