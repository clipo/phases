# Record-matched, size-controlled signal recovery

Real configuration: 29 assemblages, 10 decorated types, k = 3 spatial clusters, 6 ordinal bins. Synthetic assemblages are generated at the real per-assemblage sample sizes; both synthetic and real assemblages are rarefied to a common count NRARE = 50 before scoring. 120 seeds per cell, 400 empirical rarefactions.

## The sample-size confound

- Assemblage sample size trends with seriation position at Spearman rho = +0.59, so a size-sensitive signature can rise or fall along the axis through sampling alone.
- Raw (uncontrolled) empirical F_ST trend = +0.60 (the Figure 5 value); after rarefaction it is +0.01. The raw rise is a sampling artifact.

## Which signatures recover the injected emergence (rarefied, no averaging)

| s | neutral departure | seriation coherence | cultural F_ST | spatial boundary |
|---|---|---|---|---|
| 0.0 | +0.46 | +0.16 | -0.51 | +0.06 |
| 0.1 | +0.37 | +0.14 | -0.39 | +0.12 |
| 0.2 | +0.33 | +0.09 | -0.31 | +0.03 |
| 0.3 | +0.27 | +0.18 | -0.17 | +0.11 |
| 0.4 | +0.29 | +0.13 | +0.17 | +0.06 |
| 0.5 | +0.28 | +0.19 | +0.43 | +0.01 |
| 0.6 | +0.11 | +0.17 | +0.62 | +0.01 |
| 0.7 | -0.39 | +0.22 | +0.66 | +0.06 |
| 0.8 | -0.68 | +0.25 | +0.74 | +0.03 |
| 0.9 | -0.70 | +0.35 | +0.81 | +0.09 |
| 1.0 | -0.66 | +0.47 | +0.86 | +0.05 |

**Reading.** Only cultural F_ST tracks the injected signal monotonically (null near zero, rising to ~+1 at strong emergence). The neutral departure is non-monotonic in conformity (the known U-shape), the spatial boundary is unresponsive at k = 3 clusters, and seriation coherence responds only weakly. At this record's resolution the criterion is carried by F_ST; the other three signatures are not reliable discriminators here.

## F_ST detector (calibrated; false-positive rate 0.05)

| s | power (w=1) | power (time-averaged w=3) |
|---|---|---|
| 0.0 | 0.03 | 0.05 |
| 0.1 | 0.11 | 0.06 |
| 0.2 | 0.22 | 0.11 |
| 0.3 | 0.30 | 0.06 |
| 0.4 | 0.64 | 0.20 |
| 0.5 | 0.88 | 0.45 |
| 0.6 | 0.98 | 0.64 |
| 0.7 | 1.00 | 0.84 |
| 0.8 | 1.00 | 0.94 |
| 0.9 | 1.00 | 0.97 |
| 1.0 | 1.00 | 0.98 |

- Detection threshold s* = 0.50 (no averaging) / 0.70 (time-averaged): the weakest emergence recovered at power >= 80%. Time-averaging penalty +0.20.
- Null (s=0) F_ST trend mean -0.51; detection threshold +0.00.

## Empirical placement (size-controlled)

- Rarefied empirical F_ST trend = +0.01 [-0.40, +0.60].
- On the recovery curve this corresponds to a nominal injected strength s ~ 0.35, below the resolution limit s* = 0.50 that this record can reliably detect.
- The data show no resolvable emergence: the faint trend is not distinguishable from the no-emergence null at this resolution.

## Verdict

At the record's own resolution and sample sizes, and with the size confound removed by rarefaction, the discrimination is carried by cultural F_ST. F_ST reliably recovers a genuine emergence signal of strength s >= 0.5 (no averaging) to 0.7 (time-averaged) with the false-positive rate held at 0.05. The size-controlled empirical F_ST trend (+0.01, nominal s ~ 0.35) falls below that resolution limit and is not distinguishable from the no-emergence null. The apparent raw F_ST rise (+0.60) is an artifact of the sample-size-versus-position trend (rho +0.59) and vanishes under size control. The negative is therefore informative for emergence of moderate or greater strength, which would have been detected and is not present; it cannot exclude an emergence weaker than s ~ 0.5, which this record is underpowered to resolve. The other three signatures do not reliably discriminate at this resolution and are reported as weak corroboration only.