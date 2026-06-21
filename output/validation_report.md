# Phase 4: Criterion Validation by Simulation

This report tests the paper's core claim: that CONVERGENCE of four cultural-transmission signatures discriminates genuine group-level emergence from mimics. The harness applies an identical signature pipeline to every mechanism, blind to its name; a mechanism is flagged 'convergent' iff all four signatures show a positive ordinal trend whose standardized slope exceeds 0.1 (standard deviations of the signature per ordinal step).

The chronology is ORDINAL (relative phases), so 'trajectory' means a monotonic trend across the ordinal slice index, not a calendar rate.

## 1. Four-signature panels per mechanism

Panels shown for seed 42. Columns: neutral_departure, seriability, fst, spatial_boundary.

### group_emergence

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.293        -55.0  0.002             0.237
1              0.586        -44.0  0.624           132.876
2              0.719        -21.0  0.804           137.213
3              0.773        -22.0  0.911            77.805
4              0.872        -12.0  0.957           139.128
5              0.865         -9.0  0.973           140.745
6              0.929         -8.0  0.989           161.156
7              0.949         -7.0  0.996           161.268
```

### aggregated_conformity

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.316        -54.0  0.002            -0.072
1              0.641        -37.0  0.004             1.389
2              0.731        -21.0  0.003            -0.767
3              0.767        -14.0  0.003            -0.589
4              0.868        -12.0  0.003            -0.067
5              0.889        -10.0  0.007            -0.872
6              0.942        -10.0  0.003            -0.189
7              0.959         -8.0  0.002            -0.006
```

### patchiness

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.316        -57.0  0.317           102.612
1              0.314        -61.0  0.324           101.862
2              0.258        -66.0  0.299            98.866
3              0.333        -56.0  0.315            99.987
4              0.312        -49.0  0.325           100.425
5              0.279        -57.0  0.318           105.132
6              0.282        -58.0  0.308           104.047
7              0.314        -56.0  0.320            97.759
```

### drift_space

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.296        -44.0  0.232            32.667
1              0.280        -47.0  0.231            36.539
2              0.330        -50.0  0.253            33.406
3              0.295        -53.0  0.243            39.019
4              0.244        -46.0  0.221            33.067
5              0.304        -44.0  0.254            43.250
6              0.280        -43.0  0.236            33.075
7              0.293        -47.0  0.240            31.622
```

## 2. Discrimination verdict

Threshold on standardized ordinal slope: 0.1.

| mechanism | neutral_departure | seriability | fst | spatial_boundary | CONVERGENT |
|---|---|---|---|---|---|
| group_emergence | +0.300 | +0.340 | +0.343 | +0.292 | **True** |
| aggregated_conformity | +0.289 | +0.298 | +0.000 | -0.002 | **False** |
| patchiness | -0.006 | +0.034 | +0.000 | -0.001 | **False** |
| drift_space | -0.007 | +0.015 | +0.002 | -0.000 | **False** |

Each cell is the standardized ordinal slope of that signature. A mechanism is convergent only when all four are above threshold.

### Cross-seed robustness

Across 500 seeds: genuine emergence flagged convergent in 493/500 runs (sensitivity); mimics flagged convergent in 0 runs (false positives).

No mimic was ever flagged convergent. Specificity is the load-bearing property of the criterion and it holds at 100%.

Sensitivity is below 100%: on an occasional single realization the spatial-boundary signature carries enough sampling noise that its ordinal trend dips below threshold even though the other three rise. This is a power limit on a single noisy run, not a discrimination failure: specificity remains perfect, and emergence is flagged on the large majority of runs.

## 3. Signature-independence audit

If the four signatures were near-perfectly correlated under ALL processes, convergence would be near-automatic and the criterion trivial. The audit below shows they are not.

### Correlation of the four signatures on genuine-emergence slices

```
                   neutral_departure  seriability   fst  spatial_boundary
neutral_departure               1.00         0.91  0.97              0.85
seriability                     0.91         1.00  0.90              0.75
fst                             0.97         0.90  1.00              0.87
spatial_boundary                0.85         0.75  0.87              1.00
```
Mean absolute off-diagonal correlation: **0.88**.

### Mean absolute off-diagonal correlation under each mimic

| mechanism | mean |r| among the four signatures |
|---|---|
| group_emergence (genuine) | 0.88 |
| aggregated_conformity | 0.21 |
| patchiness | 0.29 |
| drift_space | 0.19 |

**Interpretation.**

The four signatures are strongly correlated (mean |r| = 0.88) ONLY under genuine emergence, where a single coupled process drives between-group divergence and within-group conformity together so all four co-rise. Under the mimics the same four signatures are nearly independent (mean |r| roughly 0.19-0.29). Convergence is therefore NOT a built-in artifact of correlated metrics: the signatures move together precisely when, and only when, a genuine group-forming process couples their causes. Each mimic decouples those causes (conformity without divergence; static divergence; smooth isolation-by-distance) and so fails the convergence test on at least one signature.

## 4. Verdict

**GO** on the convergence criterion.

Convergence discriminates. At the reported seed it flags only the genuine group-level emergence process and none of the three mimics, and across the seed sweep no mimic is ever a false positive. The four signatures are independent under the mimics and co-rise only under genuine emergence, so their convergence is informative rather than automatic. The one caveat is sensitivity: the spatial-boundary signature is the noisiest of the four, so genuine emergence is flagged on the large majority but not 100% of single realizations. For a single empirical assemblage this argues for reporting all four ordinal trends and their joint pattern rather than a bare pass/fail, and for treating a near-miss on one signature as weak rather than negative evidence.
