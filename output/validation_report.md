# Phase 4: Criterion Validation by Simulation

This report tests the paper's core claim: that CONVERGENCE of four cultural-transmission signatures discriminates genuine group-level emergence from mimics. The harness applies an identical signature pipeline to every mechanism, blind to its name; a mechanism is flagged 'convergent' iff all four signatures show a positive ordinal trend whose standardized slope exceeds 0.1 (standard deviations of the signature per ordinal step).

The chronology is ORDINAL (relative phases), so 'trajectory' means a monotonic trend across the ordinal slice index, not a calendar rate.

## 1. Four-signature panels per mechanism

Panels shown for seed 42. Columns: neutral_departure, seriability, fst, spatial_boundary.

### group_emergence

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.351        -30.0  0.003            -0.842
1              0.637        -20.0  0.617           132.717
2              0.739        -18.0  0.796           151.168
3              0.793        -11.0  0.902           155.443
4              0.874        -12.0  0.966           140.493
5              0.915         -5.0  0.980           189.572
6              0.942         -2.0  0.992           161.234
7              0.959         -3.0  0.997           161.528
```

### aggregated_signaling

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.303        -29.0  0.003             2.167
1              0.628        -25.0  0.005             1.452
2              0.735        -17.0  0.003            -0.738
3              0.723        -10.0  0.003            -0.190
4              0.860         -7.0  0.004             0.667
5              0.889         -5.0  0.003            -0.167
6              0.954         -6.0  0.002            -0.238
7              0.949         -5.0  0.003             0.214
```

### patchiness

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.351        -29.0  0.304            64.808
1              0.312        -27.0  0.293            62.967
2              0.329        -27.0  0.293            67.735
3              0.313        -30.0  0.295            64.935
4              0.270        -34.0  0.279            67.335
5              0.293        -32.0  0.289            65.160
6              0.324        -29.0  0.301            67.399
7              0.290        -28.0  0.297            67.321
```

### drift_space

```
   neutral_departure  seriability    fst  spatial_boundary
0              0.295        -29.0  0.239             0.747
1              0.283        -31.0  0.234             3.294
2              0.338        -31.0  0.253            -1.653
3              0.287        -34.0  0.233             0.108
4              0.313        -33.0  0.251             0.822
5              0.321        -31.0  0.252             0.936
6              0.314        -32.0  0.236             2.175
7              0.284        -28.0  0.241            -1.186
```

## 2. Discrimination verdict

Threshold on standardized ordinal slope: 0.1.

| mechanism | neutral_departure | seriability | fst | spatial_boundary | CONVERGENT |
|---|---|---|---|---|---|
| group_emergence | +0.300 | +0.325 | +0.343 | +0.292 | **True** |
| aggregated_signaling | +0.289 | +0.307 | +0.000 | -0.002 | **False** |
| patchiness | -0.006 | +0.041 | +0.000 | -0.001 | **False** |
| drift_space | -0.007 | -0.002 | +0.002 | -0.000 | **False** |

Each cell is the standardized ordinal slope of that signature. A mechanism is convergent only when all four are above threshold.

### Cross-seed robustness

Across 500 seeds: genuine emergence flagged convergent in 493/500 runs (sensitivity); mimics flagged convergent in 0 runs (false positives).

Degenerate runs: 2/500 ([187, 264]). A run is degenerate when a generator produces an ordinal slice on which a signature is undefined, which here means strong conformity fixing a single class so that cultural F_ST has no total diversity to partition. That is the mimic behaving as designed. Such a run cannot be scored as convergent and is reported here rather than folded into the rates above.

No mimic was ever flagged convergent. Specificity is the load-bearing property of the criterion and it holds at 100%.

Sensitivity is below 100%: on an occasional single realization the spatial-boundary signature carries enough sampling noise that its ordinal trend dips below threshold even though the other three rise. This is a power limit on a single noisy run, not a discrimination failure: specificity remains perfect, and emergence is flagged on the large majority of runs.

## 3. Signature-independence audit

If the four signatures were near-perfectly correlated under ALL processes, convergence would be near-automatic and the criterion trivial. The audit below shows they are not.

### Correlation of the four signatures on genuine-emergence slices

```
                   neutral_departure  seriability   fst  spatial_boundary
neutral_departure               1.00         0.92  0.97              0.85
seriability                     0.92         1.00  0.90              0.75
fst                             0.97         0.90  1.00              0.87
spatial_boundary                0.85         0.75  0.87              1.00
```
Mean absolute off-diagonal correlation: **0.88**.

### Mean absolute off-diagonal correlation under each mimic

| mechanism | mean |r| among the four signatures |
|---|---|
| group_emergence (genuine) | 0.88 |
| aggregated_signaling | 0.21 |
| patchiness | 0.28 |
| drift_space | 0.18 |

**Interpretation.**

The four signatures are strongly correlated (mean |r| = 0.88) ONLY under genuine emergence, where a single coupled process drives between-group divergence and within-group conformity together so all four co-rise. Under the mimics the same four signatures are nearly independent (mean |r| roughly 0.18-0.28). Convergence is therefore NOT a built-in artifact of correlated metrics: the signatures move together precisely when, and only when, a genuine group-forming process couples their causes. Each mimic decouples those causes (conformity without divergence; static divergence; smooth isolation-by-distance) and so fails the convergence test on at least one signature.

## 4. Verdict

**GO** on the convergence criterion.

Convergence discriminates. At the reported seed it flags only the genuine group-level emergence process and none of the three mimics, and across the seed sweep no mimic is ever a false positive. The four signatures are independent under the mimics and co-rise only under genuine emergence, so their convergence is informative rather than automatic. The one caveat is sensitivity: the spatial-boundary signature is the noisiest of the four, so genuine emergence is flagged on the large majority but not 100% of single realizations. For a single empirical assemblage this argues for reporting all four ordinal trends and their joint pattern rather than a bare pass/fail, and for treating a near-miss on one signature as weak rather than negative evidence.
