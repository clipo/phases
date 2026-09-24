# Coupling robustness of the shared-latent-assortment assumption

Phase 3 transmission-layer extension. The ceramic-copying simulator is driven by the parent monument-mls emergence engine. The replicator dynamics produce phi(t), the fraction of cooperation-signaling (monument) groups over time, and a coupling parameter sets how strongly that phi drives the latent assortment level of the ceramic style-copying process.

## Operating point

- Bistable point: sigma = 0.5, lambda_W = 0.5.
- Interior saddle phi_star = 0.4299 (finite and well inside (0,1)).
- Initial condition phi_0 = phi_star + 0.1 = 0.5299 (above the saddle, so phi rises toward 1: a genuine monument-emergence trajectory).
- Coupling: ceramic assortment a(t) = clip(coupling * phi(t), 0, 1), fed to both between-group divergence and within-group conformity on the bounded spatial rule.
- Aggregated over 20 sampling seeds; G = 12, N/group = 300, K = 10 types.

## Coupling sweep

Convergence trend is the OLS ordinal slope of the combined convergence score (mean of the four column-standardized signatures). frac_all_four_up is the fraction of seeds in which all four raw signature slopes are positive (the strict criterion).

| coupling | conv_trend_mean | conv_trend_sd | frac_all_four_up | trend_neutral | trend_seriability | trend_fst | trend_spatial | detectable |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.0 | -0.017 | 0.080 | 0.000 | -0.001 | 0.018 | 0.000 | -0.011 | False |
| 0.1 | 0.195 | 0.065 | 0.550 | 0.010 | 0.398 | 0.013 | 0.759 | True |
| 0.2 | 0.233 | 0.049 | 0.750 | 0.008 | 0.939 | 0.019 | 1.716 | True |
| 0.25 | 0.218 | 0.060 | 0.650 | 0.008 | 0.991 | 0.020 | 1.383 | True |
| 0.3 | 0.214 | 0.064 | 0.500 | 0.008 | 1.039 | 0.019 | -0.528 | True |
| 0.4 | 0.236 | 0.059 | 0.750 | 0.011 | 0.883 | 0.018 | 0.727 | True |
| 0.5 | 0.227 | 0.044 | 0.650 | 0.016 | 0.643 | 0.015 | 0.459 | True |
| 0.6 | 0.220 | 0.052 | 0.450 | 0.018 | 0.533 | 0.012 | 0.331 | True |
| 0.75 | 0.220 | 0.059 | 0.550 | 0.016 | 0.469 | 0.008 | -0.619 | True |
| 1.0 | 0.207 | 0.056 | 0.500 | 0.004 | 0.533 | 0.004 | -0.077 | True |

## Threshold

The mean convergence-score trend clears the detection threshold (0.1) at coupling >= 0.1. At coupling = 0 the ceramic record carries no emergence signal and the trend is ~0 by construction. The trend rises sharply to a plateau of roughly 0.20-0.23 for essentially any non-trivial coupling and does not increase further toward coupling = 1. The convergent signature is therefore an effectively binary function of coupling with a low threshold: detection requires only that ceramic-style assortment track cooperation-relevant assortment weakly, not perfectly.

## Interpretation of the empirical non-detection

The empirical application found no convergent signature in the central Mississippi Valley ceramic record. The coupling sweep shows the criterion detects emergence whenever ceramic-style assortment tracks cooperation-relevant assortment even weakly (coupling as low as 0.1). Because the detectable range is broad and the threshold low, a non-detection is hard to attribute to a merely weak ceramic-monument coupling. Two readings remain: (1) no group-level emergence occurred over the period, or (2) emergence occurred but ceramic style was almost entirely decoupled from the cooperation-relevant assortment carried by monuments (coupling near zero). The plateau result rules out the intermediate excuse that a moderate coupling could hide a real signal.

Honesty note: the strict four-signature criterion (frac_all_four_up) is noisier than the combined convergence-score trend. The spatial-boundary signature is the volatile component across sampling seeds, so the combined score is the more reliable detector at low coupling.
