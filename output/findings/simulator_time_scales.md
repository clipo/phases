# Simulator time scales, burn-in and stationarity (pilot step 5; D22 operating points)

Produced by `analyses/probe_latent_groups/time_scales.py`; estimator in `mls_emergence.probe.timescales` (spec 3.3, with the S/10 window rule). Exploratory design validation, not evidence about the record. Points from `operating_points.json`; seeds [1, 2, 3, 4, 5]; first stored length 5000 after burn-in 5000.

## primary: gaussian-6km, n_ind 2000, innovation 0.002, mixing 0.4

Burn-in used 5000 (repeated with 5 L_R: False); stored length at which every curve was accepted: 5000 after 0 doublings.

| scale | per seed | median | series included / excluded |
|---|---|---|---|
| L_d (local) | [194, 177, 173, 231, 182] | **182** | [269, 270, 274, 280, 273] / [11, 10, 6, 0, 7] |
| L_R (regional) | [294, 245, 261, 318, 292] | **292** | classes excluded [1, 1, 0, 0, 1] |

Probe settings that follow: stored length max(2 L_d + L_R, 10 L_R) = 2920; burn-in max(5 L_R, 5000) = 5000.

| scenario | first-half dispersion per seed | second-half | difference of means | range across seeds | stationary |
|---|---|---|---|---|---|
| S0 | [0.00123, 0.00057, 0.00065, 0.00064, 0.00049] | [0.00111, 0.00052, 0.00044, 0.00098, 0.00074] | 0.00004 | 0.00074 | yes |
| S2 | [0.00065, 0.0006, 0.00132, 0.00091, 0.00065] | [0.00079, 0.00053, 0.00111, 0.0007, 0.00049] | 0.00010 | 0.00072 | yes |

Figures: `time_scales_primary_acf.png`, `time_scales_primary_regional.png`.

## secondary: exponential-6km, n_ind 120, innovation 0.024, mixing 0.4

Burn-in used 5000 (repeated with 5 L_R: False); stored length at which every curve was accepted: 5000 after 0 doublings.

| scale | per seed | median | series included / excluded |
|---|---|---|---|
| L_d (local) | [12, 10, 10, 13, 14] | **12** | [280, 280, 280, 280, 280] / [0, 0, 0, 0, 0] |
| L_R (regional) | [40, 39, 37, 44, 42] | **40** | classes excluded [0, 0, 0, 0, 0] |

Probe settings that follow: stored length max(2 L_d + L_R, 10 L_R) = 400; burn-in max(5 L_R, 5000) = 5000.

| scenario | first-half dispersion per seed | second-half | difference of means | range across seeds | stationary |
|---|---|---|---|---|---|
| S0 | [0.00108, 0.00094, 0.00144, 0.00105, 0.001] | [0.0011, 0.00102, 0.00138, 0.00106, 0.00109] | 0.00003 | 0.00051 | yes |
| S2 | [0.00215, 0.00154, 0.00113, 0.00117, 0.00125] | [0.0021, 0.0016, 0.00146, 0.00154, 0.00136] | 0.00017 | 0.00102 | yes |

Figures: `time_scales_secondary_acf.png`, `time_scales_secondary_regional.png`.

## inherited_not_used: exponential-24km, n_ind 2000, innovation 0.001, mixing 0.02

Burn-in used 5000 (repeated with 5 L_R: False); stored length at which every curve was accepted: 5000 after 0 doublings.

| scale | per seed | median | series included / excluded |
|---|---|---|---|
| L_d (local) | [68, 95, 65, 80, 78] | **78** | [264, 280, 273, 279, 265] / [16, 0, 7, 1, 15] |
| L_R (regional) | [407, 494, 296, 360, 389] | **389** | classes excluded [2, 0, 2, 0, 1] |

Probe settings that follow: stored length max(2 L_d + L_R, 10 L_R) = 3890; burn-in max(5 L_R, 5000) = 5000.

| scenario | first-half dispersion per seed | second-half | difference of means | range across seeds | stationary |
|---|---|---|---|---|---|
| S0 | [0.00063, 0.00063, 0.00066, 0.00071, 0.00071] | [0.00072, 0.00063, 0.00068, 0.00067, 0.00062] | 0.00000 | 0.00011 | yes |
| S2 | [0.00143, 0.00109, 0.00107, 0.00095, 0.00186] | [0.00153, 0.00097, 0.00107, 0.00096, 0.00096] | 0.00018 | 0.00090 | yes |

Figures: `time_scales_inherited_not_used_acf.png`, `time_scales_inherited_not_used_regional.png`.
