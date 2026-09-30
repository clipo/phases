# Calibrated spatial drift comparison

Two neutral innovation models, each calibrated to observed within-assemblage Gini-Simpson diversity, richness, and pooled diversity before comparison. Fixed observed repertoire; full eight-slice windows; 1,200-generation burn-in.

## Calibration

| Region | Model | N | Innovation | Mixing | Matched cells | Obs H_S / rich / H_T | Achieved (baseline median) |
|---|---|---:|---:|---:|---:|---|---|
| basin | pooled | 10000 | 0.0002 | 0.005 | 61/297 | 0.510 / 5.90 / 0.522 | 0.512 / 5.31 / 0.523 |
| basin | uniform | 120 | 0.0002 | 0.200 | 1/297 | 0.510 / 5.90 / 0.522 | 0.599 / 5.83 / 0.629 |
| valley | pooled | 2000 | 0.0020 | 0.005 | 71/297 | 0.520 / 6.02 / 0.532 | 0.515 / 5.51 / 0.532 |
| valley | uniform | 120 | 0.0002 | 0.100 | 0/297 (unmatched; min loss) | 0.520 / 6.02 / 0.532 | 0.714 / 7.05 / 0.751 |
| cmv | pooled | 120 | 0.0120 | 0.005 | 10/297 | 0.308 / 5.03 / 0.328 | 0.294 / 4.95 / 0.324 |
| cmv | uniform | 120 | 0.0002 | 0.005 | 0/297 (unmatched; min loss) | 0.308 / 5.03 / 0.328 | 0.527 / 6.18 / 0.846 |

## Baseline comparisons (upper-tail Monte Carlo p, add-one)

| Model | Region | Sampling | Statistic | Observed | Median [95%] | p upper | p lower |
|---|---|---|---|---:|---|---:|---:|
| pooled | basin | contemporaneous | community_fst | 0.0430 | 0.0000 [-0.0000, 0.0107] | 0.0020 | 1.0000 |
| pooled | basin | contemporaneous | spatial_fst | 0.0179 | 0.0032 [0.0006, 0.0147] | 0.0140 | 0.9880 |
| pooled | basin | contemporaneous | ncom | 2.0000 | 1.0000 [1.0000, 2.0000] | 0.0459 | 1.0000 |
| pooled | basin | reversed | community_fst | 0.0430 | 0.0000 [0.0000, 0.0087] | 0.0020 | 1.0000 |
| pooled | basin | reversed | spatial_fst | 0.0179 | 0.0028 [0.0004, 0.0129] | 0.0080 | 0.9940 |
| pooled | basin | reversed | ncom | 2.0000 | 1.0000 [1.0000, 2.0000] | 0.0299 | 1.0000 |
| pooled | basin | time_transgressive | community_fst | 0.0430 | 0.0000 [0.0000, 0.0000] | 0.0020 | 1.0000 |
| pooled | basin | time_transgressive | spatial_fst | 0.0179 | 0.0029 [0.0005, 0.0128] | 0.0120 | 0.9900 |
| pooled | basin | time_transgressive | ncom | 2.0000 | 1.0000 [1.0000, 1.0000] | 0.0259 | 1.0000 |
| pooled | cmv | contemporaneous | community_fst | 0.0361 | 0.0290 [0.0000, 0.0573] | 0.3333 | 0.6687 |
| pooled | cmv | contemporaneous | spatial_fst | 0.0337 | 0.0125 [0.0043, 0.0353] | 0.0419 | 0.9601 |
| pooled | cmv | contemporaneous | ncom | 2.0000 | 2.0000 [1.0000, 2.0000] | 0.7485 | 1.0000 |
| pooled | cmv | reversed | community_fst | 0.0361 | 0.0288 [0.0000, 0.0598] | 0.3174 | 0.6846 |
| pooled | cmv | reversed | spatial_fst | 0.0337 | 0.0122 [0.0037, 0.0325] | 0.0220 | 0.9800 |
| pooled | cmv | reversed | ncom | 2.0000 | 2.0000 [1.0000, 2.0000] | 0.7465 | 1.0000 |
| pooled | cmv | time_transgressive | community_fst | 0.0361 | 0.0303 [0.0000, 0.0654] | 0.3413 | 0.6607 |
| pooled | cmv | time_transgressive | spatial_fst | 0.0337 | 0.0121 [0.0039, 0.0319] | 0.0240 | 0.9780 |
| pooled | cmv | time_transgressive | ncom | 2.0000 | 2.0000 [1.0000, 2.0000] | 0.7705 | 1.0000 |
| pooled | valley | contemporaneous | community_fst | 0.0413 | 0.0118 [0.0000, 0.0253] | 0.0020 | 1.0000 |
| pooled | valley | contemporaneous | spatial_fst | 0.0094 | 0.0004 [0.0000, 0.0025] | 0.0020 | 1.0000 |
| pooled | valley | contemporaneous | ncom | 2.0000 | 2.0000 [1.0000, 2.0000] | 0.9721 | 1.0000 |
| pooled | valley | contemporaneous | parkin_fst | 0.0121 | 0.0028 [0.0003, 0.0145] | 0.0479 | 0.9541 |
| pooled | valley | reversed | community_fst | 0.0413 | 0.0120 [0.0040, 0.0242] | 0.0020 | 1.0000 |
| pooled | valley | reversed | spatial_fst | 0.0094 | 0.0004 [0.0000, 0.0022] | 0.0020 | 1.0000 |
| pooled | valley | reversed | ncom | 2.0000 | 2.0000 [2.0000, 2.0000] | 0.9800 | 1.0000 |
| pooled | valley | reversed | parkin_fst | 0.0121 | 0.0027 [0.0003, 0.0146] | 0.0499 | 0.9521 |
| pooled | valley | time_transgressive | community_fst | 0.0413 | 0.0114 [0.0046, 0.0217] | 0.0020 | 1.0000 |
| pooled | valley | time_transgressive | spatial_fst | 0.0094 | 0.0004 [0.0001, 0.0022] | 0.0020 | 1.0000 |
| pooled | valley | time_transgressive | ncom | 2.0000 | 2.0000 [2.0000, 2.0000] | 0.9880 | 1.0000 |
| pooled | valley | time_transgressive | parkin_fst | 0.0121 | 0.0023 [0.0003, 0.0122] | 0.0299 | 0.9721 |
| uniform | basin | contemporaneous | community_fst | 0.0430 | 0.0000 [-0.0000, 0.0188] | 0.0020 | 1.0000 |
| uniform | basin | contemporaneous | spatial_fst | 0.0179 | 0.0065 [0.0013, 0.0215] | 0.0659 | 0.9361 |
| uniform | basin | contemporaneous | ncom | 2.0000 | 1.0000 [1.0000, 2.0000] | 0.2455 | 0.9940 |
| uniform | basin | reversed | community_fst | 0.0430 | 0.0258 [0.0000, 0.0898] | 0.1976 | 0.8044 |
| uniform | basin | reversed | spatial_fst | 0.0179 | 0.0123 [0.0024, 0.0429] | 0.3194 | 0.6826 |
| uniform | basin | reversed | ncom | 2.0000 | 2.0000 [1.0000, 3.0000] | 0.9521 | 0.9681 |
| uniform | basin | time_transgressive | community_fst | 0.0430 | 0.0277 [0.0000, 0.1109] | 0.2655 | 0.7365 |
| uniform | basin | time_transgressive | spatial_fst | 0.0179 | 0.0157 [0.0033, 0.0625] | 0.4351 | 0.5669 |
| uniform | basin | time_transgressive | ncom | 2.0000 | 2.0000 [1.0000, 3.0000] | 0.9541 | 0.9561 |
| uniform | cmv | contemporaneous | community_fst | 0.0361 | 0.1556 [0.0841, 0.2461] | 1.0000 | 0.0020 |
| uniform | cmv | contemporaneous | spatial_fst | 0.0337 | 0.0849 [0.0386, 0.1683] | 0.9880 | 0.0140 |
| uniform | cmv | contemporaneous | ncom | 2.0000 | 3.0000 [2.0000, 4.0000] | 1.0000 | 0.2136 |
| uniform | cmv | reversed | community_fst | 0.0361 | 0.1526 [0.0796, 0.2476] | 1.0000 | 0.0020 |
| uniform | cmv | reversed | spatial_fst | 0.0337 | 0.0839 [0.0398, 0.1529] | 0.9900 | 0.0120 |
| uniform | cmv | reversed | ncom | 2.0000 | 3.0000 [2.0000, 5.0000] | 1.0000 | 0.1417 |
| uniform | cmv | time_transgressive | community_fst | 0.0361 | 0.1532 [0.0802, 0.2550] | 1.0000 | 0.0020 |
| uniform | cmv | time_transgressive | spatial_fst | 0.0337 | 0.0832 [0.0443, 0.1451] | 0.9960 | 0.0060 |
| uniform | cmv | time_transgressive | ncom | 2.0000 | 3.0000 [2.0000, 5.0000] | 1.0000 | 0.1317 |
| uniform | valley | contemporaneous | community_fst | 0.0413 | 0.0125 [0.0000, 0.0447] | 0.0319 | 0.9701 |
| uniform | valley | contemporaneous | spatial_fst | 0.0094 | 0.0032 [0.0004, 0.0120] | 0.0739 | 0.9281 |
| uniform | valley | contemporaneous | ncom | 2.0000 | 2.0000 [1.0000, 2.0000] | 0.9721 | 0.9960 |
| uniform | valley | contemporaneous | parkin_fst | 0.0121 | 0.0131 [0.0022, 0.0511] | 0.5429 | 0.4591 |
| uniform | valley | reversed | community_fst | 0.0413 | 0.0209 [0.0078, 0.0543] | 0.0878 | 0.9142 |
| uniform | valley | reversed | spatial_fst | 0.0094 | 0.0037 [0.0008, 0.0121] | 0.0699 | 0.9321 |
| uniform | valley | reversed | ncom | 2.0000 | 2.0000 [2.0000, 3.0000] | 1.0000 | 0.9541 |
| uniform | valley | reversed | parkin_fst | 0.0121 | 0.0164 [0.0032, 0.0471] | 0.6886 | 0.3134 |
| uniform | valley | time_transgressive | community_fst | 0.0413 | 0.0202 [0.0067, 0.0491] | 0.0699 | 0.9321 |
| uniform | valley | time_transgressive | spatial_fst | 0.0094 | 0.0034 [0.0005, 0.0114] | 0.0659 | 0.9361 |
| uniform | valley | time_transgressive | ncom | 2.0000 | 2.0000 [2.0000, 3.0000] | 1.0000 | 0.9681 |
| uniform | valley | time_transgressive | parkin_fst | 0.0121 | 0.0138 [0.0027, 0.0480] | 0.5768 | 0.4251 |

## Boundary grid (cells bracketing observed, of 2 lengths)

| Model | Region | Leak | Bracketing / total |
|---|---|---:|---:|
| pooled | basin | 1 | 0 / 2 |
| pooled | basin | 0.5 | 2 / 2 |
| pooled | basin | 0.1 | 1 / 2 |
| pooled | basin | 0.03 | 2 / 2 |
| pooled | valley | 1 | 1 / 2 |
| pooled | valley | 0.5 | 1 / 2 |
| pooled | valley | 0.1 | 1 / 2 |
| pooled | valley | 0.03 | 1 / 2 |
| pooled | cmv | 1 | 2 / 2 |
| pooled | cmv | 0.5 | 1 / 2 |
| pooled | cmv | 0.1 | 1 / 2 |
| pooled | cmv | 0.03 | 2 / 2 |
| uniform | basin | 1 | 1 / 2 |
| uniform | basin | 0.5 | 1 / 2 |
| uniform | basin | 0.1 | 0 / 2 |
| uniform | basin | 0.03 | 1 / 2 |
| uniform | valley | 1 | 1 / 2 |
| uniform | valley | 0.5 | 1 / 2 |
| uniform | valley | 0.1 | 0 / 2 |
| uniform | valley | 0.03 | 0 / 2 |
| uniform | cmv | 1 | 0 / 2 |
| uniform | cmv | 0.5 | 0 / 2 |
| uniform | cmv | 0.1 | 0 / 2 |
| uniform | cmv | 0.03 | 0 / 2 |

## Sensitivity (50 realizations per case)

| Model | Region | Case | Median [95%] | p upper | mean H_S | mean richness |
|---|---|---|---|---:|---:|---:|
| pooled | basin | innovation_doubled | 0.0029 [0.0009, 0.0122] | 0.020 | 0.509 | 5.56 |
| pooled | basin | innovation_halved | 0.0028 [0.0007, 0.0160] | 0.039 | 0.502 | 5.04 |
| pooled | basin | longer_burnin | 0.0030 [0.0006, 0.0107] | 0.020 | 0.515 | 5.19 |
| pooled | basin | mixing_doubled | 0.0019 [0.0002, 0.0071] | 0.020 | 0.516 | 5.65 |
| pooled | basin | mixing_halved | 0.0039 [0.0006, 0.0211] | 0.059 | 0.501 | 4.92 |
| pooled | basin | monomorphic_burnin | 0.0028 [0.0006, 0.0126] | 0.020 | 0.223 | 4.23 |
| pooled | basin | monomorphic_no_burnin | 0.0033 [0.0004, 0.0141] | 0.020 | 0.106 | 3.64 |
| pooled | valley | innovation_doubled | 0.0012 [0.0002, 0.0098] | 0.039 | 0.521 | 5.90 |
| pooled | valley | innovation_halved | 0.0045 [0.0006, 0.0154] | 0.078 | 0.512 | 5.13 |
| pooled | valley | longer_burnin | 0.0024 [0.0006, 0.0152] | 0.078 | 0.517 | 5.56 |
| pooled | valley | mixing_doubled | 0.0022 [0.0004, 0.0089] | 0.020 | 0.516 | 5.88 |
| pooled | valley | mixing_halved | 0.0030 [0.0005, 0.0134] | 0.078 | 0.506 | 5.19 |
| pooled | valley | monomorphic_burnin | 0.0034 [0.0003, 0.0123] | 0.078 | 0.505 | 5.48 |
| pooled | valley | monomorphic_no_burnin | 0.0035 [0.0004, 0.0272] | 0.157 | 0.409 | 5.15 |
| uniform | basin | innovation_doubled | 0.0135 [0.0026, 0.0399] | 0.333 | 0.645 | 7.35 |
| uniform | basin | innovation_halved | 0.0143 [0.0041, 0.0655] | 0.412 | 0.490 | 4.57 |
| uniform | basin | longer_burnin | 0.0165 [0.0019, 0.0600] | 0.431 | 0.538 | 5.64 |
| uniform | basin | mixing_doubled | 0.0131 [0.0019, 0.0550] | 0.333 | 0.593 | 6.05 |
| uniform | basin | mixing_halved | 0.0203 [0.0051, 0.0655] | 0.588 | 0.566 | 5.73 |
| uniform | basin | monomorphic_burnin | 0.0134 [0.0021, 0.0549] | 0.412 | 0.364 | 5.05 |
| uniform | basin | monomorphic_no_burnin | 0.0096 [0.0012, 0.0417] | 0.373 | 0.204 | 4.26 |
| uniform | valley | innovation_doubled | 0.0135 [0.0046, 0.0402] | 0.588 | 0.746 | 8.25 |
| uniform | valley | innovation_halved | 0.0136 [0.0015, 0.0364] | 0.529 | 0.646 | 5.95 |
| uniform | valley | longer_burnin | 0.0140 [0.0024, 0.0722] | 0.529 | 0.631 | 6.62 |
| uniform | valley | mixing_doubled | 0.0079 [0.0027, 0.0370] | 0.392 | 0.697 | 7.25 |
| uniform | valley | mixing_halved | 0.0246 [0.0041, 0.0922] | 0.824 | 0.682 | 6.60 |
| uniform | valley | monomorphic_burnin | 0.0106 [0.0011, 0.0593] | 0.431 | 0.434 | 5.60 |
| uniform | valley | monomorphic_no_burnin | 0.0088 [0.0012, 0.0779] | 0.471 | 0.223 | 4.56 |

## Parkin leave-one-out (observation removed from observed and simulated counts)

| Model | Site | Observed | Median [95%] | p upper |
|---|---|---:|---|---:|
| pooled | Barton_Ranch | 0.0114 | 0.0024 [0.0003, 0.0130] | 0.038 |
| pooled | Fortune | 0.0116 | 0.0023 [0.0003, 0.0123] | 0.030 |
| pooled | Neeleys_Ferry | 0.0114 | 0.0024 [0.0003, 0.0126] | 0.034 |
| pooled | Parkin | 0.0108 | 0.0020 [0.0003, 0.0119] | 0.034 |
| pooled | Rose_Mound | 0.0134 | 0.0024 [0.0003, 0.0127] | 0.026 |
| pooled | Turnbow | 0.0127 | 0.0023 [0.0003, 0.0124] | 0.026 |
| pooled | Vernon_Paul | 0.0104 | 0.0024 [0.0003, 0.0131] | 0.052 |
| pooled | Williamson | 0.0111 | 0.0024 [0.0003, 0.0119] | 0.032 |
| uniform | Barton_Ranch | 0.0114 | 0.0130 [0.0024, 0.0463] | 0.581 |
| uniform | Fortune | 0.0116 | 0.0138 [0.0027, 0.0486] | 0.595 |
| uniform | Neeleys_Ferry | 0.0114 | 0.0119 [0.0022, 0.0411] | 0.517 |
| uniform | Parkin | 0.0108 | 0.0128 [0.0026, 0.0406] | 0.591 |
| uniform | Rose_Mound | 0.0134 | 0.0155 [0.0031, 0.0517] | 0.585 |
| uniform | Turnbow | 0.0127 | 0.0142 [0.0028, 0.0487] | 0.561 |
| uniform | Vernon_Paul | 0.0104 | 0.0128 [0.0024, 0.0466] | 0.611 |
| uniform | Williamson | 0.0111 | 0.0136 [0.0026, 0.0477] | 0.613 |
