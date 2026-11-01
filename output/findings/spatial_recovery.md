# Does the composition GP recover a known length scale?

Produced by `analyses/02_spatial/01_recovery.R` on 2026-10-01. Data simulated from the model by `analyses/49_write_stan_data.py` at the real basin geometry (29 assemblages, 10 classes, real sherd totals, max pairwise distance 127.8 km). Exponential kernel, uniform prior on log(rho) over [0.5, 300] km, true sigma 1.0, true tau 0.5, true spatial share 0.80.

## Recovery of rho

| true rho (km) | use_tau | rho median | rho 95% CI | d50 | spatial share | covers truth | rung | verdict |
|---|---|---|---|---|---|---|---|---|
| 5 | 1 | 5.2 | 3.1 to 9.6 | 3.6 | 0.91 | yes | 2 | underpowered |
| 5 | 0 | 4.7 | 2.9 to 7.2 | 3.2 | 1.00 | yes | 1 | recovered |
| 10 | 1 | 7.9 | 4.3 to 18.7 | 5.5 | 0.91 | yes | 2 | underpowered |
| 10 | 0 | 6.5 | 3.8 to 11.1 | 4.5 | 1.00 | yes | 1 | recovered |
| 20 | 1 | 18.8 | 9.0 to 72.3 | 13.0 | 0.84 | yes | 2 | ridge |
| 20 | 0 | 10.2 | 6.7 to 17.7 | 7.1 | 1.00 | **NO** | 1 | recovered |
| 40 | 1 | 16.6 | 5.9 to 135.4 | 11.5 | 0.66 | yes | 2 | underpowered |
| 40 | 0 | 5.6 | 3.0 to 9.6 | 3.9 | 1.00 | **NO** | 1 | recovered |
| 80 | 1 | 30.0 | 6.7 to 212.4 | 20.8 | 0.74 | yes | 2 | ridge |
| 80 | 0 | 6.5 | 4.1 to 10.3 | 4.5 | 1.00 | **NO** | 1 | recovered |

## Diagnostics (rule 16, every fit)

- true rho 5 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0039 | bulk ESS 327 | tail ESS 2623 | divergences 0/12000 (0.00%) | treedepth>=10 624 | E-BFMI 0.785
- true rho 5 km, use_tau 0 (rung 1): R-hat 1.0024 | bulk ESS 946 | tail ESS 1493 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.723
- true rho 10 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0213 | bulk ESS 295 | tail ESS 1764 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.763
- true rho 10 km, use_tau 0 (rung 1): R-hat 1.0060 | bulk ESS 923 | tail ESS 1707 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.751
- true rho 20 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0065 | bulk ESS 903 | tail ESS 2486 | divergences 0/12000 (0.00%) | treedepth>=10 2961 | E-BFMI 0.808
- true rho 20 km, use_tau 0 (rung 1): R-hat 1.0081 | bulk ESS 519 | tail ESS 953 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.787
- true rho 40 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0105 | bulk ESS 547 | tail ESS 1510 | divergences 0/12000 (0.00%) | treedepth>=10 80 | E-BFMI 0.797
- true rho 40 km, use_tau 0 (rung 1): R-hat 1.0034 | bulk ESS 797 | tail ESS 1415 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.742
- true rho 80 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0043 | bulk ESS 406 | tail ESS 821 | divergences 0/12000 (0.00%) | treedepth>=10 5908 | E-BFMI 0.797
- true rho 80 km, use_tau 0 (rung 1): R-hat 1.0009 | bulk ESS 1273 | tail ESS 2325 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.797

## Is the recovery USABLE? Interval widths

| true rho | use_tau | 95% CI width ratio (hi/lo) | spatial share | share 95% CI |
|---|---|---|---|---|
| 5 | 1 | **3.1x** | 0.91 | 0.60 to 1.00 (width 0.40) |
| 5 | 0 | **2.5x** | 1.00 (definitional, not estimated) | - |
| 10 | 1 | **4.4x** | 0.91 | 0.59 to 1.00 (width 0.41) |
| 10 | 0 | **2.9x** | 1.00 (definitional, not estimated) | - |
| 20 | 1 | **8.0x** | 0.84 | 0.66 to 0.98 (width 0.32) |
| 20 | 0 | **2.6x** | 1.00 (definitional, not estimated) | - |
| 40 | 1 | **22.8x** | 0.66 | 0.37 to 0.91 (width 0.54) |
| 40 | 0 | **3.2x** | 1.00 (definitional, not estimated) | - |
| 80 | 1 | **31.7x** | 0.74 | 0.47 to 0.96 (width 0.49) |
| 80 | 0 | **2.5x** | 1.00 (definitional, not estimated) | - |

## The nugget-versus-range trade-off

Across a true-rho range spanning a factor of 16, the posterior medians span a factor of 5.7 with the non-spatial term present and 2.2 without it.

../mataa's `moai_site_gp.stan` header warns that its tau-plus-field decomposition returned d50 of 7 to 13 km whatever the truth. **That warning inverts here.** The FULL model tracks the truth; it is the FIELD-ONLY model (use_tau = 0) that pins near 8 to 10 km regardless, because without a nugget the field must absorb every idiosyncratic site difference and is forced short. Do not carry their diagnostic reading across: on this design, dropping tau creates the pathology rather than exposing it.


## Saturation: does the estimate stop moving?

| truth pair | median pair | ratio | doubling recovered? |
|---|---|---|---|
| 5 -> 10 km | 5.2 -> 7.9 | 1.52 | yes |
| 10 -> 20 km | 7.9 -> 18.8 | 2.37 | yes |
| 20 -> 40 km | 18.8 -> 16.6 | 0.88 | **NO** |
| 40 -> 80 km | 16.6 -> 30.0 | 1.80 | yes |

**Ceiling at roughly 20 km.** Truths of 20 and 40 km return posterior medians of 18.8 and 16.6, which are the same answer. Above that the model cannot tell one long interaction scale from another.

**This is a second, independent witness for the ceiling.** `analyses/48_length_scale_recovery.py` located it near 32 km from summary statistics on drift simulations, with no GP anywhere in it; this locates it from a fitted latent field on data simulated from the model itself. Two methods sharing no machinery agree, which rule 10 says to treat as signal rather than coincidence. The common cause is geometric: the basin is 128 km across, so a decay longer than a few tens of km is not expressed inside the study window at all.

## Verdict

**rho is directionally recovered but only weakly identified.** Medians track the truth across a factor of 16 and every interval covers it, but the intervals span factors of 3.1 to 31.7. Coverage here is bought with vagueness. rho is reportable as an order-of-magnitude statement, not as a point estimate with a useful interval, which is the same outcome ../mataa recorded for pukao.

**`spatial_share` is recovered more usefully, but it is not sharp everywhere.** Its posterior MEDIAN sits within 0.135 of the truth of 0.80 at every length scale tested, which is genuine accuracy. Its INTERVAL, however, tightens with the length scale: width 0.40 at the shortest scale, where [0.60, 1.00] is close to no information at all on a bounded [0, 1] quantity, narrowing to 0.49 at the longest. Reporting a single width bound here would be misleading, and an earlier draft of this verdict did exactly that.

The variance partition is still the better-behaved quantity and the one the paper's argument needs, since it asks what fraction of assemblage-level compositional variance is organised by geography. But it earns a headline only where the interval is informative, which on this evidence means the model must be fitted with enough spatial signal present, and the width has to be reported alongside the median every time.

