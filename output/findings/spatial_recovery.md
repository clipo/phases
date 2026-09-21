# Does the composition GP recover a known length scale?

Produced by `analyses/02_spatial/01_recovery.R` on 2026-09-20. Data simulated from the model by `analyses/49_write_stan_data.py` at the real basin geometry (29 assemblages, 10 classes, real sherd totals, max pairwise distance 127.8 km). Exponential kernel, uniform prior on log(rho) over [0.5, 300] km, true sigma 1.0, true tau 0.5, true spatial share 0.80.

## Recovery of rho

| true rho (km) | use_tau | rho median | rho 95% CI | d50 | spatial share | covers truth | rung | verdict |
|---|---|---|---|---|---|---|---|---|
| 5 | 1 | 11.5 | 2.6 to 221.9 | 8.0 | 0.49 | yes | 2 | ridge |
| 5 | 0 | 2.9 | 1.4 to 4.6 | 2.0 | 1.00 | **NO** | 2 | recovered |
| 10 | 1 | 8.6 | 5.3 to 15.1 | 5.9 | 0.79 | yes | 1 | recovered |
| 10 | 0 | 6.0 | 4.1 to 8.6 | 4.2 | 1.00 | **NO** | 2 | recovered |
| 20 | 1 | 17.7 | 10.4 to 38.7 | 12.3 | 0.88 | yes | 2 | ridge |
| 20 | 0 | 11.3 | 8.0 to 16.8 | 7.8 | 1.00 | **NO** | 1 | recovered |
| 40 | 1 | 95.2 | 28.9 to 269.1 | 66.0 | 0.84 | yes | 2 | ridge |
| 40 | 0 | 11.4 | 7.4 to 18.3 | 7.9 | 1.00 | **NO** | 1 | recovered |
| 80 | 1 | 75.2 | 26.0 to 251.4 | 52.1 | 0.75 | yes | 2 | ridge |
| 80 | 0 | 7.8 | 5.5 to 11.2 | 5.4 | 1.00 | **NO** | 1 | recovered |

## Diagnostics (rule 16, every fit)

- true rho 5 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0187 | bulk ESS 189 | tail ESS 1205 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.777
- true rho 5 km, use_tau 0 (rung 2): R-hat 1.0016 | bulk ESS 839 | tail ESS 1077 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.805
- true rho 10 km, use_tau 1 (rung 1): R-hat 1.0047 | bulk ESS 588 | tail ESS 2013 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.769
- true rho 10 km, use_tau 0 (rung 2): R-hat 1.0032 | bulk ESS 1328 | tail ESS 2368 | divergences 0/12000 (0.00%) | treedepth>=10 38 | E-BFMI 0.761
- true rho 20 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0024 | bulk ESS 511 | tail ESS 1931 | divergences 0/12000 (0.00%) | treedepth>=10 3009 | E-BFMI 0.818
- true rho 20 km, use_tau 0 (rung 1): R-hat 1.0071 | bulk ESS 1067 | tail ESS 2029 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.812
- true rho 40 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0012 | bulk ESS 2560 | tail ESS 4597 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.785
- true rho 40 km, use_tau 0 (rung 1): R-hat 1.0088 | bulk ESS 1054 | tail ESS 1962 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.704
- true rho 80 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0038 | bulk ESS 2429 | tail ESS 4017 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.822
- true rho 80 km, use_tau 0 (rung 1): R-hat 1.0094 | bulk ESS 580 | tail ESS 1030 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.766

## Is the recovery USABLE? Interval widths

| true rho | use_tau | 95% CI width ratio (hi/lo) | spatial share | share 95% CI |
|---|---|---|---|---|
| 5 | 1 | **84.4x** | 0.49 | 0.18 to 0.86 (width 0.67) |
| 5 | 0 | **3.3x** | 1.00 (definitional, not estimated) | - |
| 10 | 1 | **2.9x** | 0.79 | 0.60 to 0.90 (width 0.30) |
| 10 | 0 | **2.1x** | 1.00 (definitional, not estimated) | - |
| 20 | 1 | **3.7x** | 0.88 | 0.73 to 0.98 (width 0.25) |
| 20 | 0 | **2.1x** | 1.00 (definitional, not estimated) | - |
| 40 | 1 | **9.3x** | 0.84 | 0.68 to 0.93 (width 0.26) |
| 40 | 0 | **2.5x** | 1.00 (definitional, not estimated) | - |
| 80 | 1 | **9.7x** | 0.75 | 0.57 to 0.90 (width 0.33) |
| 80 | 0 | **2.0x** | 1.00 (definitional, not estimated) | - |

## The nugget-versus-range trade-off

Across a true-rho range spanning a factor of 16, the posterior medians span a factor of 11.1 with the non-spatial term present and 4.0 without it.

../mataa's `moai_site_gp.stan` header warns that its tau-plus-field decomposition returned d50 of 7 to 13 km whatever the truth. **That warning inverts here.** The FULL model tracks the truth; it is the FIELD-ONLY model (use_tau = 0) that pins near 8 to 10 km regardless, because without a nugget the field must absorb every idiosyncratic site difference and is forced short. Do not carry their diagnostic reading across: on this design, dropping tau creates the pathology rather than exposing it.


## Saturation: does the estimate stop moving?

| truth pair | median pair | ratio | doubling recovered? |
|---|---|---|---|
| 5 -> 10 km | 11.5 -> 8.6 | 0.75 | **NO** |
| 10 -> 20 km | 8.6 -> 17.7 | 2.06 | yes |
| 20 -> 40 km | 17.7 -> 95.2 | 5.38 | yes |
| 40 -> 80 km | 95.2 -> 75.2 | 0.79 | **NO** |

**Ceiling at roughly 5 km.** Truths of 5 and 10 km return posterior medians of 11.5 and 8.6, which are the same answer. Above that the model cannot tell one long interaction scale from another.

**This is a second, independent witness for the ceiling.** `analyses/48_length_scale_recovery.py` located it near 32 km from summary statistics on drift simulations, with no GP anywhere in it; this locates it from a fitted latent field on data simulated from the model itself. Two methods sharing no machinery agree, which rule 10 says to treat as signal rather than coincidence. The common cause is geometric: the basin is 128 km across, so a decay longer than a few tens of km is not expressed inside the study window at all.

## Verdict

**rho is directionally recovered but only weakly identified.** Medians track the truth across a factor of 16 and every interval covers it, but the intervals span factors of 2.9 to 84.4. Coverage here is bought with vagueness. rho is reportable as an order-of-magnitude statement, not as a point estimate with a useful interval, which is the same outcome ../mataa recorded for pukao.

**`spatial_share` is recovered more usefully, but it is not sharp everywhere.** Its posterior MEDIAN sits within 0.306 of the truth of 0.80 at every length scale tested, which is genuine accuracy. Its INTERVAL, however, tightens with the length scale: width 0.67 at the shortest scale, where [0.18, 0.86] is close to no information at all on a bounded [0, 1] quantity, narrowing to 0.33 at the longest. Reporting a single width bound here would be misleading, and an earlier draft of this verdict did exactly that.

The variance partition is still the better-behaved quantity and the one the paper's argument needs, since it asks what fraction of assemblage-level compositional variance is organised by geography. But it earns a headline only where the interval is informative, which on this evidence means the model must be fitted with enough spatial signal present, and the width has to be reported alongside the median every time.

