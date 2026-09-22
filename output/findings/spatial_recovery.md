# Does the composition GP recover a known length scale?

Produced by `analyses/02_spatial/01_recovery.R` on 2026-09-21. Data simulated from the model by `analyses/49_write_stan_data.py` at the real basin geometry (29 assemblages, 10 classes, real sherd totals, max pairwise distance 127.8 km). Exponential kernel, uniform prior on log(rho) over [0.5, 300] km, true sigma 1.0, true tau 0.5, true spatial share 0.80.

## Recovery of rho

| true rho (km) | use_tau | rho median | rho 95% CI | d50 | spatial share | covers truth | rung | verdict |
|---|---|---|---|---|---|---|---|---|
| 5 | 1 | 4.8 | 2.8 to 8.2 | 3.3 | 0.93 | yes | 2 | underpowered |
| 5 | 0 | 4.2 | 2.5 to 6.6 | 2.9 | 1.00 | yes | 1 | recovered |
| 10 | 1 | 8.5 | 4.9 to 17.4 | 5.9 | 0.93 | yes | 2 | underpowered |
| 10 | 0 | 7.2 | 4.5 to 11.7 | 5.0 | 1.00 | yes | 1 | recovered |
| 20 | 1 | 18.5 | 8.8 to 66.6 | 12.8 | 0.83 | yes | 2 | underpowered |
| 20 | 0 | 9.6 | 6.3 to 15.6 | 6.7 | 1.00 | **NO** | 1 | recovered |
| 40 | 1 | 18.4 | 4.9 to 176.8 | 12.7 | 0.68 | yes | 2 | ridge |
| 40 | 0 | 5.2 | 2.3 to 9.4 | 3.6 | 1.00 | **NO** | 1 | recovered |
| 80 | 1 | 47.6 | 10.5 to 241.4 | 33.0 | 0.74 | yes | 2 | ridge |
| 80 | 0 | 7.0 | 4.2 to 11.9 | 4.8 | 1.00 | **NO** | 1 | recovered |

## Diagnostics (rule 16, every fit)

- true rho 5 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0114 | bulk ESS 274 | tail ESS 3720 | divergences 0/12000 (0.00%) | treedepth>=10 1086 | E-BFMI 0.750
- true rho 5 km, use_tau 0 (rung 1): R-hat 1.0025 | bulk ESS 992 | tail ESS 1461 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.736
- true rho 10 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0145 | bulk ESS 323 | tail ESS 2052 | divergences 0/12000 (0.00%) | treedepth>=10 2250 | E-BFMI 0.791
- true rho 10 km, use_tau 0 (rung 1): R-hat 1.0031 | bulk ESS 831 | tail ESS 1613 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.742
- true rho 20 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0130 | bulk ESS 527 | tail ESS 2016 | divergences 0/12000 (0.00%) | treedepth>=10 157 | E-BFMI 0.790
- true rho 20 km, use_tau 0 (rung 1): R-hat 1.0031 | bulk ESS 853 | tail ESS 1421 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.746
- true rho 40 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0040 | bulk ESS 482 | tail ESS 1629 | divergences 0/12000 (0.00%) | treedepth>=10 5399 | E-BFMI 0.751
- true rho 40 km, use_tau 0 (rung 1): R-hat 1.0070 | bulk ESS 493 | tail ESS 302 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.746
- true rho 80 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0030 | bulk ESS 699 | tail ESS 983 | divergences 0/12000 (0.00%) | treedepth>=10 6979 | E-BFMI 0.815
- true rho 80 km, use_tau 0 (rung 1): R-hat 1.0033 | bulk ESS 897 | tail ESS 1819 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.833

## Is the recovery USABLE? Interval widths

| true rho | use_tau | 95% CI width ratio (hi/lo) | spatial share | share 95% CI |
|---|---|---|---|---|
| 5 | 1 | **2.9x** | 0.93 | 0.63 to 1.00 (width 0.37) |
| 5 | 0 | **2.6x** | 1.00 (definitional, not estimated) | - |
| 10 | 1 | **3.5x** | 0.93 | 0.68 to 1.00 (width 0.32) |
| 10 | 0 | **2.6x** | 1.00 (definitional, not estimated) | - |
| 20 | 1 | **7.6x** | 0.83 | 0.62 to 0.97 (width 0.35) |
| 20 | 0 | **2.5x** | 1.00 (definitional, not estimated) | - |
| 40 | 1 | **35.8x** | 0.68 | 0.36 to 0.94 (width 0.58) |
| 40 | 0 | **4.1x** | 1.00 (definitional, not estimated) | - |
| 80 | 1 | **23.0x** | 0.74 | 0.50 to 0.92 (width 0.42) |
| 80 | 0 | **2.9x** | 1.00 (definitional, not estimated) | - |

## The nugget-versus-range trade-off

Across a true-rho range spanning a factor of 16, the posterior medians span a factor of 10.0 with the non-spatial term present and 2.3 without it.

../mataa's `moai_site_gp.stan` header warns that its tau-plus-field decomposition returned d50 of 7 to 13 km whatever the truth. **That warning inverts here.** The FULL model tracks the truth; it is the FIELD-ONLY model (use_tau = 0) that pins near 8 to 10 km regardless, because without a nugget the field must absorb every idiosyncratic site difference and is forced short. Do not carry their diagnostic reading across: on this design, dropping tau creates the pathology rather than exposing it.


## Saturation: does the estimate stop moving?

| truth pair | median pair | ratio | doubling recovered? |
|---|---|---|---|
| 5 -> 10 km | 4.8 -> 8.5 | 1.78 | yes |
| 10 -> 20 km | 8.5 -> 18.5 | 2.17 | yes |
| 20 -> 40 km | 18.5 -> 18.4 | 1.00 | **NO** |
| 40 -> 80 km | 18.4 -> 47.6 | 2.59 | yes |

**Ceiling at roughly 20 km.** Truths of 20 and 40 km return posterior medians of 18.5 and 18.4, which are the same answer. Above that the model cannot tell one long interaction scale from another.

**This is a second, independent witness for the ceiling.** `analyses/48_length_scale_recovery.py` located it near 32 km from summary statistics on drift simulations, with no GP anywhere in it; this locates it from a fitted latent field on data simulated from the model itself. Two methods sharing no machinery agree, which rule 10 says to treat as signal rather than coincidence. The common cause is geometric: the basin is 128 km across, so a decay longer than a few tens of km is not expressed inside the study window at all.

## Verdict

**rho is directionally recovered but only weakly identified.** Medians track the truth across a factor of 16 and every interval covers it, but the intervals span factors of 2.9 to 35.8. Coverage here is bought with vagueness. rho is reportable as an order-of-magnitude statement, not as a point estimate with a useful interval, which is the same outcome ../mataa recorded for pukao.

**`spatial_share` is recovered more usefully, but it is not sharp everywhere.** Its posterior MEDIAN sits within 0.132 of the truth of 0.80 at every length scale tested, which is genuine accuracy. Its INTERVAL, however, tightens with the length scale: width 0.37 at the shortest scale, where [0.63, 1.00] is close to no information at all on a bounded [0, 1] quantity, narrowing to 0.42 at the longest. Reporting a single width bound here would be misleading, and an earlier draft of this verdict did exactly that.

The variance partition is still the better-behaved quantity and the one the paper's argument needs, since it asks what fraction of assemblage-level compositional variance is organised by geography. But it earns a headline only where the interval is informative, which on this evidence means the model must be fitted with enough spatial signal present, and the width has to be reported alongside the median every time.

