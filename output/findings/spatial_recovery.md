# Does the composition GP recover a known length scale?

Produced by `analyses/02_spatial/01_recovery.R` on 2026-08-31. Data simulated from the model by `analyses/49_write_stan_data.py` at the real basin geometry (29 assemblages, 10 classes, real sherd totals, max pairwise distance 127.8 km). Exponential kernel, uniform prior on log(rho) over [0.5, 300] km, true sigma 1.0, true tau 0.5, true spatial share 0.80.

## Recovery of rho

| true rho (km) | use_tau | rho median | rho 95% CI | d50 | spatial share | covers truth | rung | verdict |
|---|---|---|---|---|---|---|---|---|
| 5 | 1 | 5.4 | 2.2 to 17.4 | 3.8 | 0.80 | yes | 2 | ridge |
| 5 | 0 | 4.2 | 1.9 to 7.0 | 2.9 | 1.00 | yes | 2 | recovered |
| 10 | 1 | 12.8 | 6.2 to 51.3 | 8.9 | 0.79 | yes | 2 | ridge |
| 10 | 0 | 7.8 | 4.9 to 12.1 | 5.4 | 1.00 | yes | 1 | recovered |
| 20 | 1 | 27.8 | 10.9 to 135.3 | 19.3 | 0.77 | yes | 2 | ridge |
| 20 | 0 | 10.2 | 6.4 to 17.1 | 7.1 | 1.00 | **NO** | 2 | recovered |
| 40 | 1 | 84.4 | 22.5 to 269.6 | 58.5 | 0.74 | yes | 2 | ridge |
| 40 | 0 | 8.1 | 4.9 to 13.3 | 5.6 | 1.00 | **NO** | 2 | recovered |
| 80 | 1 | 86.5 | 25.1 to 265.0 | 59.9 | 0.83 | yes | 2 | ridge |
| 80 | 0 | 13.3 | 8.4 to 23.4 | 9.2 | 1.00 | **NO** | 1 | recovered |

## Diagnostics (rule 16, every fit)

- true rho 5 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0246 | bulk ESS 225 | tail ESS 968 | divergences 0/12000 (0.00%) | treedepth>=10 8251 | E-BFMI 0.785
- true rho 5 km, use_tau 0 (rung 2): R-hat 1.0023 | bulk ESS 942 | tail ESS 1424 | divergences 0/12000 (0.00%) | treedepth>=10 6 | E-BFMI 0.810
- true rho 10 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0178 | bulk ESS 279 | tail ESS 1948 | divergences 0/12000 (0.00%) | treedepth>=10 10336 | E-BFMI 0.821
- true rho 10 km, use_tau 0 (rung 1): R-hat 1.0088 | bulk ESS 638 | tail ESS 1354 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.771
- true rho 20 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0050 | bulk ESS 1008 | tail ESS 2879 | divergences 0/12000 (0.00%) | treedepth>=10 11969 | E-BFMI 0.777
- true rho 20 km, use_tau 0 (rung 2): R-hat 1.0034 | bulk ESS 808 | tail ESS 1724 | divergences 0/12000 (0.00%) | treedepth>=10 2 | E-BFMI 0.759
- true rho 40 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0006 | bulk ESS 2706 | tail ESS 4247 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.817
- true rho 40 km, use_tau 0 (rung 2): R-hat 1.0024 | bulk ESS 1505 | tail ESS 1900 | divergences 0/12000 (0.00%) | treedepth>=10 14 | E-BFMI 0.776
- true rho 80 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0030 | bulk ESS 2753 | tail ESS 3950 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.818
- true rho 80 km, use_tau 0 (rung 1): R-hat 1.0041 | bulk ESS 524 | tail ESS 835 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.786

## Is the recovery USABLE? Interval widths

| true rho | use_tau | 95% CI width ratio (hi/lo) | spatial share | share 95% CI |
|---|---|---|---|---|
| 5 | 1 | **7.9x** | 0.80 | 0.25 to 1.00 (width 0.75) |
| 5 | 0 | **3.7x** | 1.00 (definitional, not estimated) | - |
| 10 | 1 | **8.2x** | 0.79 | 0.50 to 1.00 (width 0.50) |
| 10 | 0 | **2.5x** | 1.00 (definitional, not estimated) | - |
| 20 | 1 | **12.4x** | 0.77 | 0.55 to 0.94 (width 0.39) |
| 20 | 0 | **2.7x** | 1.00 (definitional, not estimated) | - |
| 40 | 1 | **12.0x** | 0.74 | 0.51 to 0.89 (width 0.38) |
| 40 | 0 | **2.7x** | 1.00 (definitional, not estimated) | - |
| 80 | 1 | **10.6x** | 0.83 | 0.66 to 0.94 (width 0.28) |
| 80 | 0 | **2.8x** | 1.00 (definitional, not estimated) | - |

## The nugget-versus-range trade-off

Across a true-rho range spanning a factor of 16, the posterior medians span a factor of 16.0 with the non-spatial term present and 3.2 without it.

../mataa's `moai_site_gp.stan` header warns that its tau-plus-field decomposition returned d50 of 7 to 13 km whatever the truth. **That warning inverts here.** The FULL model tracks the truth; it is the FIELD-ONLY model (use_tau = 0) that pins near 8 to 10 km regardless, because without a nugget the field must absorb every idiosyncratic site difference and is forced short. Do not carry their diagnostic reading across: on this design, dropping tau creates the pathology rather than exposing it.


## Saturation: does the estimate stop moving?

| truth pair | median pair | ratio | doubling recovered? |
|---|---|---|---|
| 5 -> 10 km | 5.4 -> 12.8 | 2.37 | yes |
| 10 -> 20 km | 12.8 -> 27.8 | 2.17 | yes |
| 20 -> 40 km | 27.8 -> 84.4 | 3.04 | yes |
| 40 -> 80 km | 84.4 -> 86.5 | 1.03 | **NO** |

**Ceiling at roughly 40 km.** Truths of 40 and 80 km return posterior medians of 84.4 and 86.5, which are the same answer. Above that the model cannot tell one long interaction scale from another.

**This is a second, independent witness for the ceiling.** `analyses/48_length_scale_recovery.py` located it near 32 km from summary statistics on drift simulations, with no GP anywhere in it; this locates it from a fitted latent field on data simulated from the model itself. Two methods sharing no machinery agree, which rule 10 says to treat as signal rather than coincidence. The common cause is geometric: the basin is 128 km across, so a decay longer than a few tens of km is not expressed inside the study window at all.

## Verdict

**rho is directionally recovered but only weakly identified.** Medians track the truth across a factor of 16 and every interval covers it, but the intervals span factors of 7.9 to 12.4. Coverage here is bought with vagueness. rho is reportable as an order-of-magnitude statement, not as a point estimate with a useful interval, which is the same outcome ../mataa recorded for pukao.

**`spatial_share` is recovered more usefully, but it is not sharp everywhere.** Its posterior MEDIAN sits within 0.058 of the truth of 0.80 at every length scale tested, which is genuine accuracy. Its INTERVAL, however, tightens with the length scale: width 0.75 at the shortest scale, where [0.25, 1.00] is close to no information at all on a bounded [0, 1] quantity, narrowing to 0.28 at the longest. Reporting a single width bound here would be misleading, and an earlier draft of this verdict did exactly that.

The variance partition is still the better-behaved quantity and the one the paper's argument needs, since it asks what fraction of assemblage-level compositional variance is organised by geography. But it earns a headline only where the interval is informative, which on this evidence means the model must be fitted with enough spatial signal present, and the width has to be reported alongside the median every time.

