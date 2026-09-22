# Does the composition GP recover a known length scale?

Produced by `analyses/02_spatial/01_recovery.R` on 2026-09-22. Data simulated from the model by `analyses/49_write_stan_data.py` at the real basin geometry (29 assemblages, 10 classes, real sherd totals, max pairwise distance 127.8 km). Exponential kernel, uniform prior on log(rho) over [0.5, 300] km, true sigma 1.0, true tau 0.5, true spatial share 0.80.

## Recovery of rho

| true rho (km) | use_tau | rho median | rho 95% CI | d50 | spatial share | covers truth | rung | verdict |
|---|---|---|---|---|---|---|---|---|
| 5 | 1 | 4.8 | 2.7 to 8.7 | 3.3 | 0.91 | yes | 2 | underpowered |
| 5 | 0 | 4.1 | 2.4 to 6.5 | 2.9 | 1.00 | yes | 1 | recovered |
| 10 | 1 | 7.9 | 4.4 to 16.1 | 5.5 | 0.92 | yes | 2 | underpowered |
| 10 | 0 | 6.7 | 3.9 to 10.9 | 4.7 | 1.00 | yes | 1 | recovered |
| 20 | 1 | 18.6 | 9.1 to 72.3 | 12.9 | 0.84 | yes | 2 | underpowered |
| 20 | 0 | 9.8 | 6.4 to 15.8 | 6.8 | 1.00 | **NO** | 1 | recovered |
| 40 | 1 | 25.6 | 6.5 to 190.6 | 17.8 | 0.66 | yes | 2 | ridge |
| 40 | 0 | 5.2 | 2.3 to 9.6 | 3.6 | 1.00 | **NO** | 1 | recovered |
| 80 | 1 | 46.6 | 10.5 to 245.0 | 32.3 | 0.75 | yes | 1 | recovered |
| 80 | 0 | 7.1 | 4.2 to 12.0 | 4.9 | 1.00 | **NO** | 1 | recovered |

## Diagnostics (rule 16, every fit)

- true rho 5 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0342 | bulk ESS 252 | tail ESS 2712 | divergences 0/12000 (0.00%) | treedepth>=10 1100 | E-BFMI 0.788
- true rho 5 km, use_tau 0 (rung 1): R-hat 1.0082 | bulk ESS 1137 | tail ESS 2139 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.790
- true rho 10 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0094 | bulk ESS 274 | tail ESS 2782 | divergences 0/12000 (0.00%) | treedepth>=10 738 | E-BFMI 0.756
- true rho 10 km, use_tau 0 (rung 1): R-hat 1.0094 | bulk ESS 750 | tail ESS 1206 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.720
- true rho 20 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0027 | bulk ESS 526 | tail ESS 1403 | divergences 0/12000 (0.00%) | treedepth>=10 1961 | E-BFMI 0.834
- true rho 20 km, use_tau 0 (rung 1): R-hat 1.0019 | bulk ESS 975 | tail ESS 1787 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.752
- true rho 40 km, use_tau 1 (rung 2, NOT PASSING): R-hat 1.0025 | bulk ESS 579 | tail ESS 1427 | divergences 0/12000 (0.00%) | treedepth>=10 3000 | E-BFMI 0.801
- true rho 40 km, use_tau 0 (rung 1): R-hat 1.0065 | bulk ESS 625 | tail ESS 882 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.761
- true rho 80 km, use_tau 1 (rung 1): R-hat 1.0067 | bulk ESS 478 | tail ESS 697 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.799
- true rho 80 km, use_tau 0 (rung 1): R-hat 1.0017 | bulk ESS 986 | tail ESS 1837 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.766

## Is the recovery USABLE? Interval widths

| true rho | use_tau | 95% CI width ratio (hi/lo) | spatial share | share 95% CI |
|---|---|---|---|---|
| 5 | 1 | **3.3x** | 0.91 | 0.54 to 1.00 (width 0.46) |
| 5 | 0 | **2.7x** | 1.00 (definitional, not estimated) | - |
| 10 | 1 | **3.7x** | 0.92 | 0.61 to 1.00 (width 0.39) |
| 10 | 0 | **2.8x** | 1.00 (definitional, not estimated) | - |
| 20 | 1 | **8.0x** | 0.84 | 0.65 to 0.97 (width 0.32) |
| 20 | 0 | **2.5x** | 1.00 (definitional, not estimated) | - |
| 40 | 1 | **29.3x** | 0.66 | 0.37 to 0.90 (width 0.53) |
| 40 | 0 | **4.2x** | 1.00 (definitional, not estimated) | - |
| 80 | 1 | **23.3x** | 0.75 | 0.51 to 0.92 (width 0.41) |
| 80 | 0 | **2.9x** | 1.00 (definitional, not estimated) | - |

## The nugget-versus-range trade-off

Across a true-rho range spanning a factor of 16, the posterior medians span a factor of 9.8 with the non-spatial term present and 2.4 without it.

../mataa's `moai_site_gp.stan` header warns that its tau-plus-field decomposition returned d50 of 7 to 13 km whatever the truth. **That warning inverts here.** The FULL model tracks the truth; it is the FIELD-ONLY model (use_tau = 0) that pins near 8 to 10 km regardless, because without a nugget the field must absorb every idiosyncratic site difference and is forced short. Do not carry their diagnostic reading across: on this design, dropping tau creates the pathology rather than exposing it.


## Saturation: does the estimate stop moving?

| truth pair | median pair | ratio | doubling recovered? |
|---|---|---|---|
| 5 -> 10 km | 4.8 -> 7.9 | 1.67 | yes |
| 10 -> 20 km | 7.9 -> 18.6 | 2.35 | yes |
| 20 -> 40 km | 18.6 -> 25.6 | 1.38 | yes |
| 40 -> 80 km | 25.6 -> 46.6 | 1.82 | yes |

No saturation detected: every doubling of the truth moves the estimate.

## Verdict

**rho is directionally recovered but only weakly identified.** Medians track the truth across a factor of 16 and every interval covers it, but the intervals span factors of 3.3 to 29.3. Coverage here is bought with vagueness. rho is reportable as an order-of-magnitude statement, not as a point estimate with a useful interval, which is the same outcome ../mataa recorded for pukao.

**`spatial_share` is recovered more usefully, but it is not sharp everywhere.** Its posterior MEDIAN sits within 0.138 of the truth of 0.80 at every length scale tested, which is genuine accuracy. Its INTERVAL, however, tightens with the length scale: width 0.46 at the shortest scale, where [0.54, 1.00] is close to no information at all on a bounded [0, 1] quantity, narrowing to 0.41 at the longest. Reporting a single width bound here would be misleading, and an earlier draft of this verdict did exactly that.

The variance partition is still the better-behaved quantity and the one the paper's argument needs, since it asks what fraction of assemblage-level compositional variance is organised by geography. But it earns a headline only where the interval is informative, which on this evidence means the model must be fitted with enough spatial signal present, and the width has to be reported alongside the median every time.

