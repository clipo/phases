# Does marginalising the nugget remove the ridge?

Produced by `analyses/02_spatial/03_geometry_repair.R` on 2026-09-22. Both models fitted to the SAME simulated data at known rho, under the resolvable-range inverse-gamma prior, so the only difference is the latent structure.

## Side by side

| true rho | model | rho median | 95% CI | width | spatial share | treedepth sat. | bulk ESS |
|---|---|---|---|---|---|---|---|
| 10 | additive | 9.4 | 5.8 to 17.7 | 3.1x | 0.90 | **1%** | 262 |
| 10 | marginal | 9.4 | 5.7 to 18.0 | 3.1x | 0.90 | **0%** | 1173 |
| 20 | additive | 17.7 | 9.8 to 43.2 | 4.4x | 0.83 | **0%** | 916 |
| 20 | marginal | 17.4 | 9.5 to 41.3 | 4.4x | 0.83 | **0%** | 922 |
| 40 | additive | 18.7 | 7.8 to 61.8 | 7.9x | 0.63 | **0%** | 466 |
| 40 | marginal | 18.7 | 8.0 to 63.6 | 8.0x | 0.63 | **0%** | 835 |

## Diagnostics (rule 16, every fit)

- rho 10, additive (rung 2, NOT PASSING): R-hat 1.0198 | bulk ESS 262 | tail ESS 3592 | divergences 0/12000 (0.00%) | treedepth>=10 137 | E-BFMI 0.798
- rho 10, marginal (rung 1): R-hat 1.0034 | bulk ESS 1173 | tail ESS 2218 | divergences 0/6000 (0.00%) | treedepth>=10 1 | E-BFMI 0.622
- rho 20, additive (rung 2): R-hat 1.0064 | bulk ESS 916 | tail ESS 3498 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.818
- rho 20, marginal (rung 1): R-hat 1.0018 | bulk ESS 922 | tail ESS 1431 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.664
- rho 40, additive (rung 1): R-hat 1.0079 | bulk ESS 466 | tail ESS 990 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.807
- rho 40, marginal (rung 1): R-hat 1.0048 | bulk ESS 835 | tail ESS 1611 | divergences 0/6000 (0.00%) | treedepth>=10 1 | E-BFMI 0.643

## Movement

- Treedepth saturation: **0% additive -> 0% marginalised** (mean).
- rho interval width: 5.1x -> 5.2x (mean).
- bulk ESS: 548 -> 977 (mean).

**Partial.** The marginalised form is materially better but still saturating. Something beyond the u/f redundancy is contributing; the sigma-rho trade-off intrinsic to a GP on 29 points is the next candidate.

