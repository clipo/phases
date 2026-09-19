# Does marginalising the nugget remove the ridge?

Produced by `analyses/02_spatial/03_geometry_repair.R` on 2026-08-31. Both models fitted to the SAME simulated data at known rho, under the resolvable-range inverse-gamma prior, so the only difference is the latent structure.

## Side by side

| true rho | model | rho median | 95% CI | width | spatial share | treedepth sat. | bulk ESS |
|---|---|---|---|---|---|---|---|
| 10 | additive | 13.8 | 7.2 to 35.4 | 4.9x | 0.76 | **75%** | 399 |
| 10 | marginal | 13.6 | 7.1 to 35.9 | 5.1x | 0.77 | **0%** | 702 |
| 20 | additive | 21.3 | 10.6 to 60.5 | 5.7x | 0.79 | **0%** | 547 |
| 20 | marginal | 21.6 | 10.6 to 60.4 | 5.7x | 0.78 | **17%** | 1249 |
| 40 | additive | 39.4 | 15.9 to 136.0 | 8.6x | 0.69 | **0%** | 1306 |
| 40 | marginal | 38.4 | 15.4 to 135.4 | 8.8x | 0.68 | **0%** | 782 |

## Diagnostics (rule 16, every fit)

- rho 10, additive (rung 2, NOT PASSING): R-hat 1.0043 | bulk ESS 399 | tail ESS 2190 | divergences 0/12000 (0.00%) | treedepth>=10 9042 | E-BFMI 0.795
- rho 10, marginal (rung 1): R-hat 1.0075 | bulk ESS 702 | tail ESS 1557 | divergences 0/6000 (0.00%) | treedepth>=10 4 | E-BFMI 0.687
- rho 20, additive (rung 1): R-hat 1.0033 | bulk ESS 547 | tail ESS 1647 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.750
- rho 20, marginal (rung 2, NOT PASSING): R-hat 1.0019 | bulk ESS 1249 | tail ESS 2047 | divergences 0/12000 (0.00%) | treedepth>=10 2081 | E-BFMI 0.709
- rho 40, additive (rung 1): R-hat 1.0020 | bulk ESS 1306 | tail ESS 2128 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.799
- rho 40, marginal (rung 1): R-hat 1.0049 | bulk ESS 782 | tail ESS 1407 | divergences 0/6000 (0.00%) | treedepth>=10 8 | E-BFMI 0.657

## Movement

- Treedepth saturation: **25% additive -> 6% marginalised** (mean).
- rho interval width: 6.4x -> 6.5x (mean).
- bulk ESS: 751 -> 911 (mean).

**Partial.** The marginalised form is materially better but still saturating. Something beyond the u/f redundancy is contributing; the sigma-rho trade-off intrinsic to a GP on 29 points is the next candidate.

