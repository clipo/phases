# Does marginalising the nugget remove the ridge?

Produced by `analyses/02_spatial/03_geometry_repair.R` on 2026-09-23. Both models fitted to the SAME simulated data at known rho, under the resolvable-range inverse-gamma prior, so the only difference is the latent structure.

## Side by side

| true rho | model | rho median | 95% CI | width | spatial share | treedepth sat. | bulk ESS |
|---|---|---|---|---|---|---|---|
| 10 | additive | 9.5 | 5.8 to 18.4 | 3.2x | 0.89 | **50%** | 357 |
| 10 | marginal | 9.5 | 5.8 to 18.2 | 3.2x | 0.90 | **0%** | 779 |
| 20 | additive | 16.7 | 9.1 to 39.1 | 4.3x | 0.84 | **0%** | 585 |
| 20 | marginal | 17.0 | 9.4 to 37.1 | 3.9x | 0.84 | **0%** | 1434 |
| 40 | additive | 19.2 | 8.3 to 61.1 | 7.3x | 0.65 | **1%** | 1064 |
| 40 | marginal | 19.6 | 8.4 to 66.9 | 7.9x | 0.65 | **0%** | 827 |

## Diagnostics (rule 16, every fit)

- rho 10, additive (rung 2, NOT PASSING): R-hat 1.0112 | bulk ESS 357 | tail ESS 3210 | divergences 0/12000 (0.00%) | treedepth>=10 5991 | E-BFMI 0.801
- rho 10, marginal (rung 1): R-hat 1.0091 | bulk ESS 779 | tail ESS 1685 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.682
- rho 20, additive (rung 2): R-hat 1.0027 | bulk ESS 585 | tail ESS 3290 | divergences 0/12000 (0.00%) | treedepth>=10 23 | E-BFMI 0.762
- rho 20, marginal (rung 1): R-hat 1.0020 | bulk ESS 1434 | tail ESS 2718 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.639
- rho 40, additive (rung 2): R-hat 1.0013 | bulk ESS 1064 | tail ESS 4005 | divergences 0/12000 (0.00%) | treedepth>=10 149 | E-BFMI 0.786
- rho 40, marginal (rung 1): R-hat 1.0027 | bulk ESS 827 | tail ESS 1702 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.696

## Movement

- Treedepth saturation: **17% additive -> 0% marginalised** (mean).
- rho interval width: 4.9x -> 5.0x (mean).
- bulk ESS: 669 -> 1013 (mean).

**Partial.** The marginalised form is materially better but still saturating. Something beyond the u/f redundancy is contributing; the sigma-rho trade-off intrinsic to a GP on 29 points is the next candidate.

