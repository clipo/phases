# Does marginalising the nugget remove the ridge?

Produced by `analyses/02_spatial/03_geometry_repair.R` on 2026-10-01. Both models fitted to the SAME simulated data at known rho, under the resolvable-range inverse-gamma prior, so the only difference is the latent structure.

## Side by side

| true rho | model | rho median | 95% CI | width | spatial share | treedepth sat. | bulk ESS |
|---|---|---|---|---|---|---|---|
| 10 | additive | 9.5 | 5.6 to 19.7 | 3.5x | 0.88 | **5%** | 362 |
| 10 | marginal | 9.4 | 5.6 to 18.8 | 3.4x | 0.88 | **0%** | 821 |
| 20 | additive | 17.1 | 9.4 to 41.3 | 4.4x | 0.84 | **0%** | 646 |
| 20 | marginal | 17.2 | 9.3 to 42.5 | 4.6x | 0.84 | **0%** | 943 |
| 40 | additive | 15.0 | 7.3 to 48.3 | 6.6x | 0.66 | **0%** | 945 |
| 40 | marginal | 15.3 | 7.1 to 46.2 | 6.5x | 0.66 | **0%** | 566 |

## Diagnostics (rule 16, every fit)

- rho 10, additive (rung 2, NOT PASSING): R-hat 1.0110 | bulk ESS 362 | tail ESS 3452 | divergences 0/12000 (0.00%) | treedepth>=10 626 | E-BFMI 0.759
- rho 10, marginal (rung 1): R-hat 1.0011 | bulk ESS 821 | tail ESS 1512 | divergences 0/6000 (0.00%) | treedepth>=10 1 | E-BFMI 0.656
- rho 20, additive (rung 2): R-hat 1.0055 | bulk ESS 646 | tail ESS 3642 | divergences 0/12000 (0.00%) | treedepth>=10 3 | E-BFMI 0.796
- rho 20, marginal (rung 1): R-hat 1.0045 | bulk ESS 943 | tail ESS 1453 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.671
- rho 40, additive (rung 2): R-hat 1.0046 | bulk ESS 945 | tail ESS 2990 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.792
- rho 40, marginal (rung 1): R-hat 1.0070 | bulk ESS 566 | tail ESS 953 | divergences 0/6000 (0.00%) | treedepth>=10 1 | E-BFMI 0.713

## Movement

- Treedepth saturation: **2% additive -> 0% marginalised** (mean).
- rho interval width: 4.8x -> 4.8x (mean).
- bulk ESS: 651 -> 777 (mean).

**Partial.** The marginalised form is materially better but still saturating. Something beyond the u/f redundancy is contributing; the sigma-rho trade-off intrinsic to a GP on 29 points is the next candidate.

