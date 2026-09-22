# Does marginalising the nugget remove the ridge?

Produced by `analyses/02_spatial/03_geometry_repair.R` on 2026-09-21. Both models fitted to the SAME simulated data at known rho, under the resolvable-range inverse-gamma prior, so the only difference is the latent structure.

## Side by side

| true rho | model | rho median | 95% CI | width | spatial share | treedepth sat. | bulk ESS |
|---|---|---|---|---|---|---|---|
| 10 | additive | 9.8 | 6.0 to 18.7 | 3.1x | 0.91 | **15%** | 261 |
| 10 | marginal | 9.8 | 6.0 to 17.5 | 2.9x | 0.92 | **0%** | 898 |
| 20 | additive | 17.0 | 9.2 to 42.0 | 4.6x | 0.82 | **0%** | 474 |
| 20 | marginal | 17.0 | 9.3 to 42.7 | 4.6x | 0.82 | **0%** | 852 |
| 40 | additive | 16.3 | 7.3 to 52.6 | 7.2x | 0.65 | **0%** | 628 |
| 40 | marginal | 16.3 | 7.4 to 52.4 | 7.1x | 0.64 | **0%** | 610 |

## Diagnostics (rule 16, every fit)

- rho 10, additive (rung 2, NOT PASSING): R-hat 1.0083 | bulk ESS 261 | tail ESS 2809 | divergences 0/12000 (0.00%) | treedepth>=10 1851 | E-BFMI 0.775
- rho 10, marginal (rung 1): R-hat 1.0049 | bulk ESS 898 | tail ESS 2101 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.660
- rho 20, additive (rung 1): R-hat 1.0049 | bulk ESS 474 | tail ESS 1802 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.781
- rho 20, marginal (rung 1): R-hat 1.0079 | bulk ESS 852 | tail ESS 1315 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.680
- rho 40, additive (rung 1): R-hat 1.0061 | bulk ESS 628 | tail ESS 2135 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.779
- rho 40, marginal (rung 1): R-hat 1.0045 | bulk ESS 610 | tail ESS 1236 | divergences 0/6000 (0.00%) | treedepth>=10 1 | E-BFMI 0.690

## Movement

- Treedepth saturation: **5% additive -> 0% marginalised** (mean).
- rho interval width: 5.0x -> 4.9x (mean).
- bulk ESS: 454 -> 786 (mean).

**Partial.** The marginalised form is materially better but still saturating. Something beyond the u/f redundancy is contributing; the sigma-rho trade-off intrinsic to a GP on 29 points is the next candidate.

