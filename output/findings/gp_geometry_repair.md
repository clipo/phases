# Does marginalising the nugget remove the ridge?

Produced by `analyses/02_spatial/03_geometry_repair.R` on 2026-09-20. Both models fitted to the SAME simulated data at known rho, under the resolvable-range inverse-gamma prior, so the only difference is the latent structure.

## Side by side

| true rho | model | rho median | 95% CI | width | spatial share | treedepth sat. | bulk ESS |
|---|---|---|---|---|---|---|---|
| 10 | additive | 9.6 | 6.1 to 16.8 | 2.7x | 0.78 | **0%** | 422 |
| 10 | marginal | 9.6 | 6.2 to 16.8 | 2.7x | 0.78 | **0%** | 597 |
| 20 | additive | 17.2 | 10.5 to 33.0 | 3.1x | 0.88 | **59%** | 686 |
| 20 | marginal | 16.9 | 10.3 to 32.8 | 3.2x | 0.88 | **0%** | 538 |
| 40 | additive | 47.4 | 20.7 to 148.7 | 7.2x | 0.80 | **0%** | 769 |
| 40 | marginal | 47.9 | 21.4 to 159.8 | 7.5x | 0.79 | **0%** | 1214 |

## Diagnostics (rule 16, every fit)

- rho 10, additive (rung 1): R-hat 1.0076 | bulk ESS 422 | tail ESS 1926 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.728
- rho 10, marginal (rung 1): R-hat 1.0058 | bulk ESS 597 | tail ESS 1228 | divergences 0/6000 (0.00%) | treedepth>=10 2 | E-BFMI 0.665
- rho 20, additive (rung 2, NOT PASSING): R-hat 1.0084 | bulk ESS 686 | tail ESS 3006 | divergences 0/12000 (0.00%) | treedepth>=10 7090 | E-BFMI 0.806
- rho 20, marginal (rung 1): R-hat 1.0039 | bulk ESS 538 | tail ESS 953 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.630
- rho 40, additive (rung 1): R-hat 1.0030 | bulk ESS 769 | tail ESS 1253 | divergences 0/6000 (0.00%) | treedepth>=10 3 | E-BFMI 0.784
- rho 40, marginal (rung 1): R-hat 1.0018 | bulk ESS 1214 | tail ESS 2042 | divergences 0/6000 (0.00%) | treedepth>=10 1 | E-BFMI 0.633

## Movement

- Treedepth saturation: **20% additive -> 0% marginalised** (mean).
- rho interval width: 4.3x -> 4.5x (mean).
- bulk ESS: 626 -> 783 (mean).

**Partial.** The marginalised form is materially better but still saturating. Something beyond the u/f redundancy is contributing; the sigma-rho trade-off intrinsic to a GP on 29 points is the next candidate.

