# Does a resolvable-range prior remove the ridge?

Produced by `analyses/02_spatial/02_prior_repair.R` on 2026-09-20.

## Why the uniform prior was the problem

The basin's median nearest-neighbour spacing is 5.03 km and its maximum pairwise distance is 127.8 km. No pair of sites can inform a correlation length shorter than the closest pair or longer than the farthest, so the likelihood is flat outside that band. Uniform on log(rho) over [0.5, 300] km places **36% of its mass below 5.03 km and 13% above 127.8 km**, so about half the prior sits where the data say nothing. That is what the sampler was exploring for 100% of its trajectories.

The replacement is inv_gamma(2.5874, 38.7141), solved for 1% tail mass at each end of that band: median 17.1 km, 90% interval [6.8, 63.1] km.

## Side by side

| true rho | prior | rho median | 95% CI | width | share width | treedepth sat. |
|---|---|---|---|---|---|---|
| 5 | uniform | 11.5 | 2.6 to 221.9 | 84.4x | 0.67 | 100% |
| 5 | invgamma | 13.9 | 5.2 to 62.1 | 11.9x | 0.58 | 61% |
| 10 | uniform | 8.6 | 5.3 to 15.1 | 2.9x | 0.30 | 0% |
| 10 | invgamma | 9.6 | 6.1 to 16.8 | 2.7x | 0.30 | 0% |
| 20 | uniform | 17.7 | 10.4 to 38.7 | 3.7x | 0.25 | 25% |
| 20 | invgamma | 17.2 | 10.5 to 33.0 | 3.1x | 0.25 | 59% |
| 40 | uniform | 95.2 | 28.9 to 269.1 | 9.3x | 0.26 | 100% |
| 40 | invgamma | 47.4 | 20.7 to 148.7 | 7.2x | 0.26 | 0% |
| 80 | uniform | 75.2 | 26.0 to 251.4 | 9.7x | 0.33 | 100% |
| 80 | invgamma | 43.4 | 20.7 to 121.5 | 5.9x | 0.31 | 0% |

## Diagnostics (rule 16, every fit)

- rho 5, uniform (rung 2, NOT PASSING): R-hat 1.0187 | bulk ESS 189 | tail ESS 1205 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.777
- rho 5, invgamma (rung 2, NOT PASSING): R-hat 1.0059 | bulk ESS 505 | tail ESS 1882 | divergences 0/12000 (0.00%) | treedepth>=10 7305 | E-BFMI 0.825
- rho 10, uniform (rung 1): R-hat 1.0047 | bulk ESS 588 | tail ESS 2013 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.769
- rho 10, invgamma (rung 1): R-hat 1.0076 | bulk ESS 422 | tail ESS 1926 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.728
- rho 20, uniform (rung 2, NOT PASSING): R-hat 1.0024 | bulk ESS 511 | tail ESS 1931 | divergences 0/12000 (0.00%) | treedepth>=10 3009 | E-BFMI 0.818
- rho 20, invgamma (rung 2, NOT PASSING): R-hat 1.0084 | bulk ESS 686 | tail ESS 3006 | divergences 0/12000 (0.00%) | treedepth>=10 7090 | E-BFMI 0.806
- rho 40, uniform (rung 2, NOT PASSING): R-hat 1.0012 | bulk ESS 2560 | tail ESS 4597 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.785
- rho 40, invgamma (rung 1): R-hat 1.0030 | bulk ESS 769 | tail ESS 1253 | divergences 0/6000 (0.00%) | treedepth>=10 3 | E-BFMI 0.784
- rho 80, uniform (rung 2, NOT PASSING): R-hat 1.0038 | bulk ESS 2429 | tail ESS 4017 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.822
- rho 80, invgamma (rung 1): R-hat 1.0018 | bulk ESS 1661 | tail ESS 2550 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.732

## Movement

- Treedepth saturation: 65% mean under uniform, **24% under inverse-gamma**.
- rho interval width: 22.0x mean under uniform, **6.2x under inverse-gamma**.
- spatial_share interval width: 0.36 mean under uniform, **0.34 under inverse-gamma**.

## Rule 20(c): recovery at values the prior disfavours

rho = 5 km and rho = 80 km sit in the prior's 1% tails. Under the inverse-gamma they return 13.9 and 43.4, covering truth: **NO** and yes. **At least one disfavoured truth is NOT recovered. The prior is doing work the data should be doing, and this prior is not reportable under rule 20(c).**

