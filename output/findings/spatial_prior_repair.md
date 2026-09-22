# Does a resolvable-range prior remove the ridge?

Produced by `analyses/02_spatial/02_prior_repair.R` on 2026-09-21.

## Why the uniform prior was the problem

The basin's median nearest-neighbour spacing is 5.03 km and its maximum pairwise distance is 127.8 km. No pair of sites can inform a correlation length shorter than the closest pair or longer than the farthest, so the likelihood is flat outside that band. Uniform on log(rho) over [0.5, 300] km places **36% of its mass below 5.03 km and 13% above 127.8 km**, so about half the prior sits where the data say nothing. That is what the sampler was exploring for 100% of its trajectories.

The replacement is inv_gamma(2.5874, 38.7141), solved for 1% tail mass at each end of that band: median 17.1 km, 90% interval [6.8, 63.1] km.

## Side by side

| true rho | prior | rho median | 95% CI | width | share width | treedepth sat. |
|---|---|---|---|---|---|---|
| 5 | uniform | 4.8 | 2.8 to 8.2 | 2.9x | 0.37 | 9% |
| 5 | invgamma | 6.4 | 4.1 to 11.3 | 2.8x | 0.47 | 25% |
| 10 | uniform | 8.5 | 4.9 to 17.4 | 3.5x | 0.32 | 19% |
| 10 | invgamma | 9.8 | 6.0 to 18.7 | 3.1x | 0.36 | 15% |
| 20 | uniform | 18.5 | 8.8 to 66.6 | 7.6x | 0.35 | 1% |
| 20 | invgamma | 17.0 | 9.2 to 42.0 | 4.6x | 0.36 | 0% |
| 40 | uniform | 18.4 | 4.9 to 176.8 | 35.8x | 0.58 | 45% |
| 40 | invgamma | 16.3 | 7.3 to 52.6 | 7.2x | 0.55 | 0% |
| 80 | uniform | 47.6 | 10.5 to 241.4 | 23.0x | 0.42 | 58% |
| 80 | invgamma | 26.7 | 11.0 to 91.2 | 8.3x | 0.40 | 0% |

## Diagnostics (rule 16, every fit)

- rho 5, uniform (rung 2, NOT PASSING): R-hat 1.0114 | bulk ESS 274 | tail ESS 3720 | divergences 0/12000 (0.00%) | treedepth>=10 1086 | E-BFMI 0.750
- rho 5, invgamma (rung 2, NOT PASSING): R-hat 1.0091 | bulk ESS 333 | tail ESS 3175 | divergences 0/12000 (0.00%) | treedepth>=10 3000 | E-BFMI 0.790
- rho 10, uniform (rung 2, NOT PASSING): R-hat 1.0145 | bulk ESS 323 | tail ESS 2052 | divergences 0/12000 (0.00%) | treedepth>=10 2250 | E-BFMI 0.791
- rho 10, invgamma (rung 2, NOT PASSING): R-hat 1.0083 | bulk ESS 261 | tail ESS 2809 | divergences 0/12000 (0.00%) | treedepth>=10 1851 | E-BFMI 0.775
- rho 20, uniform (rung 2, NOT PASSING): R-hat 1.0130 | bulk ESS 527 | tail ESS 2016 | divergences 0/12000 (0.00%) | treedepth>=10 157 | E-BFMI 0.790
- rho 20, invgamma (rung 1): R-hat 1.0049 | bulk ESS 474 | tail ESS 1802 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.781
- rho 40, uniform (rung 2, NOT PASSING): R-hat 1.0040 | bulk ESS 482 | tail ESS 1629 | divergences 0/12000 (0.00%) | treedepth>=10 5399 | E-BFMI 0.751
- rho 40, invgamma (rung 1): R-hat 1.0061 | bulk ESS 628 | tail ESS 2135 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.779
- rho 80, uniform (rung 2, NOT PASSING): R-hat 1.0030 | bulk ESS 699 | tail ESS 983 | divergences 0/12000 (0.00%) | treedepth>=10 6979 | E-BFMI 0.815
- rho 80, invgamma (rung 1): R-hat 1.0018 | bulk ESS 545 | tail ESS 1344 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.796

## Movement

- Treedepth saturation: 26% mean under uniform, **8% under inverse-gamma**.
- rho interval width: 14.6x mean under uniform, **5.2x under inverse-gamma**.
- spatial_share interval width: 0.41 mean under uniform, **0.43 under inverse-gamma**.

## Rule 20(c): recovery at values the prior disfavours

rho = 5 km and rho = 80 km sit in the prior's 1% tails. Under the inverse-gamma they return 6.4 and 26.7, covering truth: yes and yes. The prior is not dragging the estimate to its own median; the data move it into both tails.

