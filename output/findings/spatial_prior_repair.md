# Does a resolvable-range prior remove the ridge?

Produced by `analyses/02_spatial/02_prior_repair.R` on 2026-09-23.

## Why the uniform prior was the problem

The basin's median nearest-neighbour spacing is 5.03 km and its maximum pairwise distance is 127.8 km. No pair of sites can inform a correlation length shorter than the closest pair or longer than the farthest, so the likelihood is flat outside that band. Uniform on log(rho) over [0.5, 300] km places **36% of its mass below 5.03 km and 13% above 127.8 km**, so about half the prior sits where the data say nothing. That is what the sampler was exploring for 100% of its trajectories.

The replacement is inv_gamma(2.5874, 38.7141), solved for 1% tail mass at each end of that band: median 17.1 km, 90% interval [6.8, 63.1] km.

## Side by side

| true rho | prior | rho median | 95% CI | width | share width | treedepth sat. |
|---|---|---|---|---|---|---|
| 5 | uniform | 4.6 | 2.5 to 8.3 | 3.3x | 0.45 | 11% |
| 5 | invgamma | 6.3 | 4.0 to 12.1 | 3.0x | 0.58 | 13% |
| 10 | uniform | 7.9 | 4.3 to 17.1 | 3.9x | 0.38 | 0% |
| 10 | invgamma | 9.5 | 5.8 to 18.4 | 3.2x | 0.42 | 50% |
| 20 | uniform | 18.5 | 8.9 to 69.3 | 7.8x | 0.33 | 23% |
| 20 | invgamma | 16.7 | 9.1 to 39.1 | 4.3x | 0.34 | 0% |
| 40 | uniform | 26.6 | 7.0 to 188.2 | 26.7x | 0.52 | 25% |
| 40 | invgamma | 19.2 | 8.3 to 61.1 | 7.3x | 0.51 | 1% |
| 80 | uniform | 49.6 | 9.9 to 247.6 | 25.0x | 0.39 | 75% |
| 80 | invgamma | 24.3 | 9.5 to 82.2 | 8.7x | 0.40 | 0% |

## Diagnostics (rule 16, every fit)

- rho 5, uniform (rung 2, NOT PASSING): R-hat 1.0191 | bulk ESS 157 | tail ESS 2056 | divergences 0/12000 (0.00%) | treedepth>=10 1375 | E-BFMI 0.762
- rho 5, invgamma (rung 2, NOT PASSING): R-hat 1.0122 | bulk ESS 380 | tail ESS 1766 | divergences 0/12000 (0.00%) | treedepth>=10 1560 | E-BFMI 0.770
- rho 10, uniform (rung 2, NOT PASSING): R-hat 1.0113 | bulk ESS 337 | tail ESS 2171 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.752
- rho 10, invgamma (rung 2, NOT PASSING): R-hat 1.0112 | bulk ESS 357 | tail ESS 3210 | divergences 0/12000 (0.00%) | treedepth>=10 5991 | E-BFMI 0.801
- rho 20, uniform (rung 2, NOT PASSING): R-hat 1.0023 | bulk ESS 753 | tail ESS 2334 | divergences 0/12000 (0.00%) | treedepth>=10 2729 | E-BFMI 0.793
- rho 20, invgamma (rung 2): R-hat 1.0027 | bulk ESS 585 | tail ESS 3290 | divergences 0/12000 (0.00%) | treedepth>=10 23 | E-BFMI 0.762
- rho 40, uniform (rung 2, NOT PASSING): R-hat 1.0054 | bulk ESS 706 | tail ESS 1332 | divergences 0/12000 (0.00%) | treedepth>=10 3000 | E-BFMI 0.761
- rho 40, invgamma (rung 2): R-hat 1.0013 | bulk ESS 1064 | tail ESS 4005 | divergences 0/12000 (0.00%) | treedepth>=10 149 | E-BFMI 0.786
- rho 80, uniform (rung 2, NOT PASSING): R-hat 1.0033 | bulk ESS 502 | tail ESS 684 | divergences 0/12000 (0.00%) | treedepth>=10 9001 | E-BFMI 0.777
- rho 80, invgamma (rung 1): R-hat 1.0032 | bulk ESS 408 | tail ESS 877 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.749

## Movement

- Treedepth saturation: 27% mean under uniform, **13% under inverse-gamma**.
- rho interval width: 13.3x mean under uniform, **5.3x under inverse-gamma**.
- spatial_share interval width: 0.41 mean under uniform, **0.45 under inverse-gamma**.

## Rule 20(c): recovery at values the prior disfavours

rho = 5 km and rho = 80 km sit in the prior's 1% tails. Under the inverse-gamma they return 6.3 and 24.3, covering truth: yes and yes. The prior is not dragging the estimate to its own median; the data move it into both tails.

