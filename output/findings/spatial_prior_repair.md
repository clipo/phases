# Does a resolvable-range prior remove the ridge?

Produced by `analyses/02_spatial/02_prior_repair.R` on 2026-09-22.

## Why the uniform prior was the problem

The basin's median nearest-neighbour spacing is 5.03 km and its maximum pairwise distance is 127.8 km. No pair of sites can inform a correlation length shorter than the closest pair or longer than the farthest, so the likelihood is flat outside that band. Uniform on log(rho) over [0.5, 300] km places **36% of its mass below 5.03 km and 13% above 127.8 km**, so about half the prior sits where the data say nothing. That is what the sampler was exploring for 100% of its trajectories.

The replacement is inv_gamma(2.5874, 38.7141), solved for 1% tail mass at each end of that band: median 17.1 km, 90% interval [6.8, 63.1] km.

## Side by side

| true rho | prior | rho median | 95% CI | width | share width | treedepth sat. |
|---|---|---|---|---|---|---|
| 5 | uniform | 4.8 | 2.7 to 8.7 | 3.3x | 0.46 | 9% |
| 5 | invgamma | 6.3 | 4.1 to 11.3 | 2.8x | 0.51 | 19% |
| 10 | uniform | 7.9 | 4.4 to 16.1 | 3.7x | 0.39 | 6% |
| 10 | invgamma | 9.4 | 5.8 to 17.7 | 3.1x | 0.39 | 1% |
| 20 | uniform | 18.6 | 9.1 to 72.3 | 8.0x | 0.32 | 16% |
| 20 | invgamma | 17.7 | 9.8 to 43.2 | 4.4x | 0.32 | 0% |
| 40 | uniform | 25.6 | 6.5 to 190.6 | 29.3x | 0.53 | 25% |
| 40 | invgamma | 18.7 | 7.8 to 61.8 | 7.9x | 0.52 | 0% |
| 80 | uniform | 46.6 | 10.5 to 245.0 | 23.3x | 0.41 | 0% |
| 80 | invgamma | 25.1 | 10.1 to 83.1 | 8.2x | 0.42 | 0% |

## Diagnostics (rule 16, every fit)

- rho 5, uniform (rung 2, NOT PASSING): R-hat 1.0342 | bulk ESS 252 | tail ESS 2712 | divergences 0/12000 (0.00%) | treedepth>=10 1100 | E-BFMI 0.788
- rho 5, invgamma (rung 2, NOT PASSING): R-hat 1.0116 | bulk ESS 310 | tail ESS 3047 | divergences 0/12000 (0.00%) | treedepth>=10 2334 | E-BFMI 0.765
- rho 10, uniform (rung 2, NOT PASSING): R-hat 1.0094 | bulk ESS 274 | tail ESS 2782 | divergences 0/12000 (0.00%) | treedepth>=10 738 | E-BFMI 0.756
- rho 10, invgamma (rung 2, NOT PASSING): R-hat 1.0198 | bulk ESS 262 | tail ESS 3592 | divergences 0/12000 (0.00%) | treedepth>=10 137 | E-BFMI 0.798
- rho 20, uniform (rung 2, NOT PASSING): R-hat 1.0027 | bulk ESS 526 | tail ESS 1403 | divergences 0/12000 (0.00%) | treedepth>=10 1961 | E-BFMI 0.834
- rho 20, invgamma (rung 2): R-hat 1.0064 | bulk ESS 916 | tail ESS 3498 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.818
- rho 40, uniform (rung 2, NOT PASSING): R-hat 1.0025 | bulk ESS 579 | tail ESS 1427 | divergences 0/12000 (0.00%) | treedepth>=10 3000 | E-BFMI 0.801
- rho 40, invgamma (rung 1): R-hat 1.0079 | bulk ESS 466 | tail ESS 990 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.807
- rho 80, uniform (rung 1): R-hat 1.0067 | bulk ESS 478 | tail ESS 697 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.799
- rho 80, invgamma (rung 1): R-hat 1.0091 | bulk ESS 427 | tail ESS 963 | divergences 0/6000 (0.00%) | treedepth>=10 0 | E-BFMI 0.769

## Movement

- Treedepth saturation: 11% mean under uniform, **4% under inverse-gamma**.
- rho interval width: 13.5x mean under uniform, **5.3x under inverse-gamma**.
- spatial_share interval width: 0.42 mean under uniform, **0.43 under inverse-gamma**.

## Rule 20(c): recovery at values the prior disfavours

rho = 5 km and rho = 80 km sit in the prior's 1% tails. Under the inverse-gamma they return 6.3 and 25.1, covering truth: yes and yes. The prior is not dragging the estimate to its own median; the data move it into both tails.

