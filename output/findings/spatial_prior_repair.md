# Does a resolvable-range prior remove the ridge?

Produced by `analyses/02_spatial/02_prior_repair.R` on 2026-10-01.

## Why the uniform prior was the problem

The basin's median nearest-neighbour spacing is 5.03 km and its maximum pairwise distance is 127.8 km. No pair of sites can inform a correlation length shorter than the closest pair or longer than the farthest, so the likelihood is flat outside that band. Uniform on log(rho) over [0.5, 300] km places **36% of its mass below 5.03 km and 13% above 127.8 km**, so about half the prior sits where the data say nothing. That is what the sampler was exploring for 100% of its trajectories.

The replacement is inv_gamma(2.5874, 38.7141), solved for 1% tail mass at each end of that band: median 17.1 km, 90% interval [6.8, 63.1] km.

## Side by side

| true rho | prior | rho median | 95% CI | width | share width | treedepth sat. |
|---|---|---|---|---|---|---|
| 5 | uniform | 5.2 | 3.1 to 9.6 | 3.1x | 0.40 | 5% |
| 5 | invgamma | 6.7 | 4.3 to 11.7 | 2.7x | 0.46 | 4% |
| 10 | uniform | 7.9 | 4.3 to 18.7 | 4.4x | 0.41 | 0% |
| 10 | invgamma | 9.5 | 5.6 to 19.7 | 3.5x | 0.41 | 5% |
| 20 | uniform | 18.8 | 9.0 to 72.3 | 8.0x | 0.32 | 25% |
| 20 | invgamma | 17.1 | 9.4 to 41.3 | 4.4x | 0.33 | 0% |
| 40 | uniform | 16.6 | 5.9 to 135.4 | 22.8x | 0.54 | 1% |
| 40 | invgamma | 15.0 | 7.3 to 48.3 | 6.6x | 0.52 | 0% |
| 80 | uniform | 30.0 | 6.7 to 212.4 | 31.7x | 0.49 | 49% |
| 80 | invgamma | 18.3 | 7.3 to 73.4 | 10.0x | 0.51 | 0% |

## Diagnostics (rule 16, every fit)

- rho 5, uniform (rung 2, NOT PASSING): R-hat 1.0039 | bulk ESS 327 | tail ESS 2623 | divergences 0/12000 (0.00%) | treedepth>=10 624 | E-BFMI 0.785
- rho 5, invgamma (rung 2, NOT PASSING): R-hat 1.0175 | bulk ESS 321 | tail ESS 3727 | divergences 0/12000 (0.00%) | treedepth>=10 432 | E-BFMI 0.786
- rho 10, uniform (rung 2, NOT PASSING): R-hat 1.0213 | bulk ESS 295 | tail ESS 1764 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.763
- rho 10, invgamma (rung 2, NOT PASSING): R-hat 1.0110 | bulk ESS 362 | tail ESS 3452 | divergences 0/12000 (0.00%) | treedepth>=10 626 | E-BFMI 0.759
- rho 20, uniform (rung 2, NOT PASSING): R-hat 1.0065 | bulk ESS 903 | tail ESS 2486 | divergences 0/12000 (0.00%) | treedepth>=10 2961 | E-BFMI 0.808
- rho 20, invgamma (rung 2): R-hat 1.0055 | bulk ESS 646 | tail ESS 3642 | divergences 0/12000 (0.00%) | treedepth>=10 3 | E-BFMI 0.796
- rho 40, uniform (rung 2, NOT PASSING): R-hat 1.0105 | bulk ESS 547 | tail ESS 1510 | divergences 0/12000 (0.00%) | treedepth>=10 80 | E-BFMI 0.797
- rho 40, invgamma (rung 2): R-hat 1.0046 | bulk ESS 945 | tail ESS 2990 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.792
- rho 80, uniform (rung 2, NOT PASSING): R-hat 1.0043 | bulk ESS 406 | tail ESS 821 | divergences 0/12000 (0.00%) | treedepth>=10 5908 | E-BFMI 0.797
- rho 80, invgamma (rung 2): R-hat 1.0033 | bulk ESS 549 | tail ESS 1668 | divergences 0/12000 (0.00%) | treedepth>=10 0 | E-BFMI 0.809

## Movement

- Treedepth saturation: 16% mean under uniform, **2% under inverse-gamma**.
- rho interval width: 14.0x mean under uniform, **5.4x under inverse-gamma**.
- spatial_share interval width: 0.43 mean under uniform, **0.44 under inverse-gamma**.

## Rule 20(c): recovery at values the prior disfavours

rho = 5 km and rho = 80 km sit in the prior's 1% tails. Under the inverse-gamma they return 6.7 and 18.3, covering truth: yes and **NO**. **At least one disfavoured truth is NOT recovered. The prior is doing work the data should be doing, and this prior is not reportable under rule 20(c).**

