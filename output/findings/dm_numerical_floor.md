# Where does the Dirichlet-multinomial density break?

Produced by `analyses/58_dm_numerical_floor.py`. Closes F8.

Test case: a realistic bin, 10 classes and n = 200 sherds, with symmetric pi. The reference is the exact sum-of-logs form, which calls no gamma function and therefore shares no arithmetic with the implementation it audits (rule 3).

| a0 | implied F | exact log pmf | PyMC log pmf | relative error |
|---|---|---|---|---|
| 1e+00 | 5.000e-01 | -49.245325 | -49.245325 | **1.30e-15** |
| 1e+02 | 9.901e-03 | -27.379506 | -27.379506 | **2.28e-14** |
| 1e+04 | 9.999e-05 | -26.230273 | -26.230273 | **7.26e-13** |
| 1e+06 | 1.000e-06 | -26.252702 | -26.252702 | **1.13e-13** |
| 1e+09 | 1.000e-09 | -26.252942 | -26.252942 | **1.81e-08** |
| 1e+12 | 1.000e-12 | -26.252942 | -26.258789 | **2.23e-04** |
| 1e+15 | 1.000e-15 | -26.252942 | -23.500000 | **1.05e-01** |
| 1e+17 | 1.000e-17 | -26.252942 | 192.000000 | **8.31e+00** |
| 1e+18 | 1.000e-18 | -26.252942 | -2048.000000 | **7.70e+01** |

## Reading

PyMC's density agrees with the exact form to better than 1e-9 out to a0 = 1e+06, and first exceeds a relative error of 1e-6 at a0 = 1e+12.

**Does this project's sampler reach the bad region?** a0 = (1 - F)/F, so a0 = 1e6 corresponds to F = 1e-6 and a0 = 1e17 to F = 1e-17. The basin posterior under the adopted Beta(1,10) prior has a median F of 0.064 with a 95 percent interval of [0.0215, 0.1346], which puts a0 between about 6 and 45. That is fifteen orders of magnitude away from anywhere the density degrades.

**Verdict: latent, not live.** The concern was legitimate, the arithmetic does degrade eventually, and no fit this project reports goes anywhere near it. The prior change recorded in D-31 makes this safer still, since Beta(1,10) puts even less mass near F = 0 than the flat prior it replaced.

**Why it was still worth measuring.** ../mataa's own first investigation of this concluded that cancellation was NOT the problem, and their second, six days later, found it catastrophic at a0 = 1e17 and traced a frozen chain to it. Both findings were correct at their own operating points. Ours is the benign one, and now it is measured rather than assumed.
