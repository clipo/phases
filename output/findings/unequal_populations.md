# Can unequal site populations close the drift gap?

Produced by `analyses/64_unequal_populations.py`, 300 realizations per cell.

Per-site populations are lognormal with the calibrated cell's population as the arithmetic mean. CV = 0 is the uniform model the paper reports. A cell closes the gap only if it reaches the observed $F_{ST}$ AND stays inside the calibration's own 10 percent diversity tolerance; raising inequality lowers the harmonic mean and so raises $F_{ST}$ for a trivial reason, which the diversity column is there to catch.

| region | CV | harmonic N | median $F_{ST}$ | 95% | observed | reaches obs | diversity matched |
|---|---:|---:|---:|---|---:|---:|:---:|
| basin | 0.0 | 10000 | 0.0030 | [0.0005, 0.0145] | 0.0179 | 1.0% | yes |
| basin | 0.5 | 7814 | 0.0036 | [0.0004, 0.0171] | 0.0179 | 2.0% | no |
| basin | 1.0 | 5097 | 0.0053 | [0.0007, 0.0223] | 0.0179 | 6.0% | no |
| basin | 1.5 | 3361 | 0.0073 | [0.0014, 0.0363] | 0.0179 | 17.0% | no |
| basin | 2.0 | 2349 | 0.0115 | [0.0013, 0.0645] | 0.0179 | 33.0% | no |
| basin | 3.0 | 1342 | 0.0159 | [0.0016, 0.1064] | 0.0179 | 44.7% | no |
| cmv | 0.0 | 120 | 0.0117 | [0.0035, 0.0374] | 0.0337 | 3.3% | yes |
| cmv | 0.5 | 101 | 0.0127 | [0.0044, 0.0442] | 0.0337 | 5.3% | yes |
| cmv | 1.0 | 68 | 0.0155 | [0.0035, 0.0630] | 0.0337 | 14.0% | no |
| cmv | 1.5 | 45 | 0.0180 | [0.0051, 0.0626] | 0.0337 | 20.3% | no |
| cmv | 2.0 | 31 | 0.0214 | [0.0066, 0.0808] | 0.0337 | 24.7% | no |
| cmv | 3.0 | 18 | 0.0263 | [0.0059, 0.0855] | 0.0337 | 36.0% | no |

## How much inequality was there?

Maximum mound height in `data/LMVData-22March2006.xls` covers 23 of the 29 basin assemblages and is the only size field in that table with usable coverage (`Area` is a region label, `Max Mound Area` is empty). Height is not population, so the implied inequality is given as a bracket over three scalings:

- population proportional to height: CV **0.88**
- population proportional to height squared: CV **1.27**
- population proportional to height cubed: CV **1.61**

Read the sweep against that bracket rather than against its own top end.


## Reading

Reported as the shortfall, the observed value divided by the cell's median, rather than as a pass or fail on a threshold. A binary verdict here flips on where the threshold is put: at the empirical bracket the shortfall is what matters, not whether some percentage crosses five.

- **basin**: uniform populations fall short by 6.0 times (0.0030 against 0.0179). At the empirical bracket (CV 0.9 to 1.6) the shortfall is 3.4 to 2.4 times. Even at CV 3.0, far beyond anything the settlement data suggest, it is 1.1 times.
- **cmv**: uniform populations fall short by 2.9 times (0.0117 against 0.0337). At the empirical bracket (CV 0.9 to 1.6) the shortfall is 2.2 to 1.9 times. Even at CV 3.0, far beyond anything the settlement data suggest, it is 1.3 times.

**The structural point.** Population inequality does raise differentiation, monotonically and by a factor of five across the sweep. It does not close the gap, because it lowers the harmonic mean population and so raises drift everywhere, which pulls within-assemblage diversity below the observed value. Of the 10 cells with unequal populations, 1 stays inside the calibration's diversity tolerance. Unequal populations therefore move ALONG the same diversity-against-differentiation trade-off the uniform model is already on, rather than off it. That is the same obstacle the calibration sweep reports, reached from a different direction.
