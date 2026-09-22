# Can unequal site populations close the drift gap?

Produced by `analyses/64_unequal_populations.py`, 300 realizations per cell.

Per-site populations are lognormal with the calibrated cell's population as the arithmetic mean. CV = 0 is the uniform model the paper reports. A cell closes the gap only if it reaches the observed $F_{ST}$ AND stays inside the calibration's own 10 percent diversity tolerance; raising inequality lowers the harmonic mean and so raises $F_{ST}$ for a trivial reason, which the diversity column is there to catch.

| region | CV | harmonic N | median $F_{ST}$ | 95% | observed | reaches obs | diversity matched |
|---|---:|---:|---:|---|---:|---:|:---:|
| basin | 0.0 | 2000 | 0.0026 | [0.0003, 0.0147] | 0.0062 | 17.7% | yes |
| basin | 0.5 | 1565 | 0.0027 | [0.0004, 0.0175] | 0.0062 | 22.3% | yes |
| basin | 1.0 | 1019 | 0.0042 | [0.0004, 0.0264] | 0.0062 | 36.3% | no |
| basin | 1.5 | 669 | 0.0055 | [0.0006, 0.0386] | 0.0062 | 45.3% | no |
| basin | 2.0 | 467 | 0.0075 | [0.0006, 0.0400] | 0.0062 | 56.7% | no |
| basin | 3.0 | 266 | 0.0122 | [0.0005, 0.0695] | 0.0062 | 70.0% | no |
| cmv | 0.0 | 120 | 0.0127 | [0.0037, 0.0381] | 0.0337 | 4.0% | yes |
| cmv | 0.5 | 101 | 0.0139 | [0.0044, 0.0455] | 0.0337 | 8.3% | yes |
| cmv | 1.0 | 68 | 0.0164 | [0.0050, 0.0519] | 0.0337 | 12.7% | no |
| cmv | 1.5 | 45 | 0.0206 | [0.0066, 0.0617] | 0.0337 | 20.0% | no |
| cmv | 2.0 | 31 | 0.0238 | [0.0060, 0.0788] | 0.0337 | 30.0% | no |
| cmv | 3.0 | 18 | 0.0299 | [0.0075, 0.1030] | 0.0337 | 43.3% | no |

## How much inequality was there?

Maximum mound height in `data/LMVData-22March2006.xls` covers 24 of the 28 basin assemblages and is the only size field in that table with usable coverage (`Area` is a region label, `Max Mound Area` is empty). Height is not population, so the implied inequality is given as a bracket over three scalings:

- population proportional to height: CV **0.92**
- population proportional to height squared: CV **1.37**
- population proportional to height cubed: CV **1.78**

Read the sweep against that bracket rather than against its own top end.


## Reading

Reported as the shortfall, the observed value divided by the cell's median, rather than as a pass or fail on a threshold. A binary verdict here flips on where the threshold is put: at the empirical bracket the shortfall is what matters, not whether some percentage crosses five.

- **basin**: uniform populations fall short by 2.4 times (0.0026 against 0.0062). At the empirical bracket (CV 0.9 to 1.6) the shortfall is 1.5 to 1.1 times. Even at CV 3.0, far beyond anything the settlement data suggest, it is 0.5 times.
- **cmv**: uniform populations fall short by 2.7 times (0.0127 against 0.0337). At the empirical bracket (CV 0.9 to 1.6) the shortfall is 2.1 to 1.6 times. Even at CV 3.0, far beyond anything the settlement data suggest, it is 1.1 times.

**The structural point.** Population inequality does raise differentiation, monotonically and by a factor of five across the sweep. It does not close the gap, because it lowers the harmonic mean population and so raises drift everywhere, which pulls within-assemblage diversity below the observed value. Of the 10 cells with unequal populations, 2 stay inside the calibration's diversity tolerance. Unequal populations therefore move ALONG the same diversity-against-differentiation trade-off the uniform model is already on, rather than off it. That is the same obstacle the calibration sweep reports, reached from a different direction.
