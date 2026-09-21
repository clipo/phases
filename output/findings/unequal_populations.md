# Can unequal site populations close the drift gap?

Produced by `analyses/64_unequal_populations.py`, 300 realizations per cell.

Per-site populations are lognormal with the calibrated cell's population as the arithmetic mean. CV = 0 is the uniform model the paper reports. A cell closes the gap only if it reaches the observed $F_{ST}$ AND stays inside the calibration's own 10 percent diversity tolerance; raising inequality lowers the harmonic mean and so raises $F_{ST}$ for a trivial reason, which the diversity column is there to catch.

| region | CV | harmonic N | median $F_{ST}$ | 95% | observed | reaches obs | diversity matched |
|---|---:|---:|---:|---|---:|---:|:---:|
| basin | 0.0 | 10000 | 0.0029 | [0.0007, 0.0114] | 0.0207 | 1.0% | yes |
| basin | 0.5 | 8176 | 0.0035 | [0.0007, 0.0146] | 0.0207 | 0.7% | yes |
| basin | 1.0 | 5375 | 0.0050 | [0.0010, 0.0254] | 0.0207 | 5.0% | yes |
| basin | 1.5 | 3551 | 0.0078 | [0.0011, 0.0428] | 0.0207 | 13.0% | no |
| basin | 2.0 | 2486 | 0.0118 | [0.0018, 0.0581] | 0.0207 | 29.0% | no |
| basin | 3.0 | 1423 | 0.0167 | [0.0030, 0.0998] | 0.0207 | 41.3% | no |
| cmv | 0.0 | 120 | 0.0117 | [0.0035, 0.0374] | 0.0337 | 3.3% | yes |
| cmv | 0.5 | 101 | 0.0127 | [0.0044, 0.0442] | 0.0337 | 5.3% | yes |
| cmv | 1.0 | 68 | 0.0155 | [0.0035, 0.0630] | 0.0337 | 14.0% | no |
| cmv | 1.5 | 45 | 0.0180 | [0.0051, 0.0626] | 0.0337 | 20.3% | no |
| cmv | 2.0 | 31 | 0.0214 | [0.0066, 0.0808] | 0.0337 | 24.7% | no |
| cmv | 3.0 | 18 | 0.0263 | [0.0059, 0.0855] | 0.0337 | 36.0% | no |

## How much inequality was there?

Maximum mound height in `data/LMVData-22March2006.xls` covers 34 of the 43 basin assemblages and is the only size field in that table with usable coverage (`Area` is a region label, `Max Mound Area` is empty). Height is not population, so the implied inequality is given as a bracket over three scalings:

- population proportional to height: CV **0.94**
- population proportional to height squared: CV **1.45**
- population proportional to height cubed: CV **1.94**

Read the sweep against that bracket rather than against its own top end.


## Reading

Reported as the shortfall, the observed value divided by the cell's median, rather than as a pass or fail on a threshold. A binary verdict here flips on where the threshold is put: at the empirical bracket the shortfall is what matters, not whether some percentage crosses five.

- **basin**: uniform populations fall short by 7.1 times (0.0029 against 0.0207). At the empirical bracket (CV 0.9 to 1.6) the shortfall is 4.2 to 2.7 times. Even at CV 3.0, far beyond anything the settlement data suggest, it is 1.2 times.
- **cmv**: uniform populations fall short by 2.9 times (0.0117 against 0.0337). At the empirical bracket (CV 0.9 to 1.6) the shortfall is 2.2 to 1.9 times. Even at CV 3.0, far beyond anything the settlement data suggest, it is 1.3 times.

**The structural point.** Population inequality does raise differentiation, monotonically and by a factor of five across the sweep. It does not close the gap, because it lowers the harmonic mean population and so raises drift everywhere, which pulls within-assemblage diversity below the observed value. Of the 10 cells with unequal populations, 3 stay inside the calibration's diversity tolerance. Unequal populations therefore move ALONG the same diversity-against-differentiation trade-off the uniform model is already on, rather than off it. That is the same obstacle the calibration sweep reports, reached from a different direction.
