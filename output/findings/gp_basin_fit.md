# The basin spatial GP: the fit behind the reported variance share

Produced by `analyses/02_spatial/05_basin_fit.R` on 2026-09-23.

Basin contract `data/stan/basin_composition.json`: 28 assemblages, 10 decorated classes.

`spatial_share` is the fraction of assemblage-level compositional variance
carried by the spatial field rather than by the assemblage-specific term. It
is the quantity the abstract reports. The length scale `rho` is reported here
for completeness and is NOT reportable as a number: see the verdict below.

| model | rung | passed | verdict | spatial_share | 95% CI | rho median | rho 95% CI |
|---|---|---|---|---|---|---|---|
| composition_gp | 2 | **no** | ridge | 0.997 | [0.975, 1.000] | 195.4 | [94.3, 292.1] |
| composition_gp_marginal | 2 | **no** | ridge | 0.996 | [0.972, 1.000] | 195.3 | [95.6, 292.2] |

## Diagnostics (rule 16: reported for every fit, without exception)

- composition_gp (rung 2): R-hat 1.0157 | bulk ESS 438 | tail ESS 5660 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.793
- composition_gp_marginal (rung 2): R-hat 1.0020 | bulk ESS 2211 | tail ESS 7241 | divergences 0/12000 (0.00%) | treedepth>=10 12000 | E-BFMI 0.688

## Reading

The marginalised fit is the one the paper reports. It does not pass the rule-16 gate: treedepth saturation is 100 percent with 0 divergences. Under `diag_verdict` that is a **ridge**, not a step-size problem, and the ladder therefore stops escalating rather than spending compute on a geometry that will not respond to it (rule 17 distinguishes an identifiability limit from a fixable parameterisation).

The cause is established elsewhere rather than assumed here. `gp_geometry_repair.md` shows the saturation is a property of the ridge and not of the additive-versus-marginal parameterisation, and `gp_posterior_predictive.md` shows the model reproduces cultural F_ST, the distance decay and the spread of within-assemblage diversity, so misspecification is excluded as the explanation. What remains is that `rho` is flat above the design's resolving ceiling: the basin spans about 128 km and the posterior piles up against the top of its range, which is the data saying no decay in composition is detectable across the drainage.

The consequence for the manuscript is the split it already makes. `spatial_share` is identified and is reported with its interval. A length-scale NUMBER is not identified and is reported as unresolved, never as an estimate. A reader is entitled to both the share and the diagnostics above, which is why this file exists.

