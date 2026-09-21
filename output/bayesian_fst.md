# Bayesian cultural F_ST for the basin (Balding-Nichols model)

Observed St. Francis basin (canonical drainage membership, 43 assemblages), 5 spatial clusters, 10 decorated types. Balding-Nichols Dirichlet-multinomial (F ~ Beta(1,10), per-cluster frequencies marginalized), 2000 draws x 4 chains (full). The estimator is shared with a sibling project for consistency of method; no data is shared.

## 1. Cultural F_ST parameter (Balding-Nichols theta)

- Posterior median F_ST = 0.0617 (mean 0.0648), 95% credible interval [0.0315, 0.1028].
- Beta(1,10) prior on F_ST, adopted on prior-predictive grounds (analyses/57). It leans toward SMALL F_ST, which is the direction of this paper's conclusion, so the Uniform(0,1) refit below is the check that the prior is not doing the work.
- Prior sensitivity (Uniform(0,1) prior on F): median 0.0659, 95% CI [0.0339, 0.1139].

## 2. Gini-Simpson F_ST (the manuscript estimator), conjugate readout

- Posterior median = 0.0277 (mean 0.0277), 95% credible interval [0.0254, 0.0301].
- Frequentist plug-in (variance.cultural_fst) = 0.0280, which falls inside the 95% interval.

This readout is a credible interval on the exact estimator the manuscript reports, reconstructed from the same fit as p_g | x_g ~ Dirichlet((1-F)/F * pi + x_g).

## 3. Structure vs panmixia: not reported

The Bayes factor against exact panmixia that this section used to carry is **demoted and no longer computed** (F1). Two reasons.

First, its evidence band is prior-unstable exactly where it would matter. `../mataa` computed the same comparison exactly and measured the Kass-Raftery band moving across a boundary, 4.48 log units from Uniform(0,1) to Beta(1,10), for a near-null case. That is the Lindley-Bartlett effect on a nested boundary comparison and no estimator fixes it, so the quantity is robust where it is not needed and fragile where it would be load-bearing.

Second, exact panmixia is not the alternative anyone entertains. With thousands of sherds it is rejected trivially by any real assemblage set, and spatially structured neutral drift, which IS the alternative this paper tests, produces F_ST > 0 and would reject it too. The operative comparison is structure against drift, addressed by the drift-versus-groups analysis, not by a panmixia Bayes factor.

What stands in its place is section 1: the posterior for F itself, under a prior justified by prior predictive (analyses/57) and reported under two priors leaning opposite ways.

## MCMC diagnostics

- Beta(1,10) prior (primary): max R-hat = 1.0012 (want < 1.01); min bulk ESS = 5123; min tail ESS = 4419; divergences = 0/8000 (0.00%); treedepth >= 10 = 0; min E-BFMI = 0.964 (want > 0.3).
- Uniform(0,1) prior (sensitivity refit): max R-hat = 1.0010 (want < 1.01); min bulk ESS = 4375; min tail ESS = 4450; divergences = 0/8000 (0.00%); treedepth >= 10 = 0; min E-BFMI = 0.955 (want > 0.3).
