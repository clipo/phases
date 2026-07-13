# Bayesian cultural F_ST for the basin (Balding-Nichols model)

Observed St. Francis basin, 5 spatial clusters, 10 decorated types. Balding-Nichols Dirichlet-multinomial (F ~ Uniform(0,1), per-cluster frequencies marginalized), 2000 draws x 4 chains (full). This mirrors the hyperlocality project's Bayesian F_ST model.

## 1. Cultural F_ST parameter (Balding-Nichols theta)

- Posterior median F_ST = 0.0665 (mean 0.0702), 95% credible interval [0.0341, 0.1140].
- Flat prior on F_ST itself, so the estimate is not pulled toward the drift level by the prior.
- Prior sensitivity (Beta(1,3) prior on F): median 0.0657, 95% CI [0.0351, 0.1141].

## 2. Gini-Simpson F_ST (the manuscript estimator), conjugate readout

- Posterior median = 0.0327 (mean 0.0327), 95% credible interval [0.0301, 0.0352].
- Frequentist plug-in (variance.cultural_fst) = 0.0330, which falls inside the 95% interval.

This readout is a credible interval on the exact estimator the manuscript reports, reconstructed from the same fit as p_g | x_g ~ Dirichlet((1-F)/F * pi + x_g).

## 3. Structure vs panmixia (Bayes factor)

- log marginal likelihood: structure (M1) = -259.75, panmixia (M0) = -2181.61.
- 2 ln BF10 = 3843.72 (+/- 0.13 across chains); log10 BF10 = 834.65.
- Kass & Raftery (1995): very strong for structure.

CAVEAT: this Bayes factor tests structure against EXACT panmixia (every cluster sharing one frequency vector). With thousands of sherds, exact panmixia is rejected trivially for any real assemblage set, so a large value here is expected and is NOT evidence of bounded groups. It is a much weaker null than drift: spatially structured neutral drift itself produces F_ST > 0 and would also reject panmixia. The paper's operative test is structure vs DRIFT, addressed by the stochastic-drift null (analyses 21/33/35/37) and the drift-versus-groups comparison, not by this panmixia Bayes factor. It is reported here for parity with the hyperlocality analysis and as a model-adequacy check, not as the structure test.

## MCMC diagnostics (uniform-prior fit)

- max R-hat = 1.0007 (want < 1.01); min ESS = 4672; divergences = 0.
