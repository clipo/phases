# Bayesian cultural F_ST for the basin (Balding-Nichols model)

Observed St. Francis basin (canonical drainage membership, 29 assemblages), 3 spatial clusters, 10 decorated types. Balding-Nichols Dirichlet-multinomial (F ~ Uniform(0,1), per-cluster frequencies marginalized), 2000 draws x 4 chains (full). This mirrors the hyperlocality project's Bayesian F_ST model.

## 1. Cultural F_ST parameter (Balding-Nichols theta)

- Posterior median F_ST = 0.0734 (mean 0.0813), 95% credible interval [0.0260, 0.1583].
- Flat prior on F_ST itself, so the estimate is not pulled toward the drift level by the prior.
- Prior sensitivity (Beta(1,3) prior on F): median 0.0709, 95% CI [0.0240, 0.1523].

## 2. Gini-Simpson F_ST (the manuscript estimator), conjugate readout

- Posterior median = 0.0179 (mean 0.0179), 95% credible interval [0.0157, 0.0201].
- Frequentist plug-in (variance.cultural_fst) = 0.0179, which falls inside the 95% interval.

This readout is a credible interval on the exact estimator the manuscript reports, reconstructed from the same fit as p_g | x_g ~ Dirichlet((1-F)/F * pi + x_g).

## 3. Structure vs panmixia (Bayes factor)

- log marginal likelihood: structure (M1) = -157.37, panmixia (M0) = -989.63.
- 2 ln BF10 = 1664.52 (+/- 0.16 across chains); log10 BF10 = 361.44.
- Kass & Raftery (1995): very strong for structure.

CAVEAT: this Bayes factor tests structure against EXACT panmixia (every cluster sharing one frequency vector). With thousands of sherds, exact panmixia is rejected trivially for any real assemblage set, so a large value here is expected and is NOT evidence of bounded groups. It is a much weaker null than drift: spatially structured neutral drift itself produces F_ST > 0 and would also reject panmixia. The paper's operative test is structure vs DRIFT, addressed by the stochastic-drift null (analyses 21/33/35/37) and the drift-versus-groups comparison, not by this panmixia Bayes factor. It is reported here for parity with the hyperlocality analysis and as a model-adequacy check, not as the structure test.

## MCMC diagnostics

- uniform prior (primary): max R-hat = 1.0012 (want < 1.01); min bulk ESS = 4502; min tail ESS = 4789; divergences = 0/8000 (0.00%); treedepth >= 10 = 0; min E-BFMI = 0.973 (want > 0.3).
- Beta(1,3) prior (sensitivity refit): max R-hat = 1.0011 (want < 1.01); min bulk ESS = 4376; min tail ESS = 4742; divergences = 0/8000 (0.00%); treedepth >= 10 = 0; min E-BFMI = 0.982 (want > 0.3).
