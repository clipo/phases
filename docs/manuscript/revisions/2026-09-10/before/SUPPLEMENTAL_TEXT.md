# Supplemental Material

**Are the phases real? Testing bounded interaction against spatially structured drift in Mississippi Valley decorated ceramics**

**Carl P. Lipo, Robert J. DiNapoli, and Mark E. Madsen**

## Measurement definitions

For counts $n_i$ in an assemblage of total size $N$, the unbiased sample homozygosity is

$$F=\frac{\sum_i n_i(n_i-1)}{N(N-1)}.$$

The homozygosity-based estimate of the neutral transmission parameter is $\theta_F=(1-F)/F$. The richness-based estimate $\theta_E$ solves the Ewens sampling expectation for the number of observed classes $k$ [@ewens_1972; @neiman_1995]:

$$k=\sum_{i=0}^{N-1}\frac{\theta_E}{\theta_E+i}.$$

The empirical and injection analyses summarize $|1-\theta_F/\theta_E|$ within spatial clusters and average over clusters represented in each bin. The absolute value measures departure and does not identify its direction. The idealized profile experiment instead uses the signed within-row departure $1-\theta_F/\theta_E$. These are different diagnostic choices. Finite archaeological class aggregation, pooling, and time-averaging limit interpretation of either statistic against infinite-alleles expectations [@madsen_2012; @crema_bortolini_lake_2023].

The cultural variance partition uses Gini-Simpson diversity $H=1-\sum_i p_i^2$ and

$$F_{ST}=\frac{H_T-H_S}{H_T}.$$

$H_T$ is diversity in pooled counts, and $H_S$ is within-group diversity weighted by group sherd totals. This weighting is equivalent to assemblage-number weighting after equal-depth rarefaction, but not on the original unequal-size assemblages. A single group has zero between-group differentiation. Simulated single-community outcomes are retained rather than excluded when summarizing the null. Distinct partitions can yield different $F_{ST}$ values for the same count matrix; community, spatial, and named-phase estimates are labeled separately throughout.

Brainerd-Robinson similarity is $S_{ij}=200-\sum_k|100p_{ik}-100p_{jk}|$ [@brainerd_1951; @robinson_1951]. Boundary excess compares within- and between-cluster similarities among pairs at comparable distances. Its empirical implementation includes an internal spatial clustering rule, so it is a descriptive diagnostic rather than a direct observation of social membership. A positive value is not sufficient to establish a boundary, because finite spatial drift can also yield excess.

The deterministic IDSS implementation enforces unimodal columns, adjacent co-presence of at least one class, and a maximum adjacent per-class frequency jump of 0.10. It begins from valid pairs, grows orders at their ends, collapses orderings to sets, and retains maximal sets. It reports a warning if the search budget or solution cap is exceeded. Full-basin group counts describe the original assemblages. In the revised injection experiment, the statistic is the maximal group count divided by the number of assemblages in each bin. Neither operation includes the sampling-significance layer of the published IDSS method [@lipo_madsen_dunnell_2015].

## Revised injection experiment

Let $p_0$ be a Zipf profile over the ten measured classes and $q_g$ a cluster-specific permutation of that profile. For bin $t=0,\ldots,5$, define $a=st/5$ and

$$p_{gt,i}\propto\left[(1-a)p_{0,i}+a q_{g,i}\right]^{1+1.2a}.$$

This construction makes divergence and concentration increase together. It is a profile-injection experiment, not a Wright–Fisher trajectory. At $s=0$ every group has the same constant profile. Production is multinomial at the observed assemblage totals; the resulting counts are then rarefied without replacement to 50. For the averaging scenario, each group's profile is averaged over the current and up to two preceding bins before sampling. The shorter initial profile windows are an explicit feature of this injection scenario; the separate spatial model uses full windows at every observed rank.

Each strength is evaluated using 120 independent test draws. An additional 400 null draws, independent of the test draws, establish each one-sided 95th-percentile threshold. Because rank statistics are discrete and calibration is finite, the achieved false-positive rate need not be exactly 0.05. The two observed rates are 0.04 and 0.07. The first grid strengths achieving at least 80% power are 0.6 and 0.7 for the unaveraged and averaged scenarios. Values between grid points are not resolved.

The empirical mean $F_{ST}$ rank trend is +0.02, with rarefaction percentiles -0.30 to +0.60. Those percentiles are not an interval for $s$. Neither the empirical mean nor its placement against the mean response curve establishes an excluded closure strength. The injected signal is increasing divergence; stable differentiation is evaluated by the separate spatial comparison. A non-spatial or conformity-only mechanism need not be detected by this experiment.

## Idealized four-proxy experiment

The earlier framework experiment uses four constructed generators: coupled divergence and concentration, pooled conformity, static patchiness, and a smooth spatial gradient. Their abundance profiles are deliberately chosen to give a particular neutral-reading baseline. They do not sample the full distribution of evolving stochastic spatial drift. The inference rule declares convergence when all four standardized ordinal slopes exceed 0.10, with scales computed across the four generators.

One signature in this experiment is the negative number of fall-to-rise reversals in the supplied row order. The revised score correctly treats a valley as a violation and ignores flat steps. It is not the maximal IDSS group count used in the revised record-matched analysis. Correcting this implementation yields 493 convergence hits in 500 coupled-profile draws and zero in 500 draws from each other generator. The two-sided exact 95% interval for the coupled-profile hit rate is 0.971–0.994; for zero hits in 500 draws, the two-sided upper bound is 0.0074. These intervals quantify Monte Carlo performance for the specified generators, not general diagnostic sensitivity and specificity for archaeological societies.

![**Figure S1.** Revised idealized four-proxy experiment. Points show proportions satisfying the joint slope rule in 500 realizations of each specified profile generator; bars are exact two-sided 95% binomial intervals. The ordered-row seriation proxy differs from the IDSS estimator in the main recovery experiment.](../../figures/revision_idealized_validation.png){width=6in}

## Spatial update and observation model

Write $p_{it}$ for the vector of class frequencies at node $i$ and generation $t$, and $W$ for the row-normalized external-copying matrix. The next-generation probabilities are

$$q_{it}=(1-\mu)\left[(1-m)p_{it}+m\sum_j W_{ij}p_{jt}\right]+\mu/K,$$

followed independently at each node by a multinomial draw of 120 learners. This frequency-state update integrates over individual source choices and is distributionally equivalent to independent neutral learner copying under these probabilities. It does not include selection among classes or a conformist exponent. Innovation draws from the existing archaeological class repertoire; genuinely new class introduction is outside this model.

For geographic or river distances $d_{ij}$ and a proposed partition $g$, unnormalized external weights are

$$w_{ij}=\exp(-d_{ij}/L)\begin{cases}1,&g_i=g_j,\\\ell,&g_i\ne g_j,\end{cases}\qquad i\ne j.$$

Self weights are zero before normalization; local copying is represented by $1-m$. At $\ell=1$ the weight matrix is exactly the distance-only matrix. In the limiting case of an isolated singleton under a zero-leak boundary, an external-copy event defaults to self copying. The reported grid does not use zero leak.

River-network distance is constructed from the mapped LMVHydrology centerlines clipped to the site extent with a 25 km margin. Vertices within 250 m are connected, the largest connected component is retained, and sites are snapped to its nearest nodes. The path cost is shortest along-network length plus access costs at both ends. This is a reproducible transport approximation, not a reconstruction of every late pre-contact channel. The displayed regional maps use HydroRIVERS for cartographic consistency; the simulation distances use the LMVHydrology graph. Southeast-Missouri simulations use geographic distance rather than this river graph.

The main model uses each dataset's observed number of columns, 120 learners per node, mixing $m=0.02$, innovation $\mu=0.012$, and interaction length $L=24$ km. It starts from uniform frequencies and discards 1,200 generations, then records 90 evenly spaced slices over 1,800 further generations. A rank $r\in[0,1]$ maps to slice

$$u=(w-1)+\operatorname{round}\{r(R-w)\},$$

where $R=90$ and $w=8$. Each observed assemblage pools slices $u-w+1$ through $u$, and a multinomial observation draw matches its real sherd total. The earliest rank therefore receives a full window. The contemporaneous observation sets every rank to one and retains the same window length; it is not an unaveraged snapshot.

CA ranks specify relative sampling positions but do not calibrate the number of generations between archaeological assemblages. The reversed-order scenario addresses sign uncertainty, not full chronological uncertainty. Initial-condition tests compare a 2,400-generation burn-in, monomorphic initialization with the standard burn-in, and monomorphic initialization without burn-in. A 40-class sensitivity changes repertoire size while preserving the other main settings; it is a diagnostic of class dependence rather than a matched archaeological observation model.

## Baseline and matched boundary results

**Table S1.** Baseline comparisons, conditional on each stated partition and sampling scenario.

| Contrast | Sampling | Observed | Median [95%] | Upper-tail p |
|---|---|---:|---|---:|
| Basin spatial | Same time | 0.0179 | 0.0085 [0.0036, 0.0186] | 0.038 |
| Basin spatial | CA order | 0.0179 | 0.0061 [0.0023, 0.0145] | 0.004 |
| Basin spatial | Reverse | 0.0179 | 0.0062 [0.0025, 0.0140] | 0.006 |
| CMV spatial | Same time | 0.0337 | 0.0078 [0.0057, 0.0113] | 0.002 |
| CMV spatial | CA order | 0.0337 | 0.0067 [0.0049, 0.0093] | 0.002 |
| CMV spatial | Reverse | 0.0337 | 0.0067 [0.0049, 0.0093] | 0.002 |
| Parkin named | Same time | 0.0121 | 0.0061 [0.0021, 0.0142] | 0.060 |
| Parkin named | CA order | 0.0121 | 0.0042 [0.0014, 0.0113] | 0.022 |
| Parkin named | Reverse | 0.0121 | 0.0043 [0.0014, 0.0108] | 0.016 |
| Valley communities | Same time | 0.0413 | 0.0071 [0.0037, 0.0139] | 0.002 |
| Valley communities | CA order | 0.0413 | 0.0070 [0.0040, 0.0132] | 0.002 |
| Valley communities | Reverse | 0.0413 | 0.0070 [0.0039, 0.0123] | 0.002 |

Baseline intervals use 500 realizations per region. The upper-tail probability is $(e+1)/(R+1)$; no probability smaller than $1/501$ can be resolved by this ensemble. The same latent realization is observed in three sampling configurations. Therefore these columns are sensitivity comparisons, not independent replications of archaeological evidence.

**Table S2.** Number of length/innovation cells bracketing the observed partition statistic (six cells per restriction multiplier).

| Region | Multiplier | Cells bracketing observed / total |
|---|---:|---:|
| basin | 1 | 3 / 6 |
| basin | 0.5 | 2 / 6 |
| basin | 0.1 | 2 / 6 |
| basin | 0.03 | 2 / 6 |
| valley | 1 | 2 / 6 |
| valley | 0.5 | 2 / 6 |
| valley | 0.1 | 2 / 6 |
| valley | 0.03 | 2 / 6 |
| cmv | 1 | 0 / 6 |
| cmv | 0.5 | 0 / 6 |
| cmv | 0.1 | 1 / 6 |
| cmv | 0.03 | 1 / 6 |

The grid uses 30 realizations per combination of region, length, innovation, and restriction multiplier. Table entries count the cells whose 2.5th–97.5th simulation interval brackets the observed partition statistic. This count describes sensitivity to the tested settings, not a model probability or a multiple-testing decision. The parameter settings receive no prior weights. A small cell count does not establish that a mechanism is impossible outside the grid.

**Table S3.** Repertoire and initialization sensitivity (50 realizations per case).

| Region | Scenario | Median [95% interval] | Upper-tail p |
|---|---|---|---:|
| basin | fixed 40 classes | 0.0064 [0.0042, 0.0094] | 0.020 |
| basin | longer burnin | 0.0064 [0.0023, 0.0137] | 0.039 |
| basin | monomorphic burnin | 0.0065 [0.0027, 0.0128] | 0.020 |
| basin | monomorphic no burnin | 0.0069 [0.0036, 0.0135] | 0.020 |
| valley | fixed 40 classes | 0.0048 [0.0028, 0.0082] | 0.020 |
| valley | longer burnin | 0.0045 [0.0012, 0.0099] | 0.020 |
| valley | monomorphic burnin | 0.0043 [0.0016, 0.0095] | 0.020 |
| valley | monomorphic no burnin | 0.0041 [0.0016, 0.0095] | 0.039 |

The repertoire and initialization checks use 50 realizations per case. Their limited resolution is reported separately from the larger baseline. All exact numerical replicate values are retained in the accompanying CSV files.

## Exploratory pooled-population analyses

Earlier analyses fitted a frequency-dependent copying model with adoption weights proportional to $p_i^{1+b}$, where $b=0$ is neutral, $b>0$ conformist, and $b<0$ anti-conformist. That model pools production across short windows in a single population, whereas the observations pool spatial assemblages into six CA-derived bins. Its synthetic observations do not replicate the observed per-bin sherd totals, spatial heterogeneity, or the CA-ordering step. It therefore answers a different, more restrictive model-conditional question than the revised spatial comparison.

The existing ABC-SMC fit uses 800 particles and six rounds, followed by local-linear regression adjustment [@toni_etal_2009; @beaumont_zhang_balding_2002]. Its archived mean is $b=+0.006$, with an interval of $-0.020$ to $+0.035$. These are retained as historical exploratory outputs, not new estimates from the revised pipeline and not evidence excluding strong conformity.

The calibration battery itself uses 300 particles and five rounds. Across 200 prior-predictive datasets, the rank-uniformity test gives $D=0.068$, $p=0.305$. However, conditional coverage is inadequate in at least one tested region of parameter space. The nominal 95% intervals contain true $b=+0.10$ in only 13 of 20 datasets. A prior-averaged rank check cannot resolve that undercoverage. Simulation-based calibration evaluates inference under the assumed simulator; it does not verify that the simulator captures the archaeological observation process [@talts_etal_2018].

**Table S4.** Archived conditional coverage, shown to qualify the earlier inference. Each row uses 20 simulated datasets and differs from the empirical sampler configuration.

| True $b$ | Fraction covered by nominal 95% interval |
|---:|---:|
| −0.40 | 0.90 |
| −0.30 | 0.85 |
| −0.20 | 0.85 |
| −0.10 | 0.95 |
| 0.00 | 0.90 |
| +0.10 | 0.65 |
| +0.20 | 1.00 |
| +0.30 | 1.00 |
| +0.40 | 0.90 |

Before the pooled ABC analysis can carry an exclusion claim, it requires a matched observation model, calibration at the empirical inference settings, adequate local coverage, and posterior predictive checks. This revision does not rely on it to adjudicate the drift-versus-boundary question. Likewise, the six-point tempo-and-mode comparison and the hierarchical convergence analysis are retained in the research archive but do not provide an independent resolution of the observation-model limitations identified here.

## Radiocarbon chronology

The Mainfort compilation contains 40 determinations from seven basin sites assigned to the Parkin phase, including 19 determinations from Parkin [@mainfort_2001]. Calibration against IntCal20 gives a basin summed-probability median of AD 1428 and a Parkin median of AD 1483 [@reimer_etal_2020]. About 29% of the basin probability mass postdates AD 1541 and 18% postdates AD 1600. These summaries describe the collected determinations; their distribution should not be equated with an occupation or demographic trajectory without modeling sampling and calibration effects.

Five dated proveniences match curated assemblages, contributing 27 pooled determinations to the CA-calendar comparison. Their median calibrated dates correlate with CA1 at Spearman $\rho=+0.70$; Clay Hill departs from the predominant order. The limited anchor motivates the orientation sensitivity and conditional language in the main text. Adding twelve net-new Mississippian determinations from the compilation supplied by Jeffrey Alvey (personal communication, 2022) shifts the basin median by eight years, to AD 1436, without changing the reported rank correlation. The determinations and sources are listed in `radiocarbon_dates_used.csv`. These chronology summaries are retained from the existing calibration analysis; the revision changes their interpretation rather than recalibrating the determinations.

![**Figure S2.** Existing IntCal20 chronology summary. Left: summed calibrated probability of the basin determinations and Parkin dates, with AD 1541 marked. Right: CA1 against pooled median calibrated age for five matching proveniences. The sparse association supports an orientation check but does not uniquely identify CA1 as calendar time.](../../figures/figS7_chronology.png){width=6in}

## Reproducibility and revision scope

The canonical revision scripts are `21_signal_recovery.py`, `47_revision_analysis.py`, and `48_revision_validation.py`. Scripts 33–37 are earlier model experiments, not the source of the revised spatial figures. The revised simulator uses fixed measured classes, explicit burn-in, complete observation windows, and a shared update for both drift and restriction. Its output directory contains baseline, boundary-grid, initialization, and actual leave-one-out replicate tables plus a configuration and input/code fingerprint. A cache is reused only when that fingerprint matches. The `--force` option recomputes the spatial ensembles.

The pre-revision manuscript and its figure images are preserved in the revision history. The revised manuscript is assembled from explicit results summaries; no numerical result is selected to restore a prior historical conclusion. Original descriptive maps, whole-basin IDSS structure, mound-height ranks, and radiocarbon figure are retained. The revised empirical trajectories, recovery, spatial comparisons, and idealized proxy check are regenerated. The public release repository and archive must be updated separately before submission.

## References
<!-- managed in references.bib -->
