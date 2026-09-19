# Simulated critical peer review — mls-emergence

Four adversarial referees (cultural-transmission methods; CMV archaeology; multilevel-selection theory; quantitative/statistics) read the full manuscript on 2026-06-10. **All four recommend major revision.** None recommend rejection; all credit the convergence criterion and the willingness to publish an honest negative. The striking feature is convergence *among the reviewers*: the same weaknesses surface from four independent angles.

---

## Cross-cutting findings (raised by 2+ reviewers)

### A. The negative result is underpowered, and "sharp"/"no convergence" overstates it (R1-2.5, R4-M1, R4-M2)
The load-bearing convergence-score trend is rho = −0.60, **p = 0.208** (non-significant); two of the four signatures (F_ST, spatial boundary) carry no resolvable trend at n = 53 in 6 bins. "Absence of evidence" is being read as "evidence of absence" without a power analysis. Fix: power/sensitivity analysis at the realized configuration, or soften to "the resolvable signatures do not converge; two are underdetermined."

### B. Internal inconsistency: F_ST and the spatial signature are mischaracterized (R1-2.4, R3-minor, R4-M3) — FACTUAL
Abstract says "between-group variance **decline**"; Results says F_ST is "**flat** and small." The basin F_ST is rho = −0.143, CI includes 0 (flat). The spatial boundary is called "flat and noisy" but **rises** (+0.60) in the whole-LMV run and **flips sign** with bin count (+0.143 at 6 bins, −0.619 at 8). A sign-unstable signature cannot be summarized as "flat."

### C. Researcher degrees of freedom / garden of forking paths (R4-M3, R1-minor, R2-2.1)
Result depends on latitude cut (34.5), bin count (4/6/8 flip signs), IDSS continuity threshold (51/124/683 groups at cont 0.05/0.1/0.2), and k-means k (silhouette nearly flat, k=4 wins by 0.003 over k=2). Report the full sensitivity grid; present "124 groups" as threshold-conditional and lead with the *bridge topology* (stable) not the count.

### D. The chronological axis orientation is not established (R1-2.7, R2-minor, R3-minor, R4-M4)
CA↔14C Spearman = +0.50, **p = 0.391** on 5 anchors. Every "declines toward contact" statement depends on an unverified temporal polarity; flip the axis and neutrality *rises* (the emergence direction). Fix: separate strict ordinal-gradient language from temporal claims, and state explicitly that the *dissociation* conclusion is orientation-invariant while *directional* claims are not.

### E. The validation generators are weak/circular; specificity overstated (R1-2.3, R3-2.5, R4-M5)
The three mimics are parameterized to hard-zero the signature each is meant to fail, so "only genuine emergence lights all four" is close to definitional. The coupling model sets between-group divergence = within-group conformity = a = cφ, building convergence into the generator (so the 0.88 audit is an artifact). 0/20 specificity has a ~14% upper bound (rule of three). Fix: add a hard adversarial mimic (strong assortment *without* group selection), more seeds with CIs, and time-averaging matched to the empirical data.

### F. Independence is shown in simulation, never in the data (R1-2.2, R3-3.1, R4-M6)
All four signatures derive from the same 10-type × 53 matrix; θ and F_ST share the same k-means clusters. Report the **empirical** pairwise correlation matrix of the four signatures.

### G. The loudest empirical signal (neutrality −0.94) is the one most exposed to the time-averaging confound the paper itself names (R1-2.1, R2-2.7, R3-minor)
θ_F/θ_E is built from exactly the richness/evenness quantities Madsen 2012 / Premo 2014 say time-averaging distorts. A monotone decline could be a taphonomic/aggregation artifact. Fix: a time-averaging sensitivity test before resting a unit-of-selection conclusion on it.

### H. "These limits cut against a false negative rather than toward one" is backwards (R1-2.5, R2-2.7, R4-M4) — flagged by THREE reviewers
Time-averaging, small n, and ordinal chronology all *reduce* detection power, which makes a false negative *more* likely, not less. This sentence needs a real argument or removal.

---

## Reviewer-specific majors

### R1 — Cultural-transmission methods (Major revision)
- θ_F/θ_E invalid on time-averaged data without disentangling the confound (2.1, 2.8: object-mediated null not re-derived).
- Convergence/independence validated only on engineered generators, not the shared-substrate empirical matrix (2.2, 2.3).
- Empirical dissociation overstated; F_ST mislabeled "decline"; spatial signature sign-unstable (2.4).
- Underpowered negative (p = 0.208); coupling-robustness argument only holds for the *relaxed* combined score, not the strict four-signature rule (the strict rule fails 25-50% even at full coupling) (2.5, 2.6).
- CA-axis-as-time partly circular with signature 2 (both are type-frequency ordinations) (2.7).

### R2 — CMV archaeology (Major revision)
- Latitude cut is a statistical convenience, not a cultural/hydrological boundary; excludes flagged St-Francis sites (Salomon, Parchman); primacy is an artifact of where the box is drawn (2.1).
- **Most serious:** the rank-size primate claim substitutes Parkin's ~17-acre *site* area into a *mound-area* field — apples-to-oranges; primacy 41 not publishable as a measurement; rest on mound height instead (2.2).
- PFG broad set is a mound-biased *collection* sample, not a settlement census; "110 sites" rank-size and feature proportions misrepresent it (2.3).
- Parkin = Casqui / De Soto chronicles leaned on too hard as a "threshold"; Casqui-Pacaha conflict (an H1/H2-relevant point) unmentioned (2.4).
- "Three independent lines" overstated: Mainfort 2003 NMDS shares the data substrate (not independent); settlement compromised; Cobb sound but over-extended to the basin (2.5).
- **Ross-Stallings 2015 is a *Mississippi Delta* (lower Yazoo) study** — outside the study area, imported from the very region the latitude cut excluded; verify and re-site or drop (2.6). [FACTUAL]

### R3 — Multilevel-selection theory (Major revision; recommends "Path A" retreat)
- β₀/β₁ notation appears only in the Figure 2 caption and Discussion, never defined by an equation; imports the look of contextual-analysis MLS with none of its content (2.1).
- The four signatures detect *boundedness of interaction*, not a change in the *level of selection*; the group-fitness covariance term is absent, so even a positive result wouldn't license "unit of selection" (2.2).
- "Assortment becoming binding" never operationalized as a threshold; a bistable prediction is tested with a smooth-ramp rank correlation — mismatched (2.3).
- Price-equation critique is a strawman of Okasha & Otsuka (who show the partition *can* be given causal content) (2.4).
- The coupling identification (divergence = conformity = a = cφ) smuggles in the conclusion; the bistable parent model does no empirical work and is decorative — remove it and nothing changes (2.5).
- Conflates Smaldino's "emergent group-level trait" with the transitions literature's "new unit of selection"; "archaeology can observe it" overstated (model-mediated inference) (2.6, 2.7).
- **Recommends Path A:** retreat to detecting *bounded assortative interaction communities*, drop "unit of selection" and β-notation, present the signaling model as motivation not derivation.

### R4 — Quantitative / statistics (Major revision)
- Underpowered; headline rho = −0.60 is p = 0.208; no power analysis at the empirical configuration (M1).
- 6-point rank correlations over-interpreted; report p-values/CIs; make per-assemblage bootstrap slopes (n = 53) the primary inferential object, demote 6-bin rho to descriptive (M2).
- Forking paths uncontrolled; bin-count *sign instability* of F_ST and spatial boundary should be foregrounded, not buried (M3).
- Axis orientation unestablished (p = 0.391); strip temporal-direction language or defend it; state dissociation is orientation-invariant (M4).
- Validation regime (clean, G=12, N=300) doesn't match the empirical regime; 20 seeds too few; 0/20 specificity upper bound ~14% — report Clopper-Pearson, raise seeds to 200+ (M5).
- Report empirical signature correlation matrix (M6).
- Table 1 is referenced 3× but not present in the manuscript. [FACTUAL]

---

## Triage for the authors

**Quick factual/consistency fixes (do now, no science change):**
1. Abstract "between-group variance decline" → "flat" (matches Results, the data). [B]
2. Create the missing **Table 1** (signature × mechanism matrix), referenced 3×. [R4]
3. Verify/re-site **Ross-Stallings 2015** (Mississippi Delta, not St. Francis); flag or drop. [R2-2.6]
4. Remove or soften "these limits cut against a false negative" — three reviewers call it backwards. [H]
5. Soften overstated quantifiers: "sharp" negative, "confirms," "only," "primate settlement center by ... site area." [A, R4]
6. Attach p-values/CIs to every reported rho; present "124 groups" as threshold-conditional; lead with bridge topology. [C, M2]
7. Separate ordinal-gradient from temporal claims; state the dissociation is orientation-invariant. [D]

**Substantive revisions (need author decisions):**
8. Power/sensitivity analysis at the empirical configuration, or reframe the negative as "underdetermined." [A]
9. Report the empirical four-signature correlation matrix. [F]
10. Settlement: rest primacy on mound height, fix or drop the apples-to-oranges mound-area rank-size. [R2-2.2]
11. Sensitivity grid (lat cut × bins × cont × k) for the convergence verdict. [C]
12. Add a hard adversarial mimic (strong assortment without group selection); more seeds with CIs; time-averaging in the validation. [E]
13. Theory: either define β₀/β₁ with an explicit multilevel fitness model and add a group-fitness term, OR retreat to "bounded assortative interaction communities" and drop "unit of selection" (R3 Path A — recommended). [R3]
14. Justify the basin boundary on drainage/phase grounds, show robustness including Salomon/Parchman. [R2-2.1]
15. Downgrade "three independent lines" to one scoped independent line (Cobb) + two data-internal cross-checks. [R2-2.5]
16. Time-averaging sensitivity test for the −0.94 neutrality decline before any unit-of-selection conclusion. [G]
