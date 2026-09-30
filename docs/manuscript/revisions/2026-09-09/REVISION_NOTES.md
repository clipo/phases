# Outcome (2026-09-10, after author review)

The September 9 code corrections listed below were kept. The September 9 and 10 manuscript
drafts were not: the September 10 "constructing social units" reframing was set aside, the
July framing restored, and the drift comparison redone with a calibrated null (analysis 47,
figures from analysis 50). The template build (`49_build_revision.py`) and the
`before/` figure copies were removed from the working tree; the drafts remain in git
history (commit 18abdd8). The headline of the September 9 draft ("the revised baseline
underpredicts differentiation") rested on an uncalibrated innovation rate; the calibrated
result, and its limits, are in the current manuscript and in
`docs/claude_review_of_chatgpt_revision_2026-09-10.md`.

---

# September 9–10, 2026 revision

This revision responds to `docs/code_manuscript_review_2026-09-09.md`. It is a new local draft for author review. The public `phases` repository and Zenodo deposit have not been updated.

**The conclusion changed.** A flat rarefied trajectory is not treated as evidence that bounded groups or a Parkin polity were absent. The revised finite-class baseline underpredicts observed differentiation at its main settings. Some other drift and restricted-copying settings bracket the observations, so the evidence concerns model adequacy and parameter sensitivity rather than a uniquely identified political history.

| Issue | Revision |
|---|---|
| Sampling with replacement described as rarefaction | Added validated without-replacement count sampling and reran empirical/injection analyses |
| Valley accepted by violation score | Corrected fall-to-rise detection; added predicate-agreement and reversal tests |
| Recovery used an ordered-row proxy rather than IDSS | Recovery now computes actual maximal IDSS group counts per assemblage in each bin; the historical proxy experiment is separately labeled |
| Saturating new-class simulation | Added a canonical fixed-measured-class Wright–Fisher model with recurrent innovation; did not represent it as infinite new-class innovation |
| Initialization and unequal initial windows | Added 1,200-generation burn-in, full eight-slice windows, and explicit initialization/repertoire controls |
| Overstated detection/exclusion limit | Power described as sensitivity to the specified increasing-divergence injection; no empirical upper bound inferred from the power curve |
| Selected strong boundary compared with tuned drift | Matched 4 restriction multipliers × 2 lengths × 3 innovation rates in each regional set |
| Wrong interpretation of the Parkin null tail | Reported unusual contrasts as model-data discrepancies, with add-one tail probabilities and chronology sensitivity |
| “Leave-one-out” actually relabeled sites | Removed each assemblage from observed and simulated counts; distinguished observation influence from deleting a latent interaction node |
| Unvalidated cache reuse | Bound caches to numerical source, software versions, inputs, and configuration; added invalidation tests |
| ABC calibration overstated | Disclosed 65% conditional coverage, differing sampler settings, and observation mismatch; excluded archived ABC intervals from the main inferential claim |
| Overstatement about Parkin's political status | Reframed the paper around what ceramic contrasts can and cannot establish |
| CMV comparability | Explicitly identified its 39 assemblages, 24 coordinate locations, and 42 retained classes, including plain wares; treated it as auxiliary rather than an identical decorated-category replication |

**Selected revised results**

- Basin rarefied F_ST rank trend: +0.0225; conditional rarefaction percentiles −0.30 to +0.60.
- Injection 80%-power grid point: s=0.6 without averaging and s=0.7 under the three-bin scenario; independent null false-positive estimates 0.04 and 0.07.
- Basin spatial partition: observed F_ST=0.01794; baseline upper-tail Monte Carlo p=0.00399.
- Named Parkin partition: observed F_ST=0.01210; p=0.02196 under time-transgressive sampling, 0.05988 under contemporaneous time-averaged sampling, and 0.01597 under reversed order.
- Four of eight actual Parkin leave-one-out contrasts remain above their matched 97.5th percentiles. This is an observation-influence check on the same latent network.
- Matched cells bracketing observations: basin 3/6 drift and 6/18 restricted-copying; Parkin 2/6 and 6/18; southeast Missouri 0/6 and 2/18. These are sensitivity counts, not model probabilities.
- Corrected idealized proxy experiment: coupled-profile generator 493/500 hits; each of three specified alternatives 0/500. This is not a universal validation of archaeological group detection.

**Computations completed**

- 1,500 baseline spatial fields, each observed in three sampling configurations.
- 2,160 matched boundary-grid fields.
- 400 repertoire/initialization sensitivity fields.
- 4,000 leave-one-out observation contrasts on the wider-valley baseline fields.
- 2,640 injection test datasets, 800 independent null-calibration draws, and 400 empirical rarefactions, plus a separate 400-draw panel rendering.
- 2,000 idealized generator sequences (500 per generator).

Complete replicate CSVs, summaries, numerical configuration, and environment information are retained in `results/`. A stored seed was independently reproduced in all three regions and all three baseline observation configurations after the reporting changes. The full test run passed 78 tests inside the sandbox; the 13 Bayesian tests blocked by compiler-cache permissions passed when rerun with access to that cache. The additional cache-invalidation test subsequently passed, for 92 tests covered in total. Source whitespace and compilation checks passed. No numerical results were tuned to restore the earlier interpretation.

**Files to review**

- `../../MAIN_TEXT.md`, `../../MAIN_TEXT.docx`, `../../MAIN_TEXT.pdf`
- `../../SUPPLEMENTAL_TEXT.md`, `../../SUPPLEMENTAL_TEXT.docx`, `../../SUPPLEMENTAL_TEXT.pdf`
- `before/`: verbatim previous Markdown sources, bibliography, and copied figure images
- `MAIN_TEXT.template.md` and `SUPPLEMENTAL_TEXT.template.md`: revised prose with explicit result fields
- `result_values.json`: exact values and numerical paragraphs inserted into the draft

Edit the templates when changing prose, then run `analyses/49_build_revision.py`; direct edits to the assembled manuscript will be overwritten by that script. Existing maps, the full-basin deterministic seriation, mound-height ranks, and radiocarbon figure are retained. Their values were not recalculated by the spatial pipeline. The prior pooled-population ABC and hierarchy analyses remain historical exploratory work; a matched Bayesian observation model and repaired calibration are not claimed to be completed in this revision.

The baseline is still a restricted model, including uniform recurrent innovation, equal node populations, a fixed mixing rate, approximate distances, and ordinal occupation positions. The grid is limited and uses only 30 replicates per cell. These limitations are now in the manuscript. A full source-by-source archaeological citation audit and public-release synchronization remain author/submission tasks; this draft is not represented as submission-ready.
