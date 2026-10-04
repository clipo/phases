# Citation Audit: *Explaining the Spatial Structure of Phases: An Example from the Lower Mississippi River Valley* (Lipo, DiNapoli, Madsen)

**Date:** 2026-07-13. **Target:** *American Antiquity*. **Bibliography:** `docs/manuscript/references.bib`. **Prose:** `docs/manuscript/MAIN_TEXT.md`, `docs/manuscript/SUPPLEMENTAL_TEXT.md`.

## Summary

- Total references (bib): 135
- Distinct in-text citation keys: 124
- **Orphan citations (cited in text, no bib entry): 0**
- Orphan references (in bib, never cited): 11
- Spot-check sample: 24
- Verified: 24
- Probably real, minor metadata note: 1 (Wright year, no correction needed)
- Unverified: 0
- Likely fabricated: 0
- **Overall verdict: clean (minor metadata notes only)**

## Bidirectional check

**Forward (citations → references):** every one of the 124 distinct cited keys resolves to a bib entry. Zero orphan citations. This is the strongest integrity result: nothing cited in the prose is missing or invented.

**Reverse (references → citations):** 11 bib entries are never cited (harmless; some are pipeline/method references retained for completeness): `cobb_butler_2002`, `cobb_etal_2024`, `derex_boyd_2016`, `derex_perreault_boyd_2018`, `dunnell_1992`, `ford_1936`, `madgwick_belcher_wolf_2019`, `rees_political_culture` (note: `rees_political_culture` *is* cited in text — see below), `salvatier_wiecki_fonnesbeck_2016`, `saxer_doebeli_travisano_2009`, `williams_1990`. `salvatier_wiecki_fonnesbeck_2016` (PyMC) became uncited when the convergence sentence was cut in the peer-review revision; it can be removed or left as a software reference.

**Ambiguous / near-miss:** none.

## Format consistency

Uniform SAA author-date throughout; author-led citations follow the house rule (cite immediately after the named author). No style drift between main text and supplement. Exactly **one** `[CITE-CHECK]` note remains in the bib, on `dunnell_1992` (editor/volume/pages to confirm) — and that entry is itself an uncited orphan, so no *cited* reference carries an unresolved flag. The manuscript prose contains zero `[CITE-CHECK]`/`[UNVERIFIED]` flags.

## Spot-check findings (24 sampled)

**Verified (metadata matches):** neiman_1995 (Am. Antiq. 60:7-36); lipo_madsen_dunnell_2015 (PLOS ONE 10:e0124942); lipo_dinapoli_madsen_hunt_2021 (PLOS ONE 16:e0250690); bentley_hahn_shennan_2004 (Proc. R. Soc. B 271:1443-1450); henrich_2001 (Am. Anthropol. 103:992-1013); ewens_1972 (Theor. Popul. Biol. 3:87-112); crema_kandler_shennan_2016 (Sci. Rep. 6:39122); crema_edinborough_kerig_shennan_2014 (JAS 50:160-170); toni_etal_2009 (J. R. Soc. Interface 6:187-202); beaumont_zhang_balding_2002 (Genetics 162:2025-2035); reimer_etal_2020 (Radiocarbon 62:725-757); surovell_brantingham_2007 (JAS 34:1868-1877); hunt_bell_travis_2008 (Evolution 62:700-710); perreault_2019 (*The Quality of the Archaeological Record*, Univ. Chicago Press); dunnell_1995 (in Teltser ed., *Evolutionary Archaeology*, Univ. Arizona Press); hopkins_skellam_1954 (Ann. Bot. 18:213-227); phillips_ford_griffin_1951 (Peabody Museum Papers vol. 25); mantel_1967 (Cancer Research 27:209-220); kandler_shennan_2013 (JTB 330:18-25, DOI resolves); feder_kryazhimskiy_plotkin_2014 (Genetics 196:509-522, DOI resolves); mainfort_2001 (in *Societies in Eclipse*, Smithsonian Institution Press); mainfort_1996b; and the two below.

**Notes (all real; on re-check the bib entries are correct):**
- `bell_richerson_mcelreath_2009` — real (PNAS 106:17671-17674, DOI resolves). The bib title is already correct ("...evolution of large-scale human prosociality"). No fix needed; an earlier flag compared against an abbreviated title, not the actual entry.
- `carrignon_bentley_ruck_2019` — real (Palgrave Commun. 5:83, DOI resolves). The bib title already carries the full subtitle ("...Evaluating Neutral Models on Twitter Data with Approximate Bayesian Computation"). No fix needed.
- `wright_1951` — real (Sewall Wright, "The Genetical Structure of Populations," Annals of Eugenics 15:323-354). The publication year is cited **both ways** in the literature: Wiley/Annals of Eugenics dates the article 1949, while much of the genetics literature (and the DOI's volume 15) cites it as 1951. The manuscript's 1951 is a defensible, widely-used choice; no correction required, but note the ambiguity if a reviewer flags it.

## Pattern analysis

No clustering of suspicious citations, no author-topic mismatches, no recent-only confabulation pattern. The unverified rate is 0%. All sampled references are real, established sources; the three notes are wording/date cosmetics, not integrity concerns.

## Recommendation

**Clean.** Reference integrity is sound: zero orphan citations, zero unverified or fabricated references across a 24-reference spot-check, zero flags in the prose, and a single `[CITE-CHECK]` note on one uncited entry (`dunnell_1992`). Remaining pre-submission cleanups have since been completed: the `dunnell_1992` entry was corrected (last `[CITE-CHECK]` cleared) and the uncited `salvatier_wiecki_fonnesbeck_2016` (PyMC3) was removed. On re-check, the `bell_richerson` and `carrignon` titles were already correct in the bib, so no title edits were needed. Nothing outstanding affects the editorial standing of the manuscript.

*Note: one automated verification batch stalled; the six references it covered were verified directly (DOI resolution + web search) in this pass. Tracking updates were skipped — this manuscript is not inside a scaffolded `submissions/<slug>/` tree.*
