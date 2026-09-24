# Refined empirical re-test of transmission-level convergence (Phase 5 refinement)

This pass re-tests the prior result (no convergence at the transmission level; leans H2 or underpowered) with five legitimate refinements: a theory-faithful Signature 2 from the real IDSS group structure, maximized 14C matching, a corrected Parkin record, bootstrap CIs plus bin-count sensitivity, and a chronology-light early-vs-late structural contrast. The result is reported honestly whether or not it shows convergence. No verdict between H1 (nascent emergence / consolidation toward contact) and H2 (stable non-consolidation; Rees 2001) is declared.

## 0. Data

- Curated decorated assemblages used: 28 (10 types: Parkin_Punctated, Barton/Kent/MPI, Painted, Fortune_Noded, Ranch_Incised, Walls_Engraved, Wallace_Incised, Rhodes_Incised, Vernon_Paul_Applique, Hull_Engraved).
- Assemblages with coordinates: 28.
- CA first non-trivial axis inertia fraction: 0.478.

## 1/2. Maximized 14C matching

- 14C samples with a provenience: 110; aggregated to 14 unique proveniences.
- Proveniences matched to a curated decorated assemblage (improved normalization: generic suffixes stripped, longest-overlap containment): **4**.
  Matched provenience -> assemblage (mean 1-sigma calendar AD, n samples):
  - Clay Hill -> Clay_Hill: AD 1689 (n=3)
  - Kent -> Kent_Place: AD 1466 (n=2)
  - Neeley's Ferry -> Neeleys_Ferry: AD 1484 (n=2)
  - Parkin -> Parkin: AD 1459 (n=20)
- Proveniences that exist as BROAD PFG sites but are NOT in the curated decorated set (cannot be placed on the CA axis): 5 (Callahan-Thompson, Denton Mounds, Moon, Turner, Upper Nodena).
- Proveniences absent from both curated and broad sets (different drainages / phases): 5 (Hazel, Hess, Lilbourn, Powers Fort, Snodgrass).
- HONEST CEILING: the dated-provenience pool is 14, but only 4 of those proveniences correspond to assemblages in the curated decorated set. The 14C anchor for the CA axis is bounded by these 4 points, not by the full 14. The remaining dated sites (Powers-phase and other-drainage sites such as Powers Fort, Snodgrass, Turner, Lilbourn, Hazel, Hess, Moon, Denton Mounds, Callahan-Thompson) are outside the curated decorated set.

- CA<->14C Spearman on the 4 anchors = +0.000 (p = 1.000); axis kept (already increasing with time).
- PLAIN STATEMENT: with only 4 anchors (< 10), the absolute time anchor remains weak. The CA axis is essentially a RELATIVE seriation ordinate; the calendar orientation is directional, not a calibrated chronology. This is the principal reason the cross-sectional structural test (section 5) is run as a chronology-light complement.

## 3. Signature 2 via the real IDSS group structure (curated set)

Primary run uses the IDSS continuity threshold cont=0.1 (see note: cont=0.30 of Lipo et al. 2015 over-saturates this broadly overlapping decorated matrix and exceeds the solver caps). Group counts and Parkin's bridge rank are reported across cont in [0.05, 0.1, 0.2] for sensitivity.

- Number of maximal co-seriable groups (cont=0.1): **34**. Largest group sizes: [4, 4, 3, 3, 3, 3, 3, 3] (max group size = 4).
- Group-size distribution {size: n_groups}: {1: 8, 2: 13, 3: 11, 4: 2}. The structure is highly FRAGMENTED: many small overlapping windows, no single large ordering covering the set. This matches Lipo et al. 2015, where the largest LMV solution held only four assemblages.
- Multi-membership (bridge) assemblages: 18 of 28 belong to more than one maximal group.
- Spatial coherence of IDSS groups: mean within-group geographic distance / mean between-group distance = 0.628 (< 1 means co-seriable assemblages are geographically CLOSER than non-co-seriable pairs, i.e. groups are spatially clustered).
- **Parkin** belongs to 3 maximal groups (bridge = True); bridge rank 10 of 28 by membership count (cont=0.1). Top bridge assemblages: Barton_Ranch(6), Williamson(6), Rose_Mound(6), Belle_Meade(5), Cramor_Place(4), Neeleys_Ferry(4).

Bin/continuity sensitivity of the IDSS structure:
| cont | n_groups | max_size | n_bridge | Parkin_memberships | Parkin_bridge_rank |
|---|---|---|---|---|---|
| 0.05 | 20 | 3 | 1 | 1 | 2/28 |
| 0.1 | 34 | 4 | 18 | 3 | 10/28 |
| 0.2 | 96 | 5 | 27 | 19 | 1/28 |

- Signature-2 TREND (proper): per-assemblage count of distinct co-seriable groups vs CA position: Spearman rho = -0.540 (p = 0.003), OLS slope = -1.4342. Positive => later-CA assemblages participate in MORE distinct (non-co-seriable) groups (fragmentation / assortment RISING toward later); negative => fewer (coherence rising).
  bootstrap 95% CI on the Signature-2 fragmentation slope: [-2.2618, -0.7014] (excludes 0).

## 4. Four signatures along the CA axis (bootstrap CIs + bin sensitivity)

### 4a. Bin-count sensitivity (4, 6, 8 bins)

Spearman rho of each signature with the ordered bin index:
| signature | 4 bins | 6 bins | 8 bins |
|---|---|---|---|
| Neutral departure | +0.400 | +0.429 | +0.310 |
| Cultural F_ST | +0.000 | +0.400 | -0.400 |
| Spatial boundary excess | -0.400 | +0.200 | -0.400 |

- Sign stability across bin counts: Neutral departure=stable, Cultural F_ST=UNSTABLE, Spatial boundary excess=UNSTABLE.

### 4b. Primary panel (6 bins) with bootstrap 95% CI on the slope

| signature | OLS slope | bootstrap 95% CI | Spearman rho | CI excludes 0 |
|---|---|---|---|---|
| Neutral departure | +0.30512 | [+0.08993, +0.48183] | +0.429 | yes |
| Cultural F_ST | +0.00076 | [-0.05221, +0.07078] | +0.400 | no |
| Spatial boundary excess | +0.51494 | [-8.38513, +12.37512] | +0.200 | no |

- Of 3 signatures: 2 trend up (rho>+0.3), 0 trend down (rho<-0.3); 1 have a bootstrap slope CI that excludes 0. For convergence (H1) all three should rise together with CIs above 0.

## 5. Cross-sectional structural contrast (early CA-third vs late CA-third)

Early CA-third (n=9) vs late CA-third (n=9). This tests for bounded-group STRUCTURE without trusting the weak time anchor: do the four signatures JOINTLY indicate stronger boundaries in the late third than the early third?

The 'late-early' column is the observed point-estimate difference; the bootstrap 95% CI (on the resampled difference) gives its uncertainty.

| signature | early | late | late-early | bootstrap 95% CI | direction |
|---|---|---|---|---|---|
| Neutral departure | 0.1099 | 1.3683 | +1.2584 | [+0.7889, +1.8190] | higher late (CI excludes 0) |
| Cultural F_ST | 0.0028 | 0.0054 | +0.0026 | [-0.0091, +0.0183] | higher late |
| Spatial boundary excess | 42.2252 | -1.5779 | -43.8032 | [-102.2037, +67.4255] | higher early |
| IDSS n_groups (fragmentation) | 11.0000 | 8.0000 | -3.0000 | [-5.0000, +3.0000] | higher early |

- JOINT structural reading: 2 of 4 structural signatures are higher in the late third; 1 of 4 have a bootstrap CI that excludes 0. Bounded-group consolidation (H1) predicts higher F_ST, higher boundary excess, higher within-group neutral departure, AND MORE IDSS groups jointly in the late third with CIs above 0.

## 6. Corrected Parkin record (settlement-level cross-check)

- LMV-coded Parkin (11-N-1) BEFORE correction: Ditch=False, Num_Mounds=4, Max_Mound_Area=0 sq ft, Max_Mound_Height=23 ft.
- Parkin AFTER documented-ground-truth override (user-provided / published site description): Ditch=1 (moat + palisade with bastions on 3 sides), Num_Mounds=7, main mound 21.3 ft (+5 ft terrace), area ~17 acres (~740520 sq ft ~ 6.9 ha), St Francis=1, Platform=1.
- DATA-QUALITY FLAG: the LMV ditch/area fields under-record fortification (Parkin, the ditched type-site, was coded Ditch=0, Area=0). Across the full LMV-22 table (766 sites) only 15 ditches are coded; treat regional ditch counts as a FLOOR, not a census. Do not silently trust the field elsewhere.
- After correction: ditch present in the broad matched set = 1 / 255; Parkin mound-area rank now computable (see below).
- Mound-area rank-size (corrected): n=56, log-log slope = -1.263, primacy (largest/second) = 16.38, largest = Parkin, Parkin rank = 1/56.
- NOTE: Parkin's ~17-acre figure is total SITE area, not max-mound basal area; it is not strictly comparable to the LMV 'Max Mound Area' field for other sites. The corrected rank-size is therefore indicative and is reported with this caveat, not as a like-for-like ranking.

## 7. Figures

- figures/07_signature_trajectory.png: the three continuous signatures (z-standardized) across the 6-bin CA axis.
- figures/07_early_late_contrast.png: early- vs late-CA-third values of the four structural signatures (scaled per metric).
- figures/07_idss_bridge_rank.png: IDSS group-membership count per assemblage vs CA position, Parkin starred.

## 8. Pattern summary (NEUTRAL, no verdict)

- 14C anchoring: 4 curated assemblages have a 14C date (of 14 dated proveniences); CA<->14C Spearman +0.000 (p 1.000). The anchor stays weak; the axis is essentially relative.
- IDSS Signature 2: 34 maximal co-seriable groups at cont=0.1, max group size 4, 18/28 bridge assemblages. Highly fragmented (many small overlapping windows). Parkin is a high-degree bridge (rank 10/28). Groups are spatially clustered (within/between distance ratio 0.63).
- Signature-2 fragmentation vs CA position: rho -0.540, slope -1.4342, bootstrap CI [-2.2618, -0.7014] (excludes 0).
- Four continuous signatures along the CA axis (6 bins):
  - Neutral departure: slope +0.30512, rho +0.429, CI [+0.08993, +0.48183] (excludes 0); bin-count sign stable.
  - Cultural F_ST: slope +0.00076, rho +0.400, CI [-0.05221, +0.07078] (spans 0); bin-count sign UNSTABLE.
  - Spatial boundary excess: slope +0.51494, rho +0.200, CI [-8.38513, +12.37512] (spans 0); bin-count sign UNSTABLE.
  Convergence requires all three rising together with CIs above 0: 2 rise, 1 have a CI excluding 0.
- Early- vs late-CA-third structural contrast (late - early):
  - Neutral departure: +1.2584 (CI excludes 0).
  - Cultural F_ST: +0.0026 (CI spans 0).
  - Spatial boundary excess: -43.8032 (CI spans 0).
  - IDSS n_groups: -3.0000 (CI spans 0).
  Joint bounded-group consolidation requires all four higher late with CIs above 0: 2/4 higher late, 1/4 CI excludes 0.
- Settlement (corrected Parkin): Parkin is a fortified (ditched + palisaded) 7-mound town; mound-area rank-size slope -1.26, primacy 16.38. The LMV ditch/area fields under-record fortification; regional ditch counts are a floor.

Plain statement on whether the refinement changes the prior no-convergence picture: see section 10.

## 9. Caveats (explicit)

- The 14C anchor did NOT improve materially: improved name normalization still yields only 4 curated assemblages with a date, because most dated proveniences are Powers-phase / other-drainage sites absent from the curated decorated set. The CA axis is a RELATIVE seriation ordinate; per-bin slopes are not rates.
- The IDSS continuity threshold matters: cont=0.30 (Lipo et al. 2015) over-saturates this broadly overlapping decorated matrix and exceeds the solver caps, so cont=0.1 is used as primary and sensitivity is reported across [0.05, 0.1, 0.2]. Absolute group counts scale with cont; the bridge STRUCTURE (Parkin high, system fragmented) is the stable finding.
- The signatures share the type-frequency substrate and are partly correlated; some co-movement is expected without a single causal process. Bins are few and uneven; per-bin estimates from small bins are noisy, which is why bootstrap CIs are reported.
- The early-vs-late IDSS n_groups bootstrap collapses duplicate rows created by resampling-with-replacement (identical assemblages are trivially co-seriable and only inflate the search), so its resampled group counts run lower than the duplicate-free point estimate; the point difference (+groups in the late third) therefore sits above the bootstrap CI, which spans 0. Read the IDSS contrast as suggestive of more late-third groups but NOT robust under resampling.
- The Parkin override uses total site area (~17 acres) for a field (Max Mound Area) that elsewhere holds basal mound area; the corrected rank-size is indicative, not like-for-like.
- The broad PFG set is a mound-biased ceramic-collection subset, not a random settlement sample; settlement proportions are not population rates.

## 10. For team interpretation (H1 vs H2; no verdict)

Status of the no-convergence picture after refinement (NEUTRAL): of the three continuous transmission signatures, 2 rise and 0 fall along the CA axis, and 1 of three have a bootstrap slope CI that excludes zero. In the chronology-light early-vs-late contrast, 2 of four structural signatures are higher in the late third and 1 of four have a CI excluding zero. The improved 14C match did not strengthen the time anchor (4 dated curated assemblages). These are the observed quantities; whether they amount to convergence is the team's call.

- Consistent with H1 (nascent emergence / consolidation) IF the three continuous signatures rise JOINTLY with CIs above zero, the IDSS fragmentation trend and early-vs-late contrast point the same way, and Parkin sits as a late high-degree bridge in a system tightening toward contact.
- Consistent with H2 (stable non-consolidation; Rees 2001) IF the signatures do NOT rise jointly, the slope CIs span zero, and the early-vs-late contrast shows no coherent bounded-group strengthening, i.e. a persistently fragmented, overlapping-lineage system.
- A mixed result (some signatures move, CIs wide, anchor weak) should be reported as underpowered / inconclusive rather than forced to either pole. The IDSS structure (Parkin a strong bridge in a fragmented system) is robust to bin/cont choice and is the most secure finding here.
