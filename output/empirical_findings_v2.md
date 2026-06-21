# Two-level empirical application of the convergence criterion (Phase 5 v2)

Neutral pattern report on two datasets. TRANSMISSION level: a correspondence-analysis (CA) seriation axis on the curated decorated assemblages (~55), anchored and oriented by 14C calendar dates, with the four signatures computed along it. SETTLEMENT level: mound and ditch presence and a mound-area rank-size on the broad PFG/LMV set (~258). Parkin links the two. No verdict is declared between H1 (co-rise/convergence toward contact, concentrating on Parkin) and H2 (flat / non-trending; Rees 2001). See caveats and the closing section.

## Transmission level (curated decorated set)

### 1. Data and join

- Curated decorated assemblages: 55 rows, 10 decorated types (Parkin_Punctated, Barton/Kent/MPI, Painted, Fortune_Noded, Ranch_Incised, Walls_Engraved, Wallace_Incised, Rhodes_Incised, Vernon_Paul_Applique, Hull_Engraved).
- Assemblages used for CA: 55.
- Curated assemblages joined to coordinates (mainfort-pfg-cplXY.txt): 55 / 55.
- 14C samples: 110 with a provenience; 110 parsed to a 1-sigma calendar midpoint; aggregated to 14 proveniences (mean per provenience).
- 14C proveniences matched to a curated assemblage: 5 (of 14).
  Provenience -> assemblage (mean 1-sigma calendar AD, n samples):
  - Clay Hill -> Clay_Hill: AD 1689 (n=3)
  - Kent -> Kent_Place: AD 1466 (n=2)
  - Neeley's Ferry -> Neeleys_Ferry: AD 1484 (n=2)
  - Parkin -> Parkin: AD 1459 (n=20)
  - Upper Nodena -> Upper_Nodena: AD 1354 (n=1)

### 2. CA seriation axis and 14C anchor

- CA first non-trivial axis carries 0.479 of total inertia (the seriation ordinate; larger = later after orientation).
- Orientation: Spearman(CA ordinate, 14C mean date) on the 5 dated assemblages = +0.500 (p = 0.391); axis kept (already increasing with time).
- Coarse linear CA->calendar map (for reporting only): AD = 1516 + 53.9 * CA_ordinate (n=5 anchors).

### 3. Spatial clustering of curated assemblages (k-means)

- Silhouette by k: k=2:0.545, k=3:0.502, k=4:0.540, k=5:0.577, k=6:0.546; chosen k = 5.
- Cluster sizes: c0:13, c1:11, c2:14, c3:4, c4:13.
- Parkin cluster: c1.

### 4. Four signatures along the CA axis

Curated assemblages with coordinates binned into 6 ordered windows along the CA ordinate (bin 0 = earliest).

Bins (ca_bin, n assemblages, n clusters, CA range):
- bin 0: n=10, clusters=4, CA [-1.785, -0.842]
- bin 1: n=9, clusters=3, CA [-0.816, -0.560]
- bin 2: n=9, clusters=5, CA [-0.543, -0.133]
- bin 3: n=9, clusters=4, CA [-0.112, +0.043]
- bin 4: n=9, clusters=4, CA [+0.088, +0.211]
- bin 5: n=9, clusters=2, CA [+0.256, +0.404]

| ca_bin | neutral_departure | fst | spatial_boundary | sig2_n_groups |
|---|---|---|---|---|
| 0 | 1.1360 | 0.0606 | -2.720 | 18 |
| 1 | 1.0434 | 0.0235 | -4.198 | 19 |
| 2 | 0.6268 | 0.1306 | -8.113 | 11 |
| 3 | 0.1968 | 0.0109 | -2.786 | 26 |
| 4 | 0.2399 | 0.0189 | 1.457 | 32 |
| 5 | 0.3047 | 0.0140 | 7.811 | 32 |

Signature 2 (IDSS n_groups) is reported per bin where the bin is tractable (3-14 assemblages); it is NOT entered into the convergence score (it is a count of co-seriable groups, not a continuous magnitude, and bins above 14 are intractable).

Per-signature trend along the CA axis. OLS slope per bin (positive = rising toward later/contact), Spearman rank correlation with the ordered bin index, and monotonicity:

- neutral_departure: slope = -0.19992 (over 6 bins); Spearman rho = -0.771 (p = 0.072); monotone = False
- fst: slope = -0.01047 (over 6 bins); Spearman rho = -0.657 (p = 0.156); monotone = False
- spatial_boundary: slope = +2.14117 (over 6 bins); Spearman rho = +0.600 (p = 0.208); monotone = False
- convergence_score: slope = -0.11490; Spearman rho = -0.486 (p = 0.329) (over 6 complete bins).

### 5. Parkin link (transmission level)

- Parkin CA ordinate = +0.210, rank 45/55 along the (oriented) axis => **late**.
- Parkin mapped to calendar (coarse): AD 1527.
- Parkin 14C mean date: AD 1459.
- Parkin neighborhood (spatial cluster c1): 11 assemblages.
- Signature 2 in Parkin's neighborhood: 20 co-seriable groups; Parkin belongs to 7 (multi-membership / bridge = True).

## Settlement level (broad PFG/LMV set)

### 6. Settlement features

- Broad PFG assemblages: 266 rows -> 263 unique ids; matched to LMV coordinates: 255; of those, 255 joined to LMVData-22March2006 binary features by Number.
- Parkin (11-N-1) present in broad matched set: True.
- Mound present: 230 / 255 (90.2%).
- Ditch present (defensive): 0 / 255 (0.0%).
- DATA NOTE: across the full LMVData-22March table (766 sites) ditches number 15; the PFG-matched subset captures 0 of them. The PFG sites are a ceramic-collection subset, not a random settlement sample, and are mound-biased (see mound %).
- St-Francis (fortified) present: 21 / 255.
- Platform present: 63 / 255.
- Parkin features (LMVData-22March): Mound=True, Ditch=False, St-Francis=True, Num_Mounds=4, Max_Mound_Height=23 ft, Max_Mound_Area=0 sq ft.
- DATA NOTE: the LMV coding records Parkin Ditch=0 and Max Mound Area=0, even though the Parkin site is the type-site of a ditched/palisaded fortified town. The LMV ditch/area fields do not capture Parkin's known defensive ditch; this is a data-coding limitation, flagged for the team.

### 7. Mound-area rank-size

- Sites with Max Mound Area > 0: 55.
- log-log rank-size slope = -1.038 (Zipf/log-normal expectation ~= -1; shallower = convex below the top, steeper = primate).
- Primacy (largest/second) = 1.00; largest site = 20-M-5; largest is Parkin = False.
- Parkin has no recorded Max Mound Area > 0 (excluded from this curve; see height below).
- Parkin Max Mound Height = 23 ft, percentile 94.9, rank 12/255 (n with recorded height).

### 8. Spatial concentration near Parkin

- Broad k-means: chosen k=4 (max silhouette); Parkin cluster c1, size 105.
- Mound presence: 83.8% inside Parkin's cluster vs 94.7% outside.
- Ditch presence: 0.0% inside Parkin's cluster vs 0.0% outside (n ditched inside = 0, outside = 0).

## 9. Figures

- figures/06_ca_trajectory.png: the four signatures (z-standardized) plus convergence score across the CA seriation axis.
- figures/06_ca_vs_14c.png: CA ordinate vs mean 14C calendar date for the dated curated assemblages, with the orientation fit.
- figures/06_mound_rank_size.png: log-log mound-area rank-size with OLS slope; Parkin marked (Parkin has area=0 so it is absent unless recorded).
- figures/06_mound_ditch_map.png: relative jittered settlement map (no axis scale) of mound/ditch presence with Parkin starred.

## 10. Pattern summary (NEUTRAL)

Trend labels use both the OLS slope sign and the Spearman rank correlation; with only a handful of bins these are descriptive, not inferential.

- Neutral departure: slope -0.19992, Spearman rho -0.771 -> **downward but non-monotone**.
- Cultural F_ST: slope -0.01047, Spearman rho -0.657 -> **downward but non-monotone**.
- Spatial boundary excess: slope +2.14117, Spearman rho +0.600 -> **upward but non-monotone**.
- Convergence score: slope -0.11490, Spearman rho -0.486 -> **downward (non-monotone)**.

- Of the 3 continuous signatures: 1 trend upward toward later, 0 flat/non-trending, 2 downward. Convergence-score slope -0.11490 (rho -0.486).
- CA<->14C: Spearman +0.500 on 5 anchors (p 0.391); CA axis already time-ordered.
- Parkin (transmission): CA position rank 45/55 (late); 14C AD 1459; Signature-2 bridge status as reported in section 5.
- Parkin (settlement): mound present True, ditch coded False (see data note), height 23 ft.
- Settlement system: mound % = 90.2, ditch % = 0.0, mound-area rank-size slope = -1.04, primacy = 1.00, largest = 20-M-5.

Plain statement: report the panel, not a single index. Whether the transmission signatures co-rise toward contact, whether Parkin is an early/mid/late bridge, and whether the settlement system is primate or log-normal are the observed quantities above. No pole is asserted.

## 11. Caveats (explicit)

- The CA axis is RELATIVE. Even oriented and coarsely calibrated by 5 14C anchors, it is a seriation ordinate, not calendar time; per-bin slopes are not rates.
- 14C anchoring is sparse: only 5 proveniences match curated assemblages, and the name-matching used exact-then-substring matching (e.g. Kent -> Kent_Place); the CA<->date correlation rests on these few points.
- Date parsing approximates each sample by the mean of the first and last integers in the 1-sigma calibrated string (the 1-sigma midpoint); multi-intercept samples are summarized by this midpoint, and per-provenience dates are simple means over samples.
- The two levels use different datasets and ID conventions: ~55 curated decorated assemblages (name keys) vs ~258 broad PFG sites (grid-id keys). Parkin is the only guaranteed cross-level link.
- The broad PFG set is a ceramic-collection subset, mound-biased (90% have a mound) and capturing only 0 of 15 ditched sites in the full LMV table; settlement proportions are NOT a random settlement sample.
- The LMV coding records Parkin Ditch=0 and Mound-Area=0 despite Parkin being a known ditched town with mounds; the LMV ditch/area fields under-record fortification. Treat the ditch counts as a floor.
- The four signatures share the type-frequency substrate and are partly correlated; some co-movement is expected even without a single causal process.
- Signature 2 (IDSS) is a count of co-seriable groups, computed only where a bin or neighborhood is tractable (3-14 assemblages); it is shown for context but excluded from the convergence score.
- Bins are few (6) and uneven; per-bin signatures from small bins or few clusters are noisy.

## 12. For team interpretation (no verdict)

This report does not declare a winner. The numbers above are the evidence; mapping to hypotheses is the team's call. For reference:
- Consistent with H1 (nascent emergence) IF the transmission signatures co-rise toward the later end of the CA axis (positive convergence-score slope), the rise concentrates on/near Parkin, and the settlement system shows a primate/convex mound-area rank-size with Parkin dominant and differentiated.
- Consistent with H2 (stable non-consolidation; Rees 2001) IF the signatures are flat / non-trending along the CA axis, the convergence-score slope is near zero, and the settlement system is log-normal (rank-size slope near -1) without a single dominant differentiated center.
- A mixed pattern (some signatures rise, others flat; Parkin a bridge but not a runaway primate) is itself a finding and should be reported as such, not forced to one pole.
