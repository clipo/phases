# Empirical application of the convergence criterion to PFG data (Phase 5)

Neutral pattern report. The four signatures (neutral departure, cultural F_ST, spatial boundary excess, seriation structure) are computed along an ordinal trajectory and at Parkin, plus external settlement/mound cross-checks. No verdict is declared between H1 (co-rise/convergence toward contact, concentrating on Parkin) and H2 (flat / non-trending; Rees 2001). See caveats and 'For interpretation by the team'.

## 1. Data and join

- PFG assemblages loaded: 266 rows -> 263 unique assemblage ids x 26 ceramic types (3 duplicate ids collapsed by summing).
- Matched to LMV with coordinates: 255 / 263 (unmatched: 8); all in UTM Zone 15.
- Parkin (11-N-1) present and matched: True.

## 2. Ordinal axis (data note)

The seriation-frequency file (pfg-cpl-frequency.csv) is the IDSS solution table: 102 numbered ordering positions over only 12 unique decorated assemblages (each repeated, spanning nearly the full range, constant counts). It cannot supply a single per-assemblage ordinal for the 258 sites. We therefore use it for the Parkin-region decorated seriation (Signature 2, below) and adopt the LMV 'Terminal Period' phase (F earliest -> A = contact/latest) as the working ordinal trajectory axis for the 258 sites.
- Decorated assemblages in the seriation file also matched in the 258: 10 (10-P-1, 11-N-1, 11-N-4, 11-N-9, 11-O-10, 13-N-16, 13-N-5, 13-O-10, 13-O-11, 13-O-5).
- Sites with a usable Terminal Period phase: 238 / 255 (others coded '?' or blank).
- Phase counts (ordinal: 1=F earliest ... 6=A contact): 1:4, 2:18, 3:11, 4:62, 5:82, 6:61.

## 3. Spatial clustering (k-means on coordinates)

- Silhouette by k: k=4:0.499, k=5:0.455, k=6:0.425, k=7:0.422, k=8:0.392, k=9:0.380, k=10:0.393.
- Chosen k (max silhouette): **4**.
- Cluster sizes: c0:39, c1:105, c2:31, c3:80.
- Parkin (11-N-1) cluster: c1.

## 4. Ordinal trajectory panel (phase bins x signatures)

Bins (phase ordinal: 1=F ... 6=A), n sites, n clusters represented:
- bin 1: n=4, clusters=1
- bin 2: n=18, clusters=4
- bin 3: n=11, clusters=3
- bin 4: n=62, clusters=4
- bin 5: n=82, clusters=4
- bin 6: n=61, clusters=4

Panel (NaN where a signature is undefined for that bin):

| phase_ord | neutral_departure | fst | spatial_boundary |
|---|---|---|---|
| 1 | 0.7408 | nan | 1.299 |
| 2 | 0.5945 | 0.0131 | 7.792 |
| 3 | 0.1896 | 0.0049 | -19.767 |
| 4 | 0.8431 | 0.0120 | 21.054 |
| 5 | 0.8442 | 0.0133 | 20.652 |
| 6 | 0.6187 | 0.0299 | 22.574 |

Convergence score (z-averaged across the three bin-trajectory signatures; computed on complete bins only):

- phase 2: -0.151
- phase 3: -1.619
- phase 4: 0.426
- phase 5: 0.472
- phase 6: 0.871

Per-signature trend vs phase ordinal. OLS slope (positive = upward toward contact/A), Spearman rank correlation with the ordinal axis (monotone-trend test, robust to scale), and whether the bin series is monotone:

- neutral_departure: slope = +0.02263 (over 6 bins); Spearman rho = +0.314 (p = 0.544); monotone = False
- fst: slope = +0.00420 (over 5 bins); Spearman rho = +0.700 (p = 0.188); monotone = False
- spatial_boundary: slope = +5.30802 (over 6 bins); Spearman rho = +0.771 (p = 0.072); monotone = False
- convergence_score: slope = +0.41359; Spearman rho = +0.900 (p = 0.037)

## 5. Parkin focus and Signature 2 (decorated seriation)

Signature 2 is computed ONLY on the small Parkin-region decorated subset (12 unique assemblages) from the IDSS table; it is NOT a bin-trajectory signature (excluded from the convergence above).
- IDSS at continuity=0.3: 22 maximal co-seriable groups, max solution length 4.
- Parkin (11-N-1) mean ordering position = 57.5 (rank 7/12 among decorated assemblages => **mid**).
- Parkin Signature 2 multi-membership: belongs to **5** co-seriable groups (multi-membership = True; the Lipo 2001 two-lineage / bridge result).
- Robustness of Parkin multi-membership across continuity thresholds:
  - cont=1.0: n_groups=45, Parkin_in_groups=10, multi=True
  - cont=0.5: n_groups=45, Parkin_in_groups=10, multi=True
  - cont=0.3: n_groups=22, Parkin_in_groups=5, multi=True
  - cont=0.1: n_groups=14, Parkin_in_groups=2, multi=True

- Parkin Terminal Period phase = A (working ordinal 6; A=6 is contact/latest).

Per-cluster static Signature 2 was NOT computed on full spatial clusters: each exceeds the ~15-assemblage tractability limit of the deterministic IDSS solver, and the count columns differ (26 PFG types vs the 10 decorated types of the seriation file). Signature 2 is reported only on the decorated Parkin-region subset above.

## 6. Cross-checks (external settlement / mound proxies)

### Rank-size of site size (Max Mound Area, sq ft)
- NOTE: LMV 'Area' is a locality label, not site size; the size proxy used here is 'Max Mound Area (sq ft)'.
- N sites with mound area > 0: 56.
- log-log rank-size slope = -1.028 (Zipf/log-normal expectation ~= -1.0; shallower than -1 = convex below the top sites; steeper = primate).
- Primacy (largest/second) = 1.00; largest site = 20-M-5 ('Swan Lake'); is largest Parkin = False.
- Parkin has no recorded Max Mound Area > 0.

### Mound differentiation
- Max Mound Height (ft): n=255, max=54.0, median=4.0, mean=6.67, 90th pct=15.0.
- Parkin Max Mound Height = 23.0 ft, percentile 94.9, rank 12/255.
- Num_Mounds: n=200, max=15, median=2, sites with >1 mound = 112.
- Parkin Num_Mounds = 4.

### Fortification / monumental features (present count among 255)
- St Francis: 21 sites present (Parkin: 1).
- Ditch: 0 sites present (Parkin: 0.0).
- Platform: 63 sites present (Parkin: 1).

## 7. Figures

- figures/05_ordinal_trajectory.png: the three bin-trajectory signatures (z-standardized) plus the convergence score across the ordinal phase axis (1=F earliest to 6=A contact). Caption: positive trends would indicate co-rise toward contact (H1); flat lines indicate non-trending signals (H2).
- figures/05_spatial_clusters.png: k-means spatial communities in a centered, jittered frame with NO axis scale (coordinates not exposed); Parkin starred. Caption: groups used as units for F_ST and boundary excess.
- figures/05_rank_size_area.png: log-log rank-size of site size (Max Mound Area, sq ft) with OLS slope; Parkin has no recorded mound area so it is absent from this curve (its size shows in Max Mound Height, section 6). Caption: slope near -1 indicates a log-normal/Zipf settlement system; a strongly primate (convex) curve with one dominant site would indicate consolidation.

## 8. Pattern summary (NEUTRAL)

Trend label uses BOTH the OLS slope sign and the Spearman rank correlation. A signature is called 'monotone rise/decline' only when the bin series is monotone AND Spearman |rho| is large; otherwise it is 'upward/downward but non-monotone' (slope sign with volatility) or 'flat' (|rho| small). With only 5-6 bins no Spearman p-value reaches significance, so these are descriptive, not inferential, labels.

- Neutral departure: slope +0.02263, Spearman rho +0.314 -> **upward but non-monotone**.
- Cultural F_ST: slope +0.00420, Spearman rho +0.700 -> **upward but non-monotone**.
- Spatial boundary excess: slope +5.30802, Spearman rho +0.771 -> **upward but non-monotone**.
- Convergence score: slope +0.41359, Spearman rho +0.900 -> **upward (monotone)**.

- Of the 3 bin-trajectory signatures: 3 trend upward toward contact (monotone or not), 0 flat/non-trending, 0 downward. Whether they CONVERGE (co-rise together) is the joint question; the convergence-score slope is +0.41359 (Spearman rho +0.900).
- Parkin: working-axis phase A (ordinal 6); decorated-seriation position rank 7/12; Signature 2 multi-membership in 5 co-seriable groups (bridge / two-lineage status present).
- Settlement: rank-size slope -1.03, primacy 1.00, largest site is 20-M-5.

## 9. Caveats (explicit)

- The trajectory axis is ORDINAL phase (Terminal Period F->A), NOT absolute calendar time; slopes are per-phase, not rates.
- Assemblages are ~2-phase time-averaged (Period spans like B-A, D-C); the working axis uses the single Terminal (latest) phase, which compresses occupation history.
- The named primary seriation file is unusable as a 258-site ordinal (12 unique decorated assemblages, repeated solution positions); the LMV phase axis is used instead. This is a substitution, documented in section 2.
- Signature 2 (seriation) is computed only on the 12-assemblage decorated Parkin-region subset (IDSS is combinatorial) and is EXCLUDED from the bin-trajectory convergence. No per-bin Sig2 trajectory exists.
- The four signatures share the type-frequency substrate and are partly correlated; co-movement is expected to some degree even absent a single causal process.
- The convergence score z-standardizes each signature, so it is scale-free, but the underlying ABSOLUTE movements differ greatly: F_ST stays very low (~0.013 to ~0.030) and neutral departure is non-monotone, while the spatial-boundary signature is large and volatile (e.g. negative at one bin). The high convergence-score rank correlation is driven mainly by the spatial-boundary and F_ST ranks, not by a strong, smooth co-rise of all signatures. Read the panel, not just the score.
- Coverage: 255/263 assemblages matched to coordinates; only 238 carry a usable Terminal Period phase; the seriation file covers 12 unique assemblages.
- Phase bins are uneven (F and D sparse); per-bin signatures from few sites or few clusters are noisy.

## 10. For interpretation by the team

This report does not declare a winner. The pattern above is the evidence; the mapping to hypotheses is the team's call. For reference, the criteria are:
- Consistent with H1 (nascent emergence) IF the bin-trajectory signatures co-rise toward contact (positive convergence-score slope), the rise concentrates spatially on/near Parkin, AND the settlement cross-checks show a primate/convex rank-size with Parkin dominant and differentiated mounds/fortification.
- Consistent with H2 (stable non-consolidation; Rees 2001) IF the signatures are flat / non-trending, the convergence-score slope is near zero, and the settlement system is log-normal (rank-size slope near -1) without a single dominant differentiated center.
- A mixed pattern (some signatures rise, others flat; Parkin a bridge but not a runaway primate) is itself a finding and should be reported as such, not forced to one pole.
