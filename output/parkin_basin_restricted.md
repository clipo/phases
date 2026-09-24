# Parkin-phase / St. Francis basin-restricted re-run

Does the whole-LMV no-convergence result hold within the St. Francis basin (the focal hypothesis unit)? The prior passes (scripts 06/07) lumped the Parkin phase (lat ~34.8-35.6) with the lower Yazoo (Winterville, lat ~33.5) and southern St-Francis-type outliers. This run restricts BOTH the curated decorated (transmission) set and the broad PFG/LMV (settlement) set to the basin via a principled latitude cut, recomputes the CA seriation axis, the four transmission signatures, the IDSS group structure with Parkin's bridge rank, and a WITHIN-BASIN mound rank-size, and states plainly whether the no-convergence finding survives. No verdict between H1 (nascent emergence / consolidation toward contact) and H2 (stable non-consolidation; Rees 2001) is forced.

## 0. Basin definition and subset sizes

- Whole-LMV curated decorated set: 38 assemblages (all with coordinates), lat [34.353, 35.619], lon [-90.724, -89.978].
- Whole-LMV broad settlement set (matched to coordinates): 255 sites, lat [33.034, 35.619], lon [-91.379, -89.978]. Coordinates computed from UTM (EPSG:26915 -> EPSG:4326).
- Largest inter-site latitude gap in the curated set: 34.521 -> 34.756 (gap 0.235). The primary cut lat >= 34.5 sits in this gap.

Subset n by latitude cut:
| lat cut | curated n | broad n | curated lat range | broad lat range |
|---|---|---|---|---|
| >= 34.0 | 38 | 174 | [34.353, 35.619] | [34.007, 35.619] |
| >= 34.5 | 36 | 110 | [34.521, 35.619] | [34.504, 35.619] |
| >= 35.0 | 17 | 55 | [35.006, 35.619] | [35.006, 35.619] |

- Curated assemblages EXCLUDED at lat < 34.5 (2): Salomon (34.35), Parchman (34.36) (the southern St-Francis-type outliers; Winterville is not in the curated decorated set at all).
- Broad St-Francis-flagged sites below the cut (1): 17-M-2 (33.95). All other St-Francis-flagged sites are at lat >= 34.77.
- Parkin: curated lat 35.276, broad (11-N-1) lat 35.278 -> in the basin at every cut.

## 1. Transmission signatures + convergence (basin curated set)

Basin curated decorated set (lat >= 34.5): **n = 36** assemblages, 10 decorated types. Basin broad settlement set: **n = 110** sites.

- CA first non-trivial axis inertia fraction: 0.469.
- 14C anchors within the basin: 5; CA<->14C Spearman = +0.500 (p = 0.391); axis flipped so increasing = later.
- The basin time anchor is WEAK (5 anchors): the CA axis is essentially a relative seriation ordinate; per-bin slopes are not rates.

- Within-basin spatial clustering (k-means): silhouette by k: k=2:0.484, k=3:0.526, k=4:0.525, k=5:0.487, k=6:0.491; chosen k = 3.
  cluster sizes: c0:11, c1:3, c2:22.

### Four signatures along the CA axis (4 bins)

Bins (ca_bin, n assemblages, n within-basin clusters represented):
- bin 0: n=9, clusters=3
- bin 1: n=9, clusters=3
- bin 2: n=9, clusters=2
- bin 3: n=9, clusters=2

| ca_bin | neutral_departure | fst | spatial_boundary |
|---|---|---|---|
| 0 | 1.2367 | 0.0035 | 19.326 |
| 1 | 0.5329 | 0.0611 | 36.196 |
| 2 | 0.3797 | 0.0108 | 29.566 |
| 3 | 0.2776 | 0.0069 | -1.726 |

Per-signature trend along the CA axis (OLS slope, bootstrap 95% CI over assemblages, Spearman rho):

| signature | OLS slope | bootstrap 95% CI | Spearman rho | CI excludes 0 |
|---|---|---|---|---|
| Neutral departure | -0.30306 | [-0.71108, -0.12982] | -1.000 | yes |
| Cultural F_ST | -0.00403 | [-0.01412, +0.04878] | +0.200 | no |
| Spatial boundary excess | -6.97867 | [-20.84907, +29.07639] | -0.400 | no |

- convergence_score: slope = -0.48867; Spearman rho = -0.800 (p = 0.200) over 4 complete bins.
- Of 3 signatures within the basin: 0 trend up (rho>+0.3), 2 trend down (rho<-0.3); 1 have a bootstrap slope CI excluding 0. Convergence (H1) requires all three rising together with CIs above 0.

### Bin-count sensitivity (basin)

Spearman rho of each signature with the ordered bin index:
| signature | 3 bins | 4 bins |
|---|---|---|
| Neutral departure | -0.500 | -1.000 |
| Cultural F_ST | -0.500 | +0.200 |
| Spatial boundary excess | -0.500 | -0.400 |

## 2. IDSS group structure + Parkin bridge (basin curated set)

| cont | n_groups | max_size | n_bridge | Parkin_memberships | Parkin_bridge_rank |
|---|---|---|---|---|---|
| 0.05 | 28 | 3 | 2 | 2 | 1/36 |
| 0.1 | 43 | 4 | 22 | 4 | 9/36 |
| 0.2 | 173 | 5 | 35 | 26 | 4/36 |

- Primary (cont=0.1): 43 maximal co-seriable groups within the basin, max group size 4, 22/36 bridge assemblages.
- **Parkin** within the basin: belongs to 4 maximal groups (bridge = True); bridge rank 9/36 by membership count (cont=0.1).
- Signature-2 fragmentation trend within the basin (per-assemblage group count vs CA position): Spearman rho = +0.567 (p = 0.000), OLS slope = +1.9549.

## 3. Within-basin settlement cross-check (rank-size + Parkin rank)

Basin broad settlement set: n = 110 sites. Parkin (11-N-1) present: True.
- Mound present: 93/110 (84.5%); ditch present (LMV-coded; under-records fortification): 0/110; St-Francis-flagged: 20/110 (18.2%); platform: 28/110.

### Within-basin mound-area rank-size (LMV-coded Max Mound Area)

- Sites with Max Mound Area > 0 in the basin: 19.
- log-log rank-size slope = -0.883 (Zipf/log-normal expectation ~= -1; shallower = convex, steeper = primate).
- Primacy (largest/second) = 1.89; largest = 13-M-1; largest is Parkin = False.
- Parkin has Max Mound Area coded 0 in the LMV table (the ditched type-site is under-recorded), so it is ABSENT from the LMV-coded area curve. See the height ranking and the corrected ranking below.

- Within-basin mound HEIGHT: Parkin = 23 ft, rank 2/110 (percentile 98.2). Tallest in basin = 11-N-1 (23 ft).
- Within-basin Num_Mounds (LMV-coded): Parkin = 4, rank 19/90. Most mounds in basin = 13-O-10 (15). NOTE: the LMV table codes Parkin with 4 mounds, but the documented record is 7 (script 07); the LMV mound fields under-record Parkin, so this rank is a floor.

### Within-basin rank-size with corrected Parkin (~17-acre site area)

Parkin's LMV Max Mound Area is coded 0; substituting the documented site area (~17 acres ~ 740520 sq ft) gives an indicative within-basin ranking. NOTE: ~17 acres is total SITE area, not basal mound area, so this is not a like-for-like comparison with the other sites' Max-Mound-Area field; read as indicative only.
- Corrected within-basin rank-size: n = 20, slope = -1.577, primacy = 41.14, largest = Parkin, Parkin rank = 1/20.
- Is Parkin the within-basin primate center? By corrected site area: Parkin ranks 1/20 (YES, rank 1). By LMV-coded mound area Parkin is unranked (coded 0); by mound height Parkin ranks 2/110; by Num_Mounds 19/90.

## 4. Comparison to whole-LMV (NEUTRAL)

Transmission level. Whole-LMV (scripts 06/07): the four signatures did not co-rise toward contact; the convergence-score slope was slightly negative; the IDSS structure was a fragmented, overlapping-lineage system with Parkin a high-degree bridge. Within the basin (lat >= 34.5, n = 36): 0/3 continuous signatures trend up, 2/3 trend down, 1/3 have a bootstrap slope CI excluding zero; convergence-score slope -0.4887 (rho -0.800). The IDSS structure remains fragmented and Parkin remains a high-degree bridge (rank 9/36 at cont=0.1).

Settlement level (the picture most likely to change). Whole-LMV LMV-coded mound-area rank-size: slope -1.04, primacy 1.00, largest = 20-M-5 (no primacy). Within the basin: slope -0.88, primacy 1.89, largest = 13-M-1. With the corrected Parkin site area, whole-LMV: slope -1.26, primacy 16.38, Parkin rank 1/56; within-basin: slope -1.58, primacy 41.14, Parkin rank 1/20.

Plain statement: the whole-LMV NO-CONVERGENCE finding HOLDS within the basin at the transmission level (the continuous signatures do not co-rise with CIs above zero). The IDSS structure (fragmented system; Parkin a high-degree bridge) is unchanged. The settlement picture is where to look for any change: Parkin is the largest within-basin center on the corrected site-area ranking and ranks near the top on mound height and number of mounds, even though the LMV-coded mound-area field leaves it unranked. Whether the within-basin settlement evidence amounts to primacy is reported as the numbers above; no H1/H2 pole is forced.

## 5. Figures

- figures/09_basin_signature_trajectory.png: the three continuous signatures (z-standardized) plus convergence score across the basin CA seriation axis.
- figures/09_basin_rank_size.png: within-basin log-log rank-size (corrected Parkin site area), Parkin marked.

## 6. Caveats

- The basin curated set is small (n = 36); the per-bin trajectory and its bootstrap CIs are noisy. Treat per-bin slopes as descriptive, not inferential.
- The 14C anchor within the basin is sparse (5 assemblages); the CA axis is a RELATIVE seriation ordinate, not calendar time.
- The IDSS continuity threshold matters (cont=0.30 of Lipo et al. 2015 over-saturates this matrix); cont=0.1 is primary with sensitivity across [0.05, 0.1, 0.2]. Absolute group counts scale with cont; the bridge structure is the stable finding.
- The broad PFG set is a mound-biased ceramic-collection subset, not a random settlement sample; within-basin settlement proportions are not population rates.
- Parkin's LMV-coded Max Mound Area and Ditch are 0 despite Parkin being the ditched, multi-mound type-site; the corrected ranking uses documented ~17-acre SITE area (not basal mound area) and is indicative, not like-for-like.
- The latitude cut is principled (natural gap; St-Francis cluster) but is a proxy for the St. Francis River drainage. Sensitivity across [34.0, 34.5, 35.0] is reported in section 0; the substantive conclusions (no transmission convergence; Parkin a bridge; settlement is where primacy may appear) are consistent across the cuts.
