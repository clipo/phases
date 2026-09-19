# Figure Plan — mls-emergence

Target: *American Antiquity* (figure-economical; aim ~6-7 main-text figures, rest to Supplemental). House style (global CLAUDE.md): Okabe-Ito or viridis palettes, Arial/sans-serif, PNG 300 dpi, width 7 in (3.5 in single-column), no figure titles (description in captions), minimal/clean, no chartjunk. **Data sensitivity:** maps use regional scale / generalized or jittered site locations; never publish precise UTM/lat-long for protected sites (these are USGS-derived and the repo is private, but published maps should still generalize).

Status key: [READY] generatable from an existing analysis script; [SCRIPT] needs a small new plotting script; [DESIGN] conceptual diagram to be drawn (matplotlib/tikz/inkscape).

## Main text

Figures are numbered in document order. The signature-by-mechanism criterion matrix is **Table 1** (not a figure): rows = the four signatures, columns = genuine emergence + three mimics, cells shaded for which signatures each mechanism lights up.

**F1. Study area and settlement (map).** [DONE — `analyses/make_map.py` → `fig1_studyarea.png`] Central Mississippi Valley / St. Francis basin on a geographic base (LMV hydrology, rivers, geology, state/county outlines); assemblage locations, Parkin starred, markers sized by mound height, Parkin-phase delineation (buffer-union of lat≥34.5 St-Francis-type sites; Winterville excluded), North America locator inset. Orientation + settlement context for H1/H2.

**F2. The overall model (conceptual).** [DONE — `analyses/make_figures.py::fig2_model` → `fig2_model.png`] Promoted from the former S1. The shared-latent-assortment bridge (one latent assortment governs both monument signaling [parent model] and ceramic-style transmission [this study]); the two regimes along an assortment axis with the φ* threshold (β0 individually-rewarded signaling, aggregated/porous, Parkin a prominent node, vs β1 group-level emergence, bounded groups); and the four-signature readout (dissociated below φ*, convergent above). The paper's model in one figure.

**F3. Seriation / CA ordination of the decorated assemblages.** [DONE — `make_figures.py::fig3_ca_ordination` → `fig3_ca_ordination.png`] The CA seriation ordination of the curated decorated set (assemblages in CA space), colored by CA1 rank, Parkin starred; shows the **gradational continuum** (no discrete clusters), echoing Mainfort 2003 NMDS. The seriation-axis basis for the empirical test.

**F4. The criterion validates (simulation).** [DONE — `make_figures.py::fig4_validation` → `fig4_validation.png`] Small-multiples: the four-signature ordinal trajectories under each simulated mechanism, showing co-rise only for genuine emergence; inset with the independence-audit correlation (mean |r| ≈ 0.88 genuine vs ≈ 0.2 mimic). Establishes the instrument works and convergence is non-trivial.

**F5. Empirical signatures along the seriation axis (the result).** [DONE — `make_figures.py::fig5_empirical_trajectory` → `fig5_empirical_trajectory.png`] The four-signature ordinal trajectory on the real data with bootstrap CIs: **no convergence** — seriation fragmentation rises while neutral-departure and F_ST decline. The honest empirical finding.

**F6. IDSS group structure and the Parkin bridge.** [DONE — `make_figures.py::fig6_idss_structure` → `fig6_idss_structure.png`] Two-panel: histogram of co-seriable group sizes (highly fragmented; max size ~4) + lollipop of top assemblages by group-membership count, with **Parkin highlighted as a high-degree bridge node**. The headline structural result and the Lipo 2001 hook, reproduced.

**F7. Settlement/economic cross-check.** [DONE — `make_figures.py::fig7_ranksize` → `fig7_ranksize.png`] Rank-size of mound area (near-log-normal, slope ≈ -1.0, no primate center; Parkin marked), with the mound-presence / ditch / St-Francis proportions. The verticality/economic-centralization leg (pairs with the cited Cobb lithic evidence).

**F8 (optional, may merge into a table). Convergence of evidence (synthesis).** [DESIGN] A compact synthesis graphic/table: the four independent lines (ceramic convergence test; settlement rank-size; lithic exchange [Cobb]; Mainfort 2003 NMDS continuum) all → non-consolidation; Parkin = ceremonial/bridge node, β0 signaling regime not β1 fitness reorganization. Could be Table 1 instead of a figure.

## Supplemental

**S1. [PROMOTED to main-text F2.]** The two-regime / shared-latent-assortment schematic now appears up front as Figure 2 (after the map), per its role in framing the whole argument.

**S2. monument-mls bistable emergence.** [READY from `signaling.emergence`] φ* saddle and replicator trajectories (above φ* → 1, below → 0); the parent model's emergence dynamics.

**S3. Coupling robustness.** [READY from `analyses/08_coupling_robustness.py`] Convergence-detectability vs coupling strength: the criterion fires for coupling ≥ ~0.1 (flat plateau above), so the empirical non-detection is informative. The load-bearing-assumption sensitivity.

**S4. Full validation panels.** [READY from `04`] All four mechanisms × four signatures, full panels + discrimination table + cross-seed robustness.

**S5. Chronology.** [READY from `06`] CA ordinate vs the 5 Mainfort ¹⁴C anchors (the relative-axis caveat), and the phase (F→A) overlay. Honest display of the weak absolute anchor.

## Generation notes
- Figures from `04/05/06/07/08` already run; need a `figures/` export pass in house style (currently exploratory). Add a `analyses/make_figures.py` that regenerates all main+supp figures deterministically for the manuscript.
- F1, F8, S1 are conceptual diagrams (no data) — draw in matplotlib or vector tool.
- F3 (map): build from the assemblage coordinates (gitignored); export at regional scale with generalized locations. Consider a base layer (rivers/state outlines) without exposing exact site points (jitter or symbol-only).
- Caption convention: full descriptive caption (no on-figure title), state n, source, and method.
