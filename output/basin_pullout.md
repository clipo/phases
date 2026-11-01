# Does the Parkin phase pull out of the wider LMV set under drift? (n = 55)

All 55 Mainfort-PFG decorated LMV assemblages (8 Parkin-phase, 47 in other phases, after Mainfort 1996); neutral time-transgressive drift on the river network, no boundary imposed, 500 realizations. Focal node: Parkin.

## Pull-out of the Parkin phase
- Mean P(shares Parkin's community): Parkin phase 0.88 (range 0.84-0.91) vs other phases 0.68 (range 0.05-0.87).
- Separation of Parkin from other phases by P: AUC = 0.88 (Mann-Whitney p = 7.6e-04). This separation is what drift on the river network predicts: Parkin-phase assemblages lie near Parkin and so co-occur in its community more than distant assemblages do.
- Per run, Parkin's emergent community is 19% Parkin-phase (precision) and captures 89% of Parkin-phase assemblages (recall); mean 2.1 communities (range 2-3).
- Observed Parkin-vs-others F_ST = 0.012 against a drift null of 0.005 (95% 0.003 to 0.009, max 0.017); 0.4% of drift runs reach or exceed the observed value (observed at the 99.6th percentile of the null).
- Leave-one-out jackknife (matched per-drop null): relabeling single Parkin assemblages into the comparison group moves the observed F_ST to 0.007-0.013; 7 of 8 single drops stay outside their matched drift null, so the contrast is not driven by any one assemblage. Dropping Parkin gives the lowest value, F_ST = 0.007.

Interpretation: with no boundary imposed, drift on the river network already makes the Parkin-phase assemblages share Parkin's emergent community more often than the other-phase (Nodena, Kent, Walls, Tipton, Jones Bayou, Parchman) assemblages do. The pull-out (AUC) is therefore a drift-and-hydrology effect: nearby assemblages co-occur, and the co-membership is graded rather than all-or-nothing. The observed between-group F_ST sits at the extreme upper tail of the drift null and survives most single-site deletions, but it is tiny and the null's own tail reaches past it (the maximum drift value exceeds the observed). Time-averaging, which is not modeled in the null and tends to deflate between-group variance, and the small Parkin sample (n = 8) further limit the claim. We read the data as consistent with drift structured by geography, with at most a weak hint of structure beyond drift at Parkin rather than a demonstrated social boundary.

Figure: figures/fig9_parkin_pullout.png