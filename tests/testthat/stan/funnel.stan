// Test fixture: Neal's funnel in its centered form, which produces divergent
// transitions under ordinary settings. Used as the probe point at which the
// gate must FAIL with reason "divergences" (rule 3).
parameters { real v; vector[9] x; }
model { v ~ normal(0, 3); x ~ normal(0, exp(v / 2)); }
