"""Per-fit sampler diagnostics for every Stan fit of the latent-group probe (rule 16).

The fit records that R/fit_spatial.R writes beside each cached fit in
output/gp_cache/ (gitignored, 29 GB with the .rds fits) hold every numeric
diagnostic the gate read: R-hat, bulk and tail ESS with the worst variable,
divergences, treedepth saturation against the configured maximum, E-BFMI
per chain, the per-chain co-membership split for A1, the sampler rung and
its settings, and the ladder of earlier rungs. The committed findings
tables (a1_validation.md, a3_validation.md, a3_pinned_rho6.md,
a3_pinned_rho12.md, a1_prior_check.md) print the gate outcome and, for most
fits, only its reasons. This script copies the numbers out of the records
into one committed table so the supplement's promise of diagnostics for
every fit is kept by a file in the repository rather than by the cache.

One row per record. A label can have two records when the same fit was
run under two versions of the Stan source: the A3 and A1 families were
first fitted with the class-intercept prior hard-coded as Normal(0, 2) and
rerun after D25 (2026-10-18) made the prior data with default Normal(0,
10). The source version is told from the record's Stan source hash; for a
data-driven version the prior SD is read from the data file the record
names (every such file was written on 2026-10-18 with prior_alpha_sd = 10,
after D25, so it cannot predate the fit that read it).

Descriptive export; no inference. Exploratory design validation, not
evidence about the record.
Usage: .venv/bin/python analyses/probe_latent_groups/export_fit_diagnostics.py
Output: output/findings/fit_diagnostics.md, output/findings/fit_diagnostics.csv
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "output" / "gp_cache"
OUT_MD = ROOT / "output" / "findings" / "fit_diagnostics.md"
OUT_CSV = ROOT / "output" / "findings" / "fit_diagnostics.csv"

# Stan source versions, by the sha256 that R/fit_spatial.R::stan_source_hash
# computes over the model file and its includes. "hard-coded" versions carry
# `alpha ~ normal(0, 2)` in the model block; "data" versions read
# prior_alpha_mean and prior_alpha_sd from the data (D25, commit 535009b).
SOURCE_VERSIONS = {
    "83410f2962dec4d4ed0fad6dc255465fbea9568b7c9bedab683419ebb3b4364b": ("composition_gp_marginal.stan", "hard-coded"),
    "35fbf97122c9ab244be792a09756b4997094defe3fa1ea181890c703dae5b902": ("composition_gp_marginal.stan", "data"),
    "36f84cbe09100141a481f2f0ea3209e21209f059f668af9b50bff4a59ec365ba": ("composition_gp_groups.stan", "hard-coded"),
    "67327d0cbb907eb61c85bb435b4ca1316c8df7d2dbb4cbd7cb2c91d330395fe9": ("composition_gp_groups.stan", "data"),
    "2db9f963": ("composition_gp_barrier.stan", "hard-coded"),
    "feeb815ee5a75f19cca8a93f90cc7c4f1ffea413dcb07553a45f8400da1a5c25": ("composition_gp_barrier.stan", "data"),
}
FAMILIES = [
    ("a0_", "A0 spatial field: basin fit, prior predictive, sensitivity, recovery, intercept priors (refit_a0_basin.R, a0_intercept_prior.R)"),
    ("a1val_", "A1 latent-unit recovery and equivalence fits (validate_a1.R; a1_validation.md)"),
    ("a1prior_", "A1 refits under the Normal(0, 10) intercept prior (a1_prior_check.R; a1_prior_check.md)"),
    ("a3val_", "A3 barrier recovery and equivalence fits, range free (validate_a3.R; a3_validation.md and a3_validation_alphaN0-2.md)"),
    ("a3pin6_", "A3 with the range pinned at 6 km (a3_pinned_range.R 6; a3_pinned_rho6.md and _alphaN0-2.md)"),
    ("a3pin12_", "A3 with the range pinned at 12 km (a3_pinned_range.R 12; a3_pinned_rho12.md and _alphaN0-2.md)"),
    ("", "Other probe fits (mechanism, mode and nugget checks; exploratory, not reported)"),
]


def source_version(sha: str) -> tuple[str, str]:
    for key, val in SOURCE_VERSIONS.items():
        if sha.startswith(key):
            return val
    return ("unknown", f"unknown source version {sha[:8]}")


# a0_intercept_prior.R overrides the prior in memory (its line 32) and names the
# variant in the label, so for those records the label is the source of truth.
A0_ALPHA_LABELS = {"a0_alpha_P2": "Normal(0, 2)", "a0_alpha_P10": "Normal(0, 10)", "a0_alpha_PC": "Normal(pooled log-odds, 2)"}


def intercept_prior(rec: dict) -> str:
    label = rec["extra"]["label"]
    if label in A0_ALPHA_LABELS:
        return A0_ALPHA_LABELS[label]
    model, kind = source_version(rec["stan_source_sha256"])
    if kind == "hard-coded":
        return "Normal(0, 2), hard-coded"
    if kind != "data":
        return kind
    path = Path(rec["data_path"])
    if not path.exists():
        return "data (file absent)"
    d = json.loads(path.read_text())
    mean = d.get("prior_alpha_mean", [0.0])
    centre = "0" if all(abs(float(m)) < 1e-12 for m in mean) else "pooled log-odds"
    return f"Normal({centre}, {d['prior_alpha_sd']:g})"


def row(rec: dict, path: Path) -> dict:
    g = rec["diagnostics"]
    td = g["treedepth"]
    gate = rec["gate"]
    reasons = gate["reasons"]
    if isinstance(reasons, str):
        reasons = [reasons]
    cfg = rec["sampler"]["rung_cfg"]
    cm = g.get("comembership_split")
    fam = next(p for p, _ in FAMILIES if rec["extra"]["label"].startswith(p))
    return {
        "family": fam.rstrip("_") or "other",
        "label": rec["extra"]["label"],
        "model": rec["model_file"],
        "intercept_prior": intercept_prior(rec),
        "written": rec["written"][:10],
        "seed": rec["sampler"]["seed"],
        "rung": rec["sampler"]["rung"],
        "iter": cfg["iter"], "warmup": cfg["warmup"], "adapt_delta": cfg["adapt_delta"], "max_treedepth": cfg["max_treedepth"],
        "rhat_max": g["rhat_max"], "rhat_worst": g["rhat_worst"],
        "ess_bulk_min": g["ess_bulk_min"], "ess_bulk_worst": g["ess_bulk_worst"],
        "ess_tail_min": g["ess_tail_min"], "ess_tail_worst": g["ess_tail_worst"],
        "n_div": g["n_div"], "n_total": g["n_total"], "pct_div": g["pct_div"],
        "treedepth_at_max": td["n_at_max"], "treedepth_pct": td["pct_at_max"], "treedepth_max_observed": td["max_observed"],
        "ebfmi_min": g["ebfmi_min"],
        "comembership_split": None if cm is None else cm["max_abs_diff"],
        "n_stuck": g["n_stuck"],
        "gate": "pass" if gate["pass"] else "fail",
        "reasons": ", ".join(reasons) if reasons else "",
        "earlier_rungs": "; ".join(f"rung {h['rung']}: {', '.join(h['gate']['reasons']) if isinstance(h['gate']['reasons'], list) else h['gate']['reasons']}"
                                   for h in rec["ladder_history"][:-1]) if len(rec["ladder_history"]) > 1 else "",
        "seconds": rec["sampler"]["seconds"],
        "software": rec["software"]["R"].replace("R version ", "R "),
        "record": path.name,
    }


def main() -> None:
    if not CACHE.exists():
        raise FileNotFoundError(f"{CACHE} missing: the fit cache is gitignored and lives on the fitting machine")
    recs = []
    for p in sorted(CACHE.glob("*.json")):
        rec = json.loads(p.read_text())
        if "diagnostics" not in rec or "extra" not in rec:
            continue
        recs.append(row(rec, p))
    if not recs:
        raise RuntimeError("no fit records found")
    df = pd.DataFrame(recs)
    order = {p.rstrip("_") or "other": i for i, (p, _) in enumerate(FAMILIES)}
    df["_o"] = df.family.map(order)
    df = df.sort_values(["_o", "label", "intercept_prior", "written"]).drop(columns="_o").reset_index(drop=True)
    df.to_csv(OUT_CSV, index=False)

    def fmt(r) -> str:
        cm = "" if pd.isna(r.comembership_split) else f"{r.comembership_split:.3f}"
        return (f"| {r.label} | {r.intercept_prior} | {r.written} | {r.rung} ({r.iter}/{r.warmup}, {r.adapt_delta}, {r.max_treedepth}) | "
                f"{r.rhat_max:.4f} ({r.rhat_worst}) | {r.ess_bulk_min:.0f} ({r.ess_bulk_worst}) | {r.ess_tail_min:.0f} ({r.ess_tail_worst}) | "
                f"{r.n_div}/{r.n_total} ({r.pct_div:.2f}%) | {r.treedepth_at_max}/{r.n_total} ({r.treedepth_pct:.1f}%; deepest {r.treedepth_max_observed}) | "
                f"{r.ebfmi_min:.3f} | {cm} | {r.gate}{(': ' + r.reasons) if r.reasons else ''} | {r.seconds:.0f} |")

    L = ["# Per-fit sampler diagnostics for the latent-group probe's Stan fits (rule 16)", "",
         "Produced by `analyses/probe_latent_groups/export_fit_diagnostics.py` from the fit records in `output/gp_cache/` "
         f"(gitignored; {len(df)} records with diagnostics). Every value is the one the gate read at the final rung; "
         "the gate thresholds are R-hat at most 1.01 on gated variables, bulk and tail ESS at least 400, no divergent transitions, "
         "E-BFMI at least 0.3 and, for A1, a per-chain co-membership split of at most 0.5. Treedepth saturation is reported, not gated. "
         "Rung settings are iterations per chain / warmup, adapt_delta and max_treedepth; four chains throughout. "
         "A label with two rows was fitted under two versions of the Stan source, with the class-intercept prior hard-coded as Normal(0, 2) "
         "and, after D25 (2026-10-18), read from the data (default Normal(0, 10)); the findings tables named per family report the "
         "data-driven version under the plain name and the hard-coded version as `*_alphaN0-2.md`. "
         "Descriptive export; exploratory design validation, not evidence about the record. Columns in `fit_diagnostics.csv`.", ""]
    for prefix, desc in FAMILIES:
        fam = prefix.rstrip("_") or "other"
        sub = df[df.family == fam]
        if sub.empty:
            continue
        n_pass = int((sub.gate == "pass").sum())
        L += [f"## {desc}", "", f"{len(sub)} records, {n_pass} pass the gate.", "",
              "| label | intercept prior | written | rung (iter/warmup, adapt_delta, max_treedepth) | R-hat max (worst) | bulk ESS min (worst) | tail ESS min (worst) | divergences | treedepth at max | E-BFMI min | co-membership split | gate | s |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        L += [fmt(r) for r in sub.itertuples()]
        L.append("")
    OUT_MD.write_text("\n".join(L))
    print(f"wrote {OUT_MD} and {OUT_CSV}: {len(df)} records")
    print(df.groupby(["family", "gate"]).size())


if __name__ == "__main__":
    main()
