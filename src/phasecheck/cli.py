"""Command line: `phasecheck TABLE --out DIR [options]`.

Parses the options, reads the table with `data.read_table`, and hands it to
`report.run`, which answers questions 1 to 5 and, with `--model` or
`--rates`, questions 6 to 9. Bad input is reported as one line on standard
error with exit status 2, never as a traceback; success prints where the
report was written and exits 0.

Seeds: `--seed` (default 0) drives the posterior draws and alternative
divisions of questions 2 and 3 (see `questions.compare_division` for the
offsets). Question 4 uses its own fixed seed (96) and the copying model its
own fixed seed families (see `report.run`), so `--seed` does not move them.
"""
from __future__ import annotations

import argparse
import sys

from .data import read_table
from .report import run


def main(argv=None) -> int:
    """Run phasecheck from the command line.

    Parameters
    ----------
    argv : list of str, optional
        The arguments without the program name; `None` reads `sys.argv`.

    Returns
    -------
    int
        0 when the report was written; 2 when the input or options were
        refused (the reason is printed to standard error). argparse itself
        exits with status 2 on an unknown option.
    """
    ap = argparse.ArgumentParser(
        prog="phasecheck",
        description="Ask whether archaeological phases behave as bounded groups. TABLE has one row per "
                    "assemblage: a name, latitude and longitude in decimal degrees, a phase, and one "
                    "column of sherd counts per class.")
    ap.add_argument("table", help="spreadsheet or text table (.xlsx, .xls, .csv, .tsv)")
    ap.add_argument("--out", default="phasecheck_report", help="directory for the report (default: %(default)s)")
    ap.add_argument("--sheet", default=None, help="worksheet name, for a spreadsheet with several")
    ap.add_argument("--name-col", default="name", help="column holding assemblage names")
    ap.add_argument("--lat-col", default="latitude", help="column holding latitude, decimal degrees")
    ap.add_argument("--lon-col", default="longitude", help="column holding longitude, decimal degrees")
    ap.add_argument("--phase-col", default="phase", help="column holding phase labels")
    ap.add_argument("--class-cols", default=None,
                    help="comma-separated count columns (default: every column not named above)")
    ap.add_argument("--ignore-cols", default=None,
                    help="comma-separated columns to leave out (notes, identifiers, totals, measurements)")
    ap.add_argument("--min-count", type=int, default=0,
                    help="drop assemblages with fewer sherds than this (default: keep all)")
    ap.add_argument("--distance", default=None,
                    help="file with a square matrix of distances in km between assemblages, names as header "
                         "and first column (default: straight-line distance)")
    ap.add_argument("--totals-of-100-are-counts", action="store_true",
                    help="accept a table in which every row sums to about 100 as counts, not percentages")
    ap.add_argument("--model", action="store_true",
                    help="also calibrate the copying model with no groups to this record and answer "
                         "questions 6 to 9 (takes minutes to tens of minutes)")
    ap.add_argument("--length-km", type=float, default=None,
                    help="distance over which copying falls off (default: a fifth of the largest distance)")
    ap.add_argument("--order-col", default=None,
                    help="column placing assemblages along a sequence (a date or seriation score; larger is "
                         "later). Only the order is used, not the spacing. Without it the model treats "
                         "them as contemporaneous")
    ap.add_argument("--rates", default=None,
                    help="skip calibration and run the model at 'learners,innovation,mixing', for example "
                         "'2000,0.001,0.02'")
    ap.add_argument("--jobs", type=int, default=1, help="processor cores for the calibration")
    ap.add_argument("--model-runs", type=int, default=300,
                    help="runs of the calibrated model, and simulated records per setting in question 9")
    ap.add_argument("--draws", type=int, default=2000, help="posterior draws of the counts")
    ap.add_argument("--alternatives", type=int, default=300, help="alternative divisions of each kind")
    ap.add_argument("--seed", type=int, default=0,
                    help="random seed for the comparisons of questions 2 and 3 (default: %(default)s)")
    a = ap.parse_args(argv)
    try:
        data = read_table(a.table, name_col=a.name_col, lat_col=a.lat_col, lon_col=a.lon_col,
                          phase_col=a.phase_col, sheet=a.sheet, min_count=a.min_count, distance=a.distance,
                          totals_of_100_are_counts=a.totals_of_100_are_counts, order_col=a.order_col,
                          class_cols=[c.strip() for c in a.class_cols.split(",")] if a.class_cols else None,
                          ignore_cols=[c.strip() for c in a.ignore_cols.split(",")] if a.ignore_cols else None)
        rates = None
        if a.rates:
            # "2000,0.001,0.02" -> learners, innovation, mixing. The learner count
            # goes through float so that "2e3" is accepted.
            parts = [p.strip() for p in a.rates.split(",")]
            if len(parts) != 3:
                raise ValueError("--rates takes three numbers: learners,innovation,mixing")
            rates = {"n_ind": int(float(parts[0])), "innovation": float(parts[1]), "mixing": float(parts[2])}
        # The upper bound keeps run seeds inside their seed family (copying.MAX_RUNS).
        if not 50 <= a.model_runs <= 5000:
            raise ValueError("--model-runs must lie between 50 and 5000")
        res = run(data, a.out, draws=a.draws, n_alt=a.alternatives, seed=a.seed, source=a.table,
                  model=a.model or rates is not None, rates=rates, length=a.length_km, jobs=a.jobs,
                  # One option sets both: runs of question 7 and records per setting of question 9.
                  model_runs=a.model_runs, power_records=a.model_runs,
                  progress=lambda msg: print(msg, flush=True))
    except (ValueError, OSError) as e:
        print(f"phasecheck: {e}", file=sys.stderr)
        return 2
    print(f"report written to {res['report']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
