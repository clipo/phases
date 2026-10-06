"""Read and check the input table: one row per assemblage."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import geo


@dataclass
class Dataset:
    """Assemblages by classes, with a location and a phase for each assemblage."""
    names: list
    lat: np.ndarray
    lon: np.ndarray
    phases: np.ndarray            # phase name of each assemblage
    classes: list
    counts: np.ndarray            # assemblages x classes, non-negative integers
    dist: np.ndarray              # pairwise distance in km
    dist_kind: str = "straight-line"
    order: np.ndarray | None = None   # position along a sequence, 0 earliest to 1 latest; None = contemporaneous
    dropped: list = field(default_factory=list)   # (name, total) removed by the minimum count
    notes: list = field(default_factory=list)     # things a reader of the report must be told

    @property
    def phase_names(self) -> list:
        return sorted(set(self.phases.tolist()))

    @property
    def phase_index(self) -> np.ndarray:
        order = self.phase_names
        return np.array([order.index(p) for p in self.phases.tolist()])

    @property
    def pts(self) -> np.ndarray:
        return geo.planar_km(self.lat, self.lon)


def _read_any(path: Path, sheet) -> pd.DataFrame:
    if not path.is_file():
        raise ValueError(f"no such file: {path}")
    suffix = path.suffix.lower()
    # Only an empty cell is missing. "NA", "None" and the like are legitimate
    # phase and site names, which pandas would otherwise read as blanks.
    na = dict(keep_default_na=False, na_values=[""])
    try:
        if suffix in (".xlsx", ".xlsm", ".xls"):
            return pd.read_excel(path, sheet_name=0 if sheet is None else sheet, **na)
        if suffix in (".csv", ".txt", ".tsv", ".tab"):
            sep = "\t" if suffix in (".tsv", ".tab") else ","
            try:
                return pd.read_csv(path, sep=sep, **na)
            except UnicodeDecodeError:
                return pd.read_csv(path, sep=sep, encoding="cp1252", **na)   # the usual Excel export
    except ValueError:
        raise
    except Exception as e:                       # a damaged or mislabeled file
        raise ValueError(f"could not read {path.name}: {e}") from e
    raise ValueError(f"unrecognized file type {suffix!r}; use .xlsx, .xls, .csv or .tsv")


NOT_A_CLASS = ("id", "no", "num", "number", "year", "date", "elev", "elevation", "area", "size", "notes",
               "note", "comment", "comments", "utm", "easting", "northing", "x", "y", "total", "sum", "n",
               "county", "state", "quad", "trinomial")


def _check_class_columns(raw: pd.DataFrame, chosen_by_user: bool) -> None:
    """Refuse columns that are not sherd counts. A wrong class column gives a plausible, wrong report."""
    hint = (" Name the count columns with class_cols (--class-cols), or leave columns out with "
            "ignore_cols (--ignore-cols).")
    for c in raw.columns:
        v = pd.to_numeric(raw[c], errors="coerce")
        if v.isna().all():
            raise ValueError(f"column {c!r} holds text, not sherd counts." + hint)
        if v.isna().any():
            raise ValueError(f"column {c!r} has blank or non-numeric cells; enter 0 where a class is absent."
                             + ("" if chosen_by_user else hint))
        if np.any(v < 0):
            raise ValueError(f"column {c!r} has negative values; counts must be non-negative")
        if not np.allclose(v, np.round(v)):
            raise ValueError(f"column {c!r} has fractional values; counts must be whole numbers of sherds, "
                             "not percentages, proportions or measurements." + ("" if chosen_by_user else hint))
        if chosen_by_user:
            continue
        x = v.to_numpy(float)
        if str(c).strip().lower().replace("_", " ").split(" ")[0] in NOT_A_CLASS or str(c).strip().lower() in NOT_A_CLASS:
            raise ValueError(f"column {c!r} is named like something other than a class count." + hint)
        if len(x) > 3 and len(set(x)) == len(x) and np.all(np.diff(np.sort(x)) == 1):
            raise ValueError(f"column {c!r} is a run of consecutive numbers and looks like an identifier." + hint)
        if x.min() > 0 and (x.max() - x.min()) < 0.2 * x.min() and x.min() >= 100:
            raise ValueError(f"column {c!r} is nearly constant ({x.min():g} to {x.max():g}) and looks like a "
                             "year, elevation or other measurement." + hint)
        if set(np.unique(x)) <= {0.0, 1.0} and len(x) > 3:
            raise ValueError(f"column {c!r} holds only 0 and 1 and looks like a presence flag." + hint)
    if not chosen_by_user and raw.shape[1] > 2:
        vals = raw.apply(pd.to_numeric, errors="coerce").to_numpy(float)
        for j, c in enumerate(raw.columns):
            rest = np.delete(vals, j, axis=1).sum(1)
            if vals[:, j].sum() > 0 and np.allclose(vals[:, j], rest):
                raise ValueError(f"column {c!r} equals the sum of the other count columns and looks like a "
                                 "total." + hint)


def read_distance(path, names: list) -> np.ndarray:
    """A square distance matrix in km with assemblage names as header and first column."""
    d = _read_any(Path(path), None)
    d = d.set_index(d.columns[0])
    d.index = d.index.astype(str)
    d.columns = d.columns.astype(str)
    missing = [n for n in names if n not in d.index or n not in d.columns]
    if missing:
        raise ValueError(f"the distance matrix has no row or column for: {', '.join(missing[:8])}")
    m = d.loc[names, names].to_numpy(float)
    if not np.all(np.isfinite(m)) or np.any(m < 0):
        raise ValueError("the distance matrix must hold finite, non-negative numbers")
    if not np.allclose(m, m.T, rtol=1e-6, atol=1e-6):
        raise ValueError("the distance matrix is not symmetric")
    if not np.allclose(np.diag(m), 0.0, atol=1e-9):
        raise ValueError("the distance matrix must have zeros on its diagonal")
    return m


def read_table(path, *, name_col: str = "name", lat_col: str = "latitude", lon_col: str = "longitude",
               phase_col: str = "phase", class_cols=None, ignore_cols=None, sheet=None, min_count: int = 0,
               distance=None, totals_of_100_are_counts: bool = False, order_col=None) -> Dataset:
    """Read a table of assemblages and refuse anything the tests cannot use.

    One row per assemblage. Required columns: a name, latitude and longitude in
    decimal degrees, and a phase. Every other column is read as a class count
    unless `class_cols` names the count columns or `ignore_cols` leaves some
    out; a column that does not look like counts is refused, because a year,
    an identifier or a total read as a class gives a plausible, wrong report.
    Assemblages with fewer than
    `min_count` sherds in those classes are dropped and listed. `distance` is
    an optional file holding a square matrix of distances in km (for example
    along rivers or trails). It is used for the boundary excess, the one
    measure that matches pairs by distance; the alternative divisions of the
    map are always built from the coordinates. Without it, straight-line
    distance is used throughout.
    """
    df = _read_any(Path(path), sheet)
    df.columns = [str(c).strip() for c in df.columns]
    need = {"name": name_col, "latitude": lat_col, "longitude": lon_col, "phase": phase_col}
    absent = [f"{k} (looked for a column named {v!r})" for k, v in need.items() if v not in df.columns]
    if absent:
        raise ValueError("missing column(s): " + "; ".join(absent) + ". Columns found: " + ", ".join(df.columns))
    chosen = class_cols is not None
    ignore_cols = list(ignore_cols or [])
    if order_col is not None:
        if order_col not in df.columns:
            raise ValueError(f"order column {order_col!r} is not in the table")
        ignore_cols.append(order_col)
    gone = [c for c in ignore_cols if c not in df.columns]
    if gone:
        raise ValueError("ignored column(s) not in the table: " + ", ".join(gone))
    if not chosen:
        class_cols = [c for c in df.columns if c not in need.values() and c not in ignore_cols
                      and not str(c).startswith("Unnamed:")]
    else:
        class_cols = [str(c).strip() for c in class_cols]
        gone = [c for c in class_cols if c not in df.columns]
        if gone:
            raise ValueError("class column(s) not in the table: " + ", ".join(gone))
        twice = sorted({c for c in class_cols if class_cols.count(c) > 1})
        if twice:
            raise ValueError("class column(s) named more than once: " + ", ".join(twice))
        clash = [c for c in class_cols if c in need.values()]
        if clash:
            raise ValueError("these are the name, coordinate or phase columns, not classes: " + ", ".join(clash))
    if len(set(df.columns)) != len(df.columns):
        raise ValueError("two columns share a name; column names must be unique")
    if len(class_cols) < 2:
        raise ValueError("at least two class columns are needed")

    df = df.dropna(how="all")
    for col in need.values():
        if df[col].isna().any():
            raise ValueError(f"column {col!r} has empty cells in row(s) "
                             f"{[int(i) + 2 for i in np.flatnonzero(df[col].isna().to_numpy())][:8]}")
    names = df[name_col].astype(str).str.strip().tolist()
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise ValueError("assemblage names must be unique; repeated: " + ", ".join(dup[:8]))
    try:
        lat, lon = df[lat_col].astype(float).to_numpy(), df[lon_col].astype(float).to_numpy()
    except (TypeError, ValueError) as e:
        raise ValueError("latitude and longitude must be decimal degrees") from e
    if np.any(np.abs(lat) > 90) or np.any(np.abs(lon) > 180):
        raise ValueError("latitude must lie in -90..90 and longitude in -180..180 (decimal degrees); "
                         "check that the two columns are not swapped or in projected units")
    _check_class_columns(df[class_cols], chosen)
    counts = np.round(df[class_cols].apply(pd.to_numeric).to_numpy(float)).astype(int)
    row_totals = counts.sum(1)
    slack = max(1, counts.shape[1] // 2)          # whole-number percentages need not sum to exactly 100
    if (not totals_of_100_are_counts and len(row_totals) > 3
            and np.all(np.abs(row_totals - 100) <= slack)):
        raise ValueError("every row sums to about 100, so the table looks like percentages. The method needs "
                         "the sherd counts, because it carries their sampling uncertainty. If these really "
                         "are counts of about 100 sherds each, say so with totals_of_100_are_counts "
                         "(--totals-of-100-are-counts).")
    phases = df[phase_col].astype(str).str.strip().to_numpy()
    notes = []
    lowered = {}
    for p in sorted(set(phases.tolist())):
        lowered.setdefault(p.lower(), []).append(p)
    same = [v for v in lowered.values() if len(v) > 1]
    if same:
        notes.append("Phase names that differ only in capitalization were kept as separate phases: "
                     + "; ".join(" / ".join(v) for v in same) + ".")

    order = None
    if order_col is not None:
        o = pd.to_numeric(df[order_col], errors="coerce").to_numpy(float)
        if not np.isfinite(o).all():
            raise ValueError(f"order column {order_col!r} must hold a number for every assemblage "
                             "(a date, a seriation score or a rank; larger is later)")
        order = o
    totals = counts.sum(1)
    keep = totals >= max(int(min_count), 1)
    dropped = [(n, int(t)) for n, t, k in zip(names, totals, keep) if not k]
    names = [n for n, k in zip(names, keep) if k]
    lat, lon, phases, counts = lat[keep], lon[keep], phases[keep], counts[keep]
    if order is not None:
        from scipy.stats import rankdata
        o = order[keep]
        # Evenly spaced by rank: only the order is used, not the spacing of the
        # values. Tied values share a position, so row order cannot matter.
        order = (rankdata(o, method="average") - 1.0) / (len(o) - 1.0)
    live = counts.sum(0) > 0
    classes = [c for c, k in zip(class_cols, live) if k]
    counts = counts[:, live]

    if len(names) < 6:
        raise ValueError(f"only {len(names)} assemblages remain; the comparisons need at least 6")
    if len(classes) < 2:
        raise ValueError("fewer than two classes have any sherds")
    sizes = pd.Series(phases).value_counts()
    if len(sizes) < 2:
        raise ValueError("the table names a single phase; at least two are needed")
    if (sizes < 2).any():
        raise ValueError("every phase needs at least two assemblages; too few in: "
                         + ", ".join(f"{p} ({n})" for p, n in sizes[sizes < 2].items()))
    if len(np.unique(np.column_stack([lat, lon]), axis=0)) < len(sizes):
        raise ValueError("there are fewer distinct locations than phases")

    geo.planar_km(lat, lon)                       # raises if the region is too large
    if distance is None:
        dist, kind = geo.great_circle_km(lat, lon), "straight-line"
    else:
        dist, kind = read_distance(distance, names), f"supplied ({Path(distance).name})"
    top = counts.sum(0).max() / counts.sum()
    if top > 0.95:
        notes.append(f"One class holds {100 * top:.0f} percent of all sherds. With so little variation to "
                     "divide, the prior on the proportions supplies much of the difference between "
                     "assemblages, and every result below should be read as weak.")
    empty = [c for c, k in zip(class_cols, live) if not k]
    if empty:
        notes.append("Columns with no sherds were left out: " + ", ".join(map(str, empty)) + ".")
    return Dataset(names=names, lat=lat, lon=lon, phases=phases, classes=classes, counts=counts,
                   dist=dist, dist_kind=kind, order=order, dropped=dropped, notes=notes)
