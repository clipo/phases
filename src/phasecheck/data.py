"""Read and check the input table: one row per assemblage.

Every question starts from the `Dataset` built here. `read_table` refuses
input the tests cannot use rather than returning a plausible, wrong report:
columns that do not look like sherd counts (identifiers, years, totals,
flags, text, fractions), tables that look like percentages, duplicate names,
too few assemblages, phases or classes, and study areas too large for the
flat-map approximation of `geo.planar_km`. Every refusal is a `ValueError`
whose message says what to change; the command line prints it and exits 2.

Units: latitude and longitude in decimal degrees; distances in km. Counts are
whole sherds. The optional sequence order is converted to evenly spaced
ranks, 0 earliest to 1 latest, and is used only by the copying model.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import geo


@dataclass
class Dataset:
    """Assemblages by classes, with a location and a phase for each assemblage.

    Built by `read_table`, which guarantees the invariants below; the
    questions assume them.

    Attributes
    ----------
    names : list of str
        Unique assemblage names, in table order.
    lat, lon : numpy.ndarray, shape (n,)
        Decimal degrees.
    phases : numpy.ndarray of str, shape (n,)
        Phase of each assemblage; at least two phases, each with at least
        two assemblages.
    classes : list
        Class (column) names, only those with at least one sherd.
    counts : numpy.ndarray of int, shape (n, classes)
        Sherd counts; every row has at least one sherd.
    dist : numpy.ndarray, shape (n, n)
        Pairwise distance in km: great-circle, or the supplied matrix.
    dist_kind : str
        "straight-line", or "supplied (<file name>)".
    order : numpy.ndarray or None
        Position along a sequence, 0 earliest to 1 latest, evenly spaced by
        rank with ties sharing a position; None means contemporaneous.
    dropped : list of (str, int)
        Assemblages removed for falling below the minimum count, with their
        sherd totals.
    notes : list of str
        Things a reader of the report must be told.
    """

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
        """Distinct phase names in sorted order; this order numbers the phases everywhere."""
        return sorted(set(self.phases.tolist()))

    @property
    def phase_index(self) -> np.ndarray:
        """Phase of each assemblage as an integer, its position in `phase_names`."""
        order = self.phase_names
        return np.array([order.index(p) for p in self.phases.tolist()])

    @property
    def pts(self) -> np.ndarray:
        """Flat east/north coordinates in km (`geo.planar_km`), recomputed on each access."""
        return geo.planar_km(self.lat, self.lon)


def _read_any(path: Path, sheet) -> pd.DataFrame:
    """Read a spreadsheet or delimited text file by its extension; raise ValueError on any failure."""
    if not path.is_file():
        raise ValueError(f"no such file: {path}")
    suffix = path.suffix.lower()
    # Only an empty cell is missing. "NA", "None" and the like are legitimate
    # phase and site names, which pandas would otherwise read as blanks.
    na = {"keep_default_na": False, "na_values": [""]}
    try:
        if suffix in (".xlsx", ".xlsm", ".xls"):
            return pd.read_excel(path, sheet_name=0 if sheet is None else sheet, **na)
        if suffix in (".csv", ".txt", ".tsv", ".tab"):
            sep = "\t" if suffix in (".tsv", ".tab") else ","
            try:
                return pd.read_csv(path, sep=sep, **na)
            except UnicodeDecodeError:
                return pd.read_csv(path, sep=sep, encoding="cp1252", **na)   # the usual Excel export
    except ValueError:                           # already a message for the reader; pass it on
        raise
    except Exception as e:                       # a damaged or mislabeled file
        raise ValueError(f"could not read {path.name}: {e}") from e
    raise ValueError(f"unrecognized file type {suffix!r}; use .xlsx, .xls, .csv or .tsv")


# Column names (or first words of names) that mark a column as something other
# than a class count. Checked only when the class columns were not named.
NOT_A_CLASS = ("id", "no", "num", "number", "year", "date", "elev", "elevation", "area", "size", "notes",
               "note", "comment", "comments", "utm", "easting", "northing", "x", "y", "total", "sum", "n",
               "county", "state", "quad", "trinomial")


def _check_class_columns(raw: pd.DataFrame, chosen_by_user: bool) -> None:
    """Refuse columns that are not sherd counts. A wrong class column gives a plausible, wrong report.

    Every column must be numeric, complete, non-negative and whole. When the
    user did not name the class columns (`chosen_by_user` False), columns are
    also refused that are named like a non-count (`NOT_A_CLASS`), run through
    consecutive numbers, are nearly constant around a large value, hold only
    0 and 1, or equal the sum of the others. Raises ValueError naming the
    column; returns None when every column passes.
    """
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
        # The name, or its first word ("year built", "elev_m" -> "year", "elev"), marks a non-count.
        if str(c).strip().lower().replace("_", " ").split(" ")[0] in NOT_A_CLASS or str(c).strip().lower() in NOT_A_CLASS:
            raise ValueError(f"column {c!r} is named like something other than a class count." + hint)
        # A permutation of consecutive values (1, 2, 3, ...) is an identifier. More than 3 rows,
        # so that a short table of genuine counts is not refused by chance.
        if len(x) > 3 and len(set(x)) == len(x) and np.all(np.diff(np.sort(x)) == 1):
            raise ValueError(f"column {c!r} is a run of consecutive numbers and looks like an identifier." + hint)
        # Range under a fifth of the minimum, every value at least 100: a year, an elevation or
        # a measurement. Sherd counts vary far more than that from assemblage to assemblage.
        if x.min() > 0 and (x.max() - x.min()) < 0.2 * x.min() and x.min() >= 100:
            raise ValueError(f"column {c!r} is nearly constant ({x.min():g} to {x.max():g}) and looks like a "
                             "year, elevation or other measurement." + hint)
        if set(np.unique(x)) <= {0.0, 1.0} and len(x) > 3:
            raise ValueError(f"column {c!r} holds only 0 and 1 and looks like a presence flag." + hint)
    # A column equal to the sum of all the others is a total. Needs at least three columns,
    # since with two, each one "totals" the other whenever they are equal.
    if not chosen_by_user and raw.shape[1] > 2:
        vals = raw.apply(pd.to_numeric, errors="coerce").to_numpy(float)
        for j, c in enumerate(raw.columns):
            rest = np.delete(vals, j, axis=1).sum(1)
            if vals[:, j].sum() > 0 and np.allclose(vals[:, j], rest):
                raise ValueError(f"column {c!r} equals the sum of the other count columns and looks like a "
                                 "total." + hint)


def read_distance(path, names: list) -> np.ndarray:
    """A square distance matrix in km with assemblage names as header and first column.

    Parameters
    ----------
    path : str or path-like
        A table readable by `_read_any`; extra rows and columns are ignored.
    names : list of str
        Assemblage names; the result follows this order.

    Returns
    -------
    numpy.ndarray, shape (len(names), len(names))

    Raises
    ------
    ValueError
        If a name has no row or column, or the matrix is not finite,
        non-negative, symmetric (to 1e-6) and zero on the diagonal.
    """
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

    Parameters
    ----------
    path : str or path-like
        .xlsx, .xlsm, .xls, .csv, .txt, .tsv or .tab. Only an empty cell
        counts as missing, so "NA" can be a phase name.
    name_col, lat_col, lon_col, phase_col : str
        Names of the required columns.
    class_cols : list of str, optional
        The count columns. When given, the name-based and shape-based
        checks for non-count columns are skipped.
    ignore_cols : list of str, optional
        Columns to leave out when the count columns are not named.
    sheet : str or int, optional
        Worksheet of a spreadsheet; default the first.
    min_count : int
        Assemblages with fewer sherds are dropped and listed (default 0;
        an assemblage with no sherds is always dropped). The paper used 75
        for its own record; another record needs its own check.
    distance : str or path-like, optional
        Square distance matrix in km (see `read_distance`).
    totals_of_100_are_counts : bool
        Accept a table whose rows all sum to about 100 as counts.
    order_col : str, optional
        Column placing assemblages along a sequence (larger is later). It is
        not read as a class. Only the order is kept.

    Returns
    -------
    Dataset

    Raises
    ------
    ValueError
        On any input the tests cannot use: a missing or duplicated column, a
        column that is not counts, empty required cells, duplicate names,
        coordinates out of range, rows that look like percentages, fewer
        than 6 assemblages or 2 classes after dropping, one phase or a phase
        with one assemblage, fewer distinct locations than phases, a study
        area beyond `geo.MAX_EXTENT_KM`, or a bad distance matrix.
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
                             # +2: spreadsheet row numbers count the header and start at 1;
                             # the index keeps the original positions after dropna
                             f"{[int(i) + 2 for i in df.index[df[col].isna().to_numpy()]][:8]}")
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
    # Whole-number percentages need not sum to exactly 100: each class can round by up to 0.5,
    # so allow half a sherd per class. More than 3 rows, so a tiny table is not refused by chance.
    slack = max(1, counts.shape[1] // 2)
    if (not totals_of_100_are_counts and len(row_totals) > 3
            and np.all(np.abs(row_totals - 100) <= slack)):
        raise ValueError("every row sums to about 100, so the table looks like percentages. The method needs "
                         "the sherd counts, because it carries their sampling uncertainty. If these really "
                         "are counts of about 100 sherds each, say so with totals_of_100_are_counts "
                         "(--totals-of-100-are-counts).")
    phases = df[phase_col].astype(str).str.strip().to_numpy()
    notes = []
    lowered = {}                                  # phase names that differ only in case: kept apart, noted
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
    keep = totals >= max(int(min_count), 1)       # an empty assemblage is always dropped
    dropped = [(n, int(t)) for n, t, k in zip(names, totals, keep) if not k]
    names = [n for n, k in zip(names, keep) if k]
    lat, lon, phases, counts = lat[keep], lon[keep], phases[keep], counts[keep]
    if order is not None:
        from scipy.stats import rankdata
        o = order[keep]
        # Evenly spaced by rank: only the order is used, not the spacing of the
        # values. Tied values share a position, so row order cannot matter.
        order = (rankdata(o, method="average") - 1.0) / (len(o) - 1.0)
    live = counts.sum(0) > 0                      # drop classes left with no sherds
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
    # Above 95 percent in one class the Jeffreys prior of the posterior draws, not the data,
    # supplies much of the variation between assemblages; warn, do not refuse.
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
