# phasecheck: do your phases behave as bounded groups?

`phasecheck` applies the phase tests of this paper to any table of assemblages.
It needs class counts, a latitude and longitude, and a phase for each
assemblage. It answers a series of questions in order and writes a report.

A plain run answers five questions that compare the phase lines with other
lines on the same map. With `--model` it also calibrates a model of copying
across distance with no groups to your record and answers four more: how much
difference distance alone produces on your map, and whether the comparisons
could detect a boundary on a record like yours. Without the model, a
difference between phases is not, by itself, evidence of a boundary.

## Input

One row per assemblage, in `.xlsx`, `.xls`, `.csv` or `.tsv`:

| name | latitude | longitude | phase | class A | class B | ... |
|---|---|---|---|---|---|---|
| Site 1 | 35.379 | -90.390 | Parkin | 404 | 208 | ... |

- `latitude` and `longitude` are decimal degrees.
- Every other column is read as a count of sherds (whole numbers; enter 0 where
  a class is absent). Percentages are refused when they can be recognized
  (fractions, or every row summing to about 100), because the method carries
  the sampling uncertainty of the counts. Enter counts.
- Columns that are not counts must be left out. A year, an identifier, a total
  or an elevation read as a class gives a plausible, wrong report. The reader
  refuses the ones it can recognize and the report lists the classes it used;
  check that list.
- Every phase needs at least two assemblages, and the table at least six
  assemblages, two phases and two classes with sherds.
- The study area must be under 1,000 km across.

If your columns have other names, say so: `--name-col Site --lat-col Lat
--lon-col Long --phase-col Phase`. When the table holds other columns as well,
name the count columns with `--class-cols "A,B,C"` or leave the others out with
`--ignore-cols "notes,elevation"`.

## Run

```bash
pip install -e .
phasecheck my_assemblages.xlsx --out my_report
```

`python -m phasecheck ...` does the same without the installed command. The
package itself needs only numpy, scipy, pandas and openpyxl; installing this
repository also installs what the paper's own analyses need.

Options:

- `--min-count 75` drops assemblages with fewer sherds and lists them. The
  paper justified 75 for its own record; yours needs its own check.
- `--distance river_km.csv` supplies a square matrix of distances in km
  (assemblage names as header and first column), for example along rivers or
  trails. It is used for the boundary excess, the one measure that matches
  pairs by distance. The alternative divisions of the map are always built
  from the coordinates. Without it, straight-line distance is used throughout.
- `--draws`, `--alternatives`, `--seed` set the posterior draws of the counts
  (default 2000), the number of alternative divisions of each kind (default
  300), and the random seed (default 0). The seed moves questions 2 and 3
  only. Question 4 uses a fixed seed of its own, and the copying model uses
  fixed seed families of its own, so `--seed` does not change them.
- `--sheet` names the worksheet of a spreadsheet with several.
- `--totals-of-100-are-counts` accepts a table whose rows all sum to about
  100 as counts, when they really are counts.

The paper's own record is the worked example:

```bash
phasecheck examples/phasecheck/st_francis_basin.csv \
  --distance examples/phasecheck/st_francis_river_km.csv \
  --order-col seriation_order --seed 86000 --out basin_report
```

That run reproduces the paper's phase-boundary comparison to the printed digit
(`tests/phasecheck/test_phasecheck.py` checks it).

## The copying model (`--model`)

```bash
phasecheck my_assemblages.xlsx --model --jobs 8 --out my_report
```

- **Calibration.** The model's three settings (learners per assemblage,
  innovation rate, mixing rate) are tuned to your record's diversity: 297
  combinations are screened, the matches are confirmed on fresh runs, and the
  one that gives the phases the largest difference is used, so the model with
  no groups is compared at its most favorable for that comparison. If nothing
  reproduces your record's diversity, the tool stops and says so. On the
  paper's 28 assemblages calibration takes under a minute on 12 cores and
  about five minutes on one; the whole `--model` run takes about four minutes
  on 12 cores.
- **Settings that fit equally well can disagree.** Diversity often constrains
  the model weakly. The report therefore shows the between-phase difference at
  three confirmed settings (the largest, the middle and the smallest) when at
  least three were confirmed, and says so plainly when the answer depends on
  which is used, or when the model makes assemblages more different than they
  are and so does not fit. With `--rates` only the supplied setting is shown.
- **`--length-km`** is the distance over which copying falls off. The default is
  a fifth of the largest distance between assemblages (along the supplied
  distances, when `--distance` is given). It is an assumption,
  not an estimate; rerun at half and double to see what depends on it.
- **`--order-col`** names a column that places assemblages along a sequence (a
  date or a seriation score; larger is later). Only the order is used, not the
  spacing. The column is never read as a class, so name it even in a run
  without the model when the table holds one (the worked example does).
  Without it the model treats the assemblages as contemporaneous.
- **`--rates "2000,0.001,0.02"`** skips calibration and runs the model at
  settings you supply (learners, innovation, mixing). It implies `--model`.
- **`--model-runs`** (default 300, 50 to 5,000) sets the runs of the model in
  questions 7 and 8 and the simulated records per setting in question 9.
- **`--jobs`** sets the processor cores used for the calibration.

The paper's settings on the worked example:

```bash
phasecheck examples/phasecheck/st_francis_basin.csv \
  --distance examples/phasecheck/st_francis_river_km.csv --order-col seriation_order \
  --rates 2000,0.001,0.02 --length-km 24 --out basin_model_report
```

That run gives values close to the paper's. The paper's own values are
reproduced exactly from Python with the paper's seeds, where 7 of 300 runs
reach the observed between-phase F_ST and none reaches the observed difference
in sherds (`cp.compare_with_model(data, rates, seed=84000)`;
`tests/phasecheck/test_copying_model.py` checks it, and checks that
calibration finds the paper's 51 matched and 48 confirmed combinations). The
command line uses its own seeds, kept apart from the calibration's. The paper
chose among the confirmed combinations by the difference between two spatial
clusters; the tool chooses by the difference between the phases, so a fresh
calibration of this record at `--length-km 24` selects another of the paper's
confirmed combinations (10,000 learners, 0.0005, 0.002).

## What the report answers

1. **Does location alone recover the phases?** Clustering on coordinates and
   composition at varying weights, scored against the phases.
2. **Do the phase boundaries separate the assemblages better than other
   lines?** Cultural F_ST and the boundary excess for the phases, beside
   divisions of the same sizes around random centers and made compact.
3. **Does each phase stand apart from the rest?** The same comparison for each
   phase against the others together.
4. **Does each assemblage fit its own phase?** Each profile against the pooled
   profile of every phase.
5. **How different are the phases, in sherds?** The share of sherds that would
   have to change class to make two phases alike.

With `--model`:

6. **Can copying with no groups reproduce this record's diversity?** The
   calibration.
7. **Does the record differ between phases by more than copying across
   distance produces?** The record beside 300 runs of the model, at the
   calibrated settings or the ones you supply.
8. **Would assemblages made with no groups sort into these phases?** Simulated
   records clustered as the real one was.
9. **Could this comparison detect a restriction at the phase boundaries?**
   Records simulated with and without a restriction (300 per setting), scored
   as the real record is against 200 alternative divisions of each kind. These
   figures change by several hundredths with the random seed.

It writes `report.md`, `recovery.csv`, `assemblage_fit.csv` and `results.json`,
and with the model `calibration.csv` and `model_comparison.csv`.

What the tool does not include from the paper: local innovation as an
alternative to a restriction, the map of compositional change, the
settlement-gap test, sensitivity to the copying length, a full account of
every confirmed combination (three are shown), and figures.

## How to read the probabilities

Each probability is the phase lines' standing among alternative lines: the
share of comparisons, over posterior draws of the counts, in which the phase
lines separate the assemblages better. It is not a significance level. An
arbitrary division of a map with any spatial pattern can land anywhere between
0 and 1, depending on how its lines happen to run. A value near one half says
the phase lines are unremarkable. A value near 1 is necessary for the lines to
stand apart but is not sufficient: arbitrary lines reach it too, most easily
against compact divisions, which are few (the report counts the distinct ones).
The paper found that the comparison discriminates only moderately on a record
of 28 assemblages. The second decimal place changes from seed to seed.

The tests ask whether learning was bounded where the phase lines fall. They do
not test whether a society or polity existed.

## From Python

```python
import phasecheck as pc
data = pc.read_table("my_assemblages.xlsx", min_count=75)
result = pc.boundary_comparison(data, draws=2000, n_alt=300, seed=1)
print(result["p_fst_around"], result["p_fst_compact"])
```
