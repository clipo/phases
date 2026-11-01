# The canonical basin membership cannot be regenerated from the repository

Date: 2026-08-31. Found while attempting item 0 task 1 step 4 (regenerate the
baseline). Not a modeling defect; a provenance defect.

**Addendum 2026-09-22.** The counts below (29 curated, 92 broad) are the
sets as they stood on 2026-09-04. Since the 2026-09-21 rulings (Parchman out on
geography, 75-decorated-sherd minimum, phase membership rule) the curated set
holds 28 assemblages; `tests/inference/test_basin_scope.py` pins it and
`analyses/16_basin_membership.py` prints what it drops. The provenance chain
described here is unchanged.

**Status: CLOSED 2026-08-31.** Author ruled option 1 ("you can track that .shp
file"), the share was mounted, and the shapefile set is now in the repository
and verified.

The file was at `//128.226.22.203/LMV`, under
`LMS Survey/Mississippi Site Architecture/Shapefiles/`, dated March 2005. Seven
files copied (`.shp .shx .dbf .prj .sbn .sbx .shp.xml`), 728 KB total, each
md5-verified against the source. No git-lfs needed at that size.

**The decisive check: the membership regenerates byte-identically.** Running
`analyses/16_basin_membership.py` against the copied shapefile reproduces both
committed lists exactly, 29 curated and 92 broad, with `diff` clean and no
working-tree change. The 29-assemblage set is no longer taken on trust; it is
re-derivable from the repository.

**Consequence worth stating: the NAS mount is now a convenience, not a
dependency.** The pipeline runs from a clean clone without it.

## What happened

`analyses/16_basin_membership.py` is the script that defines the study's
analytical unit. It reads `data/Shapefiles/LMVHydrology.shp`, computes distance
from each assemblage to the St. Francis / Tyronza / L'Anguille drainage, and
writes `data/processed/basin_members_curated.txt` (29 assemblages, the
transmission analysis set) and `basin_members_broad.txt` (92 settlement sites).

**`data/Shapefiles/` is not in the repository.** `.gitignore` ignores `data/**`
and allowlists specific files; the shapefile directory is not among them. The
script therefore cannot run from a clean clone.

**It is also not in `data/DATA_INVENTORY.md`.** The inventory lists shapefiles
under `data/CMV/` (GEOFILE, PHILLIPS, phillips1970, PFGOUT, FORD36) and makes no
mention of `LMVHydrology`. So the input to the basin definition is neither
tracked nor inventoried, and a reader of the inventory would not know it is
needed.

## Why it stayed invisible

The two output files **are** allowlisted and committed. Every downstream
analysis reads the committed `.txt`, not the shapefile, so nothing fails. The
artifact is present and the chain that produced it is absent. This is the same
shape as the undeclared-dependency defect found the same day
(`docs/CODE_REVIEW_2026-08-31.md` F14): a committed output masking a broken
generator.

## Why it matters more than a missing file usually would

Rule 19 says not to take a committed result at face value, and this is a
committed result that currently **cannot** be re-derived. The 29-assemblage set
is not an incidental parameter: it is the unit every reported number in the
paper is defined on, including the spatial variance share that carries the
conclusion and the size-controlled F_ST trajectory trend that corroborates it.
(That trend was quoted as +0.01 when this was written; it is +0.004 with a
Monte Carlo interval since `9b8673d`, and the conclusion now leads with the
partition-free spatial model rather than with F_ST.) `docs/METHODS_DECISIONS.md` D-01 and D-02 record the basin
definition as a project choice with a measured reason; that reason is currently
unverifiable from the repository.

It also defeats rule 15 in the specific case that matters most. A `MANIFEST.md`
that lists `16_basin_membership.py` as the entry point for the basin would
today be listing a step that cannot be executed.

## Disposition: option 1, authorized 2026-08-31

**AUTHOR RULING** (Lipo, 2026-08-31): track the shapefile.

Done, to the extent it can be done without the file:
`.gitignore` now carries `!data/Shapefiles/` and
`!data/Shapefiles/LMVHydrology.*`, verified with `git add --dry-run` to admit
every sidecar (`.shp`, `.shx`, `.dbf`, `.prj`) while still ignoring anything
else dropped in that directory. `data/DATA_INVENTORY.md` now lists the file,
which was a gap independent of tracking.

**Resolved 2026-08-31, re-verified 2026-09-04.** The file was copied from the
NAS and committed. `analyses/16_basin_membership.py` runs from the repository
and reproduces both committed lists byte-identically (md5 unchanged, `git
status` clean on `data/processed/`). The "Why it matters" section above records
why this mattered; it no longer describes the current state.

Two things have been added since. The basemap layers the figure scripts need
(`LMVMajorRivers`, `LMVstates`, `LMVcounties4`, `LMVgeology3`) were tracked on
the same reasoning in `9b8673d`, because `make_map.basin_basemap` loads all five
and Figures S4 and S5 cannot be produced without them. And the stricter
watershed subset the main text cites as a robustness check (19 assemblages) is
now derived by the same script rather than asserted in its docstring, and
written to `data/processed/basin_members_watershed.txt` (`0739577`).

### The options as originally written, retained for the record

1. **Track the shapefile.** It is the smallest fix and it restores the chain.
   Add a `!data/Shapefiles/LMVHydrology.*` allowlist, with git-lfs if the file
   is large. Record its source in `DATA_INVENTORY.md`, which is a gap
   independent of whether the file is tracked.
2. **Commit a derived intermediate.** If the shapefile cannot be redistributed,
   commit the per-assemblage distance-to-drainage table, which is small, is the
   only thing the script actually uses from the shapefile, and makes the
   membership reproducible from the repo with the cut clearly visible.
3. **Document the break honestly.** Record in `DATA_INVENTORY.md` and the
   methods supplement that the membership is an external input, and state where
   the shapefile came from so a reader can obtain it.

Option 3 alone is not sufficient for the data-availability statement of a paper
whose conclusion is defined on this unit.

The options above are retained as written for the record. Option 1 was taken;
options 2 and 3 were not needed.
