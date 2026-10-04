# September 2026 revision

`MAIN_TEXT.md` and `SUPPLEMENTAL_TEXT.md` are edited directly (the template
assembly used on 2026-09-09 was retired). Regenerate figures with analyses 21,
47, 48, 50, and 29 before building; then use the pandoc commands below.

# Building the manuscript Word/PDF files

The `.pdf` outputs are gitignored build products. The `.docx` outputs ARE
tracked, because docx is the format American Antiquity requires and the file
the coauthors receive; regenerate both from the tracked `.md` sources with
pandoc. (This paragraph said "`.docx`/`.pdf` ... are gitignored" until
2026-09-21, which was wrong about the docx.) `reference.docx` (tracked) supplies the Word
styling: **Times New Roman, 12 pt, black, with page numbers** (per American
Antiquity manuscript format). Do not drop `--reference-doc`, or the Word files
revert to pandoc's default Aptos/theme fonts with no page numbers.

From `docs/manuscript/`:

```bash
# Main text (docx + pdf)
pandoc MAIN_TEXT.md --citeproc --bibliography=references.bib \
  --csl=american-antiquity.csl --reference-doc=reference.docx -o MAIN_TEXT.docx
pandoc MAIN_TEXT.md --citeproc --bibliography=references.bib \
  --csl=american-antiquity.csl --include-in-header=_pdf_header.tex -o MAIN_TEXT.pdf

# Supplemental (docx + pdf)
pandoc SUPPLEMENTAL_TEXT.md --citeproc --bibliography=references.bib \
  --csl=american-antiquity.csl --reference-doc=reference.docx -o SUPPLEMENTAL_TEXT.docx
pandoc SUPPLEMENTAL_TEXT.md --citeproc --bibliography=references.bib \
  --csl=american-antiquity.csl --include-in-header=_pdf_header.tex -o SUPPLEMENTAL_TEXT.pdf

# Cover letter (docx)
pandoc COVER_LETTER.md --reference-doc=reference.docx -o COVER_LETTER.docx
```

`reference.docx` was built from `pandoc --print-default-data-file reference.docx`
with the theme Latin fonts set to Times New Roman, `docDefaults` forced to
Times New Roman / 12 pt / black, every style size normalized to 12 pt, and a
centered `PAGE`-field footer added.

# Keeping Google Drive current

The coauthors read and edit the paper in the Drive folder
`Explaining Phase Structure - Parkin Paper`, so it drifts from the repo in both
directions. `scripts/sync_drive.sh` builds the docx/pdf and pushes the paper --
both `.md` files, the four build products, `references.bib`, this file, and
exactly the figures the two manuscripts reference, read out of them rather than
listed anywhere. It verifies by checksum after pushing, because a copy that
silently did nothing also exits 0.

```bash
scripts/sync_drive.sh --check     # report differences, change nothing
scripts/sync_drive.sh             # build, push, verify
scripts/sync_drive.sh --force     # push although the copies differ
```

It stops when the markdown differs, shows the diff, and makes you say which
side is newer, since a checksum cannot. The working habit is to `rclone cat`
the Drive copy first, edit that, and write both.

On 2026-09-21 the Drive `figures/` folder was four days stale and did not hold
Figure 12 at all, while the `.md` files were current: the manuscript on Drive
referenced a figure that was not there. That is what this script exists to stop.
