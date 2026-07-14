# Building the manuscript Word/PDF files

The `.docx`/`.pdf` outputs are gitignored build products; regenerate them from the
tracked `.md` sources with pandoc. `reference.docx` (tracked) supplies the Word
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
