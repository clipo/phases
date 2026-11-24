#!/usr/bin/env python3
"""Build the submission Word and PDF files and the upload archive.

Writes, under docs/manuscript/submission/:

    LMV-Phases.docx / .pdf                      from MAIN_TEXT.md
    LMV-Phases_Supplemental_Text.docx / .pdf    from SUPPLEMENTAL_TEXT.md
    COVER_LETTER.docx                           from COVER_LETTER.md

and docs/manuscript/submission.zip, holding the submission folder (the two
manuscripts in both formats, the cover letter and the figure package).

Each Word file is the pandoc build described in docs/manuscript/BUILD.md
(Times New Roman, 12 pt, page numbers, from reference.docx) with one change:
the text inside tables is set smaller, 9 pt in tables of six or more columns
and 10 pt in the others. Those are the sizes the hand-styled submission files
of 2026-11-19 used; at 12 pt the wider tables break words across lines. No
text is changed, and the script checks that by comparing the plain text of
each file before and after the change.

The PDFs are converted from the Word files with LibreOffice, so they follow
the Word layout.

The cover letter is built from COVER_LETTER.md, which is the authoritative
copy by author ruling of 2026-11-22; edits made to the Word letter by hand are
overwritten on the next run. It has no tables and no citations, so it gets no
table sizing and no bibliography. The figure package is built by
scripts/build_submission_figures.py; run that first if a figure changed.

Usage: .venv/bin/python scripts/build_submission_docs.py
Requires pandoc and soffice (LibreOffice) on the PATH.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "docs" / "manuscript"
SUB = MS / "submission"
ARCHIVE = MS / "submission.zip"
DOCS = {"MAIN_TEXT.md": "LMV-Phases", "SUPPLEMENTAL_TEXT.md": "LMV-Phases_Supplemental_Text"}
WIDE_COLUMNS = 6          # tables with at least this many columns take the smaller size
HALF_POINTS_WIDE = 18     # 9 pt
HALF_POINTS_OTHER = 20    # 10 pt
SKIP = {".DS_Store"}


def skipped(name: str) -> bool:
    """Files in the submission folder that are not part of the upload archive.

    "~$" files are the lock files Word leaves beside an open document. "AQ-S-*.pdf"
    is the PDF the journal's submission system builds from the uploaded files.
    """
    return name in SKIP or name.startswith("~$") or (name.startswith("AQ-S-") and name.endswith(".pdf"))


def run(cmd: list[str], cwd: Path, strict: bool = False) -> None:
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed:\n{r.stdout}\n{r.stderr}")
    if strict and r.stderr.strip():           # pandoc reports an unresolved citation here
        raise RuntimeError(f"{cmd[0]} printed a warning, which is treated as an error:\n{r.stderr}")


def size_table_text(xml: str) -> tuple[str, list[tuple[int, int]]]:
    """Set the run size inside every table; returns the new XML and (columns, half-points) per table."""
    done: list[tuple[int, int]] = []

    def one_table(m: re.Match) -> str:
        tbl = m.group(0)
        if tbl.count("<w:tbl>") != 1:
            raise ValueError("a nested table was found; this script sizes only flat tables")
        cols = len(re.findall(r"<w:gridCol\b", tbl))
        if cols == 0:
            raise ValueError("a table without a column grid was found")
        hp = HALF_POINTS_WIDE if cols >= WIDE_COLUMNS else HALF_POINTS_OTHER
        size = f'<w:sz w:val="{hp}" /><w:szCs w:val="{hp}" />'

        def one_run(r: re.Match) -> str:
            run_xml = r.group(0)
            if "<w:sz " in run_xml:
                raise ValueError("a table run already carries a size; refusing to guess which to keep")
            if "<w:rPr>" in run_xml:
                return run_xml.replace("<w:rPr>", "<w:rPr>" + size, 1)
            return re.sub(r"(<w:r(?: [^>]*)?>)", r"\1<w:rPr>" + size + "</w:rPr>", run_xml, count=1)

        done.append((cols, hp))
        return re.sub(r"<w:r(?: [^>]*)?>.*?</w:r>", one_run, tbl, flags=re.S)

    return re.sub(r"<w:tbl>.*?</w:tbl>", one_table, xml, flags=re.S), done


def plain(docx: Path) -> str:
    r = subprocess.run(["pandoc", str(docx), "-t", "plain", "--wrap=none"], capture_output=True, text=True, check=True)
    return r.stdout


def build(md: str, stem: str) -> dict:
    out = SUB / f"{stem}.docx"
    with tempfile.TemporaryDirectory() as tmp:
        raw = Path(tmp) / "raw.docx"
        run(["pandoc", md, "--citeproc", "--bibliography=references.bib", "--csl=american-antiquity.csl",
             "--reference-doc=reference.docx", "-o", str(raw)], cwd=MS, strict=True)
        with zipfile.ZipFile(raw) as z:
            members = [(i, z.read(i.filename)) for i in z.infolist()]
        tables: list[tuple[int, int]] = []
        styled = Path(tmp) / "styled.docx"
        with zipfile.ZipFile(styled, "w", zipfile.ZIP_DEFLATED) as z:
            for info, data in members:
                if info.filename == "word/document.xml":
                    text, tables = size_table_text(data.decode("utf-8"))
                    data = text.encode("utf-8")
                z.writestr(info, data)
        if plain(raw) != plain(styled):
            raise RuntimeError(f"{stem}: the table sizing changed the text of the document")
        shutil.copyfile(styled, out)
    pdf = SUB / f"{stem}.pdf"
    pdf.unlink(missing_ok=True)
    run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", str(SUB), str(out)], cwd=SUB)
    if not pdf.is_file():
        raise RuntimeError(f"{pdf} was not written")
    return {"stem": stem, "words": len(plain(out).split()), "tables": tables}


def build_letter() -> int:
    out = SUB / "COVER_LETTER.docx"
    run(["pandoc", "COVER_LETTER.md", "--reference-doc=reference.docx", "-o", str(out)], cwd=MS, strict=True)
    return len(plain(out).split())


def archive() -> int:
    files = sorted(p for p in SUB.rglob("*") if p.is_file() and not skipped(p.name))
    ARCHIVE.unlink(missing_ok=True)
    with zipfile.ZipFile(ARCHIVE, "w", zipfile.ZIP_DEFLATED) as z:
        for p in files:
            z.write(p, Path("submission") / p.relative_to(SUB))
    return len(files)


def main() -> int:
    for tool in ("pandoc", "soffice"):
        if shutil.which(tool) is None:
            print(f"{tool} is not on the PATH", file=sys.stderr)
            return 2
    for md, stem in DOCS.items():
        r = build(md, stem)
        print(f"{r['stem']}.docx and .pdf written: {r['words']} words by whitespace count of the plain text; "
              f"tables (columns, half-points): {r['tables']}")
    if (MS / "COVER_LETTER.md").is_file():
        print(f"COVER_LETTER.docx written: {build_letter()} words by whitespace count of the plain text")
    else:                                     # the public release does not carry the letter
        print("COVER_LETTER.md is absent, so no cover letter was built")
    print(f"{ARCHIVE.relative_to(ROOT)} written with {archive()} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
