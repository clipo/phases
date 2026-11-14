"""Build the submission figure package for American Antiquity.

Reads the figure files the manuscripts embed (figures/*.png and their vector
companions) and writes, under docs/manuscript/submission/figures/, one TIFF
and one PDF per figure named as the Cambridge Journals Artwork Guide asks
([First author surname]FigN), plus MANIFEST.md listing each file's source,
producing script, pixel size, printed width at 300 dpi and SHA-256.

Formats (Cambridge Journals Artwork Guide, read 2026-11-01): TIFF with LZW
compression for halftones at 300 dpi or more; PDF with embedded fonts
recommended for line and combination artwork. The TIFF is the matplotlib
render the manuscript shows, flattened onto white (the PNGs carry an alpha
channel) and tagged at 300 dpi. The PDF is matplotlib's own vector export
where the figure script wrote one; for the three figures whose scripts write
only SVG, the PDF is converted from the SVG with rsvg-convert. Matplotlib
writes the SVG text as paths, so the PDF needs no font embedding, and the
rsvg render of the SVG was checked pixel-identical to the PNG for Figure 1.

Colour is left in RGB. The guide asks for CMYK; a numerical RGB-to-CMYK
conversion shifts the Okabe-Ito hues the figures depend on, so that
conversion is left to the press, and the manifest says so.

No analysis is run here. The figures are produced by the scripts named in
FIGURES; regenerate them first if their findings files changed.

Usage: .venv/bin/python scripts/build_submission_figures.py
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
OUT_DIR = ROOT / "docs" / "manuscript" / "submission" / "figures"
SURNAME = "Lipo"
DPI = 300

# (label, source stem, producing script, manuscript, short title)
FIGURES = [
    ("Fig1", "fig1_phillips1970_phases", "analyses/93_phillips1970_phase_map.py", "Main text",
     "Study area and the late pre-contact phases as Phillips (1970) defined them"),
    ("Fig2", "fig2_two_histories", "analyses/102_fig2_two_histories.py", "Main text",
     "Two histories, one record (illustrative schematic)"),
    ("Fig3", "fig3_model_schematic", "analyses/103_fig3_model_schematic.py", "Main text",
     "The three models on the site locations (schematic)"),
    ("Fig4", "fig4_phase_recovery_map", "analyses/104_fig4_phase_recovery_map.py", "Main text",
     "What recovers the phases"),
    ("Fig5", "fig5_spatial_field", "analyses/105_fig5_spatial_field.py", "Main text",
     "The spatial field"),
    ("Fig6", "fig6_per_class_drift", "analyses/106_fig6_per_class_drift.py", "Main text",
     "Common decorated classes against neutral copying with regional drift"),
    ("Fig7", "fig7_rim_attributes", "analyses/107_fig7_rim_attributes.py", "Main text",
     "Rim attributes beside decoration on the 19 sites with rim tabulations"),
    ("Fig8", "fig8_latent_unit_recovery", "analyses/108_fig8_latent_unit_recovery.py", "Main text",
     "Recovery of known structure on this site geometry"),
    ("FigS1", "figS1_chronology", "analyses/11_chronology_14c.py", "Supplemental Text",
     "Radiocarbon chronology of the St. Francis basin"),
    ("FigS2", "figS2_ranksize", "analyses/make_figures.py", "Supplemental Text",
     "Settlement cross-check (mound height, hoe exchange)"),
]


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def write_tiff(src: Path, dst: Path) -> tuple[int, int]:
    im = Image.open(src)
    if im.mode in ("RGBA", "LA", "P"):
        rgba = im.convert("RGBA")
        flat = Image.new("RGB", rgba.size, (255, 255, 255))
        flat.paste(rgba, mask=rgba.split()[3])
        im = flat
    else:
        im = im.convert("RGB")
    im.save(dst, format="TIFF", compression="tiff_lzw", dpi=(DPI, DPI))
    return im.size


def write_pdf(stem: str, dst: Path) -> str:
    mpl_pdf = FIG_DIR / f"{stem}.pdf"
    svg = FIG_DIR / f"{stem}.svg"
    if mpl_pdf.is_file():
        shutil.copyfile(mpl_pdf, dst)
        return "matplotlib PDF export"
    if not svg.is_file():
        raise FileNotFoundError(f"no PDF or SVG for {stem}")
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(dst), str(svg)], check=True)
    return "rsvg-convert from the matplotlib SVG"


def main() -> None:
    if OUT_DIR.exists():
        for old in OUT_DIR.glob(f"{SURNAME}Fig*"):
            old.unlink()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, stem, script, where, title in FIGURES:
        png = FIG_DIR / f"{stem}.png"
        if not png.is_file():
            raise FileNotFoundError(png)
        tif = OUT_DIR / f"{SURNAME}{label}.tif"
        pdf = OUT_DIR / f"{SURNAME}{label}.pdf"
        w, h = write_tiff(png, tif)
        pdf_src = write_pdf(stem, pdf)
        mtime = datetime.fromtimestamp(png.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        rows.append((label, where, title, stem, script, mtime, w, h, tif, pdf, pdf_src))

    n_main = sum(1 for row in FIGURES if row[3] == "Main text")
    lines = [
        "# Submission figures",
        "",
        f"Built {datetime.now().strftime('%Y-%m-%d %H:%M')} by `scripts/build_submission_figures.py` "
        "from the figure files the manuscripts embed. Rebuild after any figure script runs.",
        "",
        "Each figure is supplied twice. The TIFF (LZW, RGB, 300 dpi) is the exact render shown in "
        "the manuscript, flattened onto white. The PDF is the vector version, with text as paths, "
        "and is the file to prefer for line and combination artwork at any resolution. Colour is "
        "RGB; the Cambridge Journals Artwork Guide asks for CMYK, and the conversion is left to "
        "the press because a numerical conversion shifts the Okabe-Ito palette. Captions are in "
        "the manuscript text, not in the files.",
        "",
        f"Main-text figures: {n_main} (journal maximum 10). Supplemental figures S1 and S2 are embedded "
        "in SUPPLEMENTAL_TEXT.docx, which the journal posts as supplied; the files here are for "
        "the record.",
        "",
        "| Figure | Manuscript | Short title | Source figure file | Producing script | Source written | Pixels (w x h) | Width at 300 dpi (in) | TIFF | PDF (origin) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for label, where, title, stem, script, mtime, w, h, tif, pdf, pdf_src in rows:
        lines.append(
            f"| {label.replace('Fig', 'Figure ')} | {where} | {title} | `figures/{stem}.png` | `{script}` | {mtime} "
            f"| {w} x {h} | {w / DPI:.2f} | `{tif.name}` ({tif.stat().st_size // 1024} KB, sha256 {sha256(tif)}) "
            f"| `{pdf.name}` ({pdf.stat().st_size // 1024} KB, sha256 {sha256(pdf)}; {pdf_src}) |"
        )
    (OUT_DIR / "MANIFEST.md").write_text("\n".join(lines) + "\n")
    for r in rows:
        print(f"{r[0]:6s} {r[6]}x{r[7]} {r[6]/DPI:.2f} in  {r[8].name}  {r[9].name}")
    print(f"wrote {OUT_DIR / 'MANIFEST.md'}")


if __name__ == "__main__":
    main()
