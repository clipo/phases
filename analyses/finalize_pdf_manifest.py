"""Retry two known PMCIDs, then write present/missing manifests.

Produces docs/references/pdfs/_to_acquire_manually.md listing every DOI without a
local PDF, so they can be fetched manually and dropped into the pdfs/ directory.
"""
from __future__ import annotations

import os
import time
import urllib.request

from fetch_oa_pdfs import DOIS, OUTDIR, UA

# Tier-1 high-priority keys (from docs/references/acquisition-list.md)
TIER1 = {
    "smaldino_2014", "gardner_grafen_2009", "bourrat_2022", "vanveelen_2012",
    "west_fisher_gardner_kiers_2015", "powers_2016", "hauert_2002_science",
    "biernaskie_perry_west_2019", "neiman_1995", "crema_kandler_shennan_2016",
    "kandler_wilder_fortunato_2017", "lipo_madsen_dunnell_2015",
    "shennan_wilkinson_2001", "bell_richerson_mcelreath_2009", "sullivan_2024",
    "cobb_krus_2025", "holland-lulewicz_2019", "thompson_2022", "bird_2022_p3k14c",
}

RETRY_PMC = {"frank_2012": "PMC3354028", "nowak_2006": "PMC3279745"}


def fetch(url: str, timeout: int = 60) -> bytes | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except Exception:  # noqa: BLE001
        return None


def main() -> None:
    for key, pmcid in RETRY_PMC.items():
        path = os.path.join(OUTDIR, f"{key}.pdf")
        if os.path.exists(path):
            continue
        for _ in range(3):
            data = fetch(f"https://europepmc.org/articles/{pmcid}?pdf=render")
            if data and data[:4] == b"%PDF":
                with open(path, "wb") as f:
                    f.write(data)
                print(f"OK retry {key} ({len(data)//1024}kb)")
                break
            time.sleep(1.0)
        else:
            print(f"still failing: {key}")

    present, missing = [], []
    for key, doi in DOIS.items():
        if os.path.exists(os.path.join(OUTDIR, f"{key}.pdf")):
            present.append((key, doi))
        else:
            missing.append((key, doi))

    out = os.path.join(OUTDIR, "_to_acquire_manually.md")
    with open(out, "w") as f:
        f.write("# PDFs to acquire manually\n\n")
        f.write(f"Present locally: {len(present)} / {len(DOIS)}. "
                f"Missing (no legal OA copy found): {len(missing)}.\n\n")
        f.write("Fetch each via institutional access and save as "
                "`docs/references/pdfs/<key>.pdf` (use the key in column 1).\n\n")
        f.write("## High priority (Tier 1)\n\n| key | DOI | link |\n|---|---|---|\n")
        for key, doi in missing:
            if key in TIER1:
                f.write(f"| {key} | {doi} | https://doi.org/{doi} |\n")
        f.write("\n## Remaining\n\n| key | DOI | link |\n|---|---|---|\n")
        for key, doi in missing:
            if key not in TIER1:
                f.write(f"| {key} | {doi} | https://doi.org/{doi} |\n")
        f.write("\n## Already present\n\n")
        for key, doi in sorted(present):
            f.write(f"- {key} ({doi})\n")
    print(f"\nPresent {len(present)}, missing {len(missing)}. List: {out}")


if __name__ == "__main__":
    main()
