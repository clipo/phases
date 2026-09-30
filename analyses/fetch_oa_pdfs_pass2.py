"""Second-pass OA fetch via Europe PMC (and arXiv direct) for DOIs missed by Unpaywall.

Europe PMC hosts OA full-text PDFs for most PNAS / Royal Society / eLife / Nature-OA
articles and bypasses publisher bot-blocks. Run after fetch_oa_pdfs.py.
"""
from __future__ import annotations

import json
import os
import time
import urllib.parse
import urllib.request

from fetch_oa_pdfs import DOIS, OUTDIR, UA

EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EPMC_PDF = "https://europepmc.org/articles/{pmcid}?pdf=render"

# arXiv items (DOI not indexed by Unpaywall): key -> arxiv id
ARXIV = {"madsen_2012": "1204.2043"}


def fetch(url: str, timeout: int = 60) -> bytes | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read()
    except Exception:  # noqa: BLE001
        return None


def epmc_ids(doi: str):
    q = urllib.parse.quote(f'DOI:"{doi}"')
    url = f"{EPMC_SEARCH}?query={q}&format=json&resultType=lite"
    raw = fetch(url, timeout=30)
    if not raw:
        return None
    try:
        results = json.loads(raw)["resultList"]["result"]
    except Exception:  # noqa: BLE001
        return None
    # Prefer an explicitly-OA record, but accept any PMC-hosted copy.
    for res in results:
        if res.get("isOpenAccess") == "Y" and res.get("pmcid"):
            return res["pmcid"]
    for res in results:
        if res.get("pmcid"):
            return res["pmcid"]
    return None


def save_pdf(data: bytes, path: str) -> str:
    if not data or data[:4] != b"%PDF":
        return f"not_pdf_{0 if not data else len(data)}b"
    with open(path, "wb") as f:
        f.write(data)
    return f"ok_{len(data)//1024}kb"


def main() -> None:
    got, tried = 0, 0
    for key, doi in DOIS.items():
        path = os.path.join(OUTDIR, f"{key}.pdf")
        if os.path.exists(path):
            continue
        tried += 1
        # arXiv direct
        if key in ARXIV:
            data = fetch(f"https://arxiv.org/pdf/{ARXIV[key]}")
            msg = save_pdf(data, path) if data else "arxiv_fail"
            print(f"{('OK' if msg.startswith('ok') else 'FAIL'):5s} {key:32s} arxiv {msg}")
            if msg.startswith("ok"):
                got += 1
            continue
        pmcid = epmc_ids(doi)
        if not pmcid:
            print(f"FAIL  {key:32s} no_epmc_pmc")
            continue
        data = fetch(EPMC_PDF.format(pmcid=pmcid))
        msg = save_pdf(data, path) if data else "epmc_fetch_fail"
        print(f"{('OK' if msg.startswith('ok') else 'FAIL'):5s} {key:32s} {pmcid} {msg}")
        if msg.startswith("ok"):
            got += 1
        time.sleep(0.4)
    print(f"\nPass 2: recovered {got} of {tried} still-missing.")


if __name__ == "__main__":
    main()
