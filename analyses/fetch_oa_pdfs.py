"""Fetch open-access PDFs for the lit-review DOIs via Unpaywall.

Unpaywall indexes only legally hosted OA copies. Paywalled items resolve to
no OA location and are reported for institutional retrieval. Run from repo root:
    python analyses/fetch_oa_pdfs.py
"""
from __future__ import annotations

import json
import os
import time
import urllib.request
import urllib.error

EMAIL = "clipo@binghamton.edu"
OUTDIR = "docs/references/pdfs"
UA = "Mozilla/5.0 (mls-emergence reference pipeline; mailto:clipo@binghamton.edu)"

# key (author_year, pipeline filename stem) -> DOI
DOIS: dict[str, str] = {
    # Strand A: transmission / neutral-theory methods
    "neiman_1995": "10.2307/282074",
    "shennan_wilkinson_2001": "10.2307/2694174",
    "steele_glatz_kandler_2010": "10.1016/j.jas.2009.12.039",
    "madsen_2012": "10.48550/arxiv.1204.2043",
    "premo_2014": "10.1086/674873",
    "bentley_hahn_shennan_2004": "10.1098/rspb.2004.2746",
    "bentley_shennan_2003": "10.2307/3557104",
    "bentley_shennan_2005": "10.1126/science.309.5736.877",
    "mesoudi_lycett_2009": "10.1016/j.evolhumbehav.2008.07.005",
    "brantingham_perreault_2010": "10.1016/j.jas.2010.07.021",
    "bentley_2008": "10.1371/journal.pone.0003057",
    "acerbi_bentley_2014": "10.1016/j.evolhumbehav.2014.02.003",
    "carrignon_bentley_ruck_2019": "10.1057/s41599-019-0295-9",
    "kandler_shennan_2013": "10.1016/j.jtbi.2013.03.006",
    "crema_kandler_shennan_2016": "10.1038/srep39122",
    "kandler_shennan_2015": "10.1098/rsif.2015.0905",
    "crema_edinborough_2014": "10.1016/j.jas.2014.07.014",
    "kandler_wilder_fortunato_2017": "10.1101/111575",
    "newberry_plotkin_2022": "10.1038/s41562-022-01342-6",
    "carrignon_bentley_obrien_2023": "10.1016/j.jaa.2023.101545",
    "crema_bortolini_lake_2023": "10.1007/s10816-022-09599-x",
    "fogarty_kandler_creanza_2024": "10.1073/pnas.2418106121",
    "mesoudi_2021": "10.1098/rstb.2020.0053",
    "lipo_madsen_dunnell_2015": "10.1371/journal.pone.0124942",
    "madsen_lipo_2014": "10.6084/m9.figshare.1269322.v2",
    "lyman_2019": "10.1016/j.shpsc.2019.101178",
    "porcic_2018": "10.1002/9781119188230.saseas0531",
    "bell_richerson_mcelreath_2009": "10.1073/pnas.0903232106",
    "richerson_etal_2016": "10.1017/s0140525x1400106x",
    "riede_hoggard_shennan_2019": "10.1057/s41599-019-0260-7",
    "derex_perreault_boyd_2018": "10.1098/rstb.2017.0062",
    "premo_2020": "10.1007/s10816-020-09464-9",
    # Strand B: central MS Valley / chronology
    "cobb_legg_smith_2021": "10.1017/aaq.2021.17",
    "cobb_butler_2002": "10.2307/1593795",
    "cobb_krus_2023": "10.1007/s10816-023-09613-w",
    "cobb_krus_2025": "10.1017/aaq.2025.27",
    "sullivan_2024": "10.1017/aaq.2024.1",
    "alvey_2019": "10.32469/10355/70059",
    "pompeani_2021": "10.1038/s41598-021-92900-x",
    "white_2020": "10.1017/aaq.2019.103",
    "bird_2017": "10.1038/srep41628",
    "krus_2023": "10.1016/j.jaa.2023.101486",
    "kessler_2022": "10.1017/rdc.2022.84",
    "holland-lulewicz_ritchison_2021": "10.1017/aap.2021.10",
    "holland-lulewicz_2019": "10.1073/pnas.1818346116",
    "thompson_2022": "10.1017/aaq.2022.31",
    "gavrilets_anderson_turchin_2010": "10.21237/c7clio11193",
    "mehta_2012": "10.1155/2012/192923",
    "lankford_house_2006": "10.2307/40038297",
    "ritchison_2025": "10.1017/aaq.2025.22",
    "ogorman_conner_2023": "10.1017/aaq.2022.95",
    "bird_2022_p3k14c": "10.1038/s41597-022-01118-7",
    # Strand C: MLS / assortment / transitions
    "price_1972": "10.1111/j.1469-1809.1957.tb01874.x",
    "gardner_2009": "10.1098/rstb.2009.0108",
    "frank_2012": "10.1111/j.1420-9101.2012.02498.x",
    "vanveelen_2012": "10.1016/j.jtbi.2011.07.025",
    "okasha_otsuka_2020": "10.1098/rstb.2019.0365",
    "birch_2019": "10.1016/j.cub.2019.01.065",
    "birch_okasha_2023": "10.1007/s11229-023-04285-1",
    "riolo_cohen_axelrod_2001": "10.1038/35106555",
    "roberts_sherratt_2002": "10.1038/418499b",
    "jansen_vanbaalen_2006": "10.1038/nature04387",
    "traulsen_nowak_2007": "10.1371/journal.pone.0000270",
    "gardner_west_2010": "10.1111/j.1558-5646.2009.00842.x",
    "west_gardner_2010": "10.1126/science.1178332",
    "biernaskie_perry_west_2019": "10.1016/j.tree.2019.08.001",
    "nowak_2006": "10.1126/science.1133755",
    "noe_hammerstein_1994": "10.1007/BF00167053",
    "barclay_2013": "10.1016/j.evolhumbehav.2013.02.002",
    "baumard_andre_sperber_2013": "10.1017/S0140525X11002202",
    "geffroy_andre_2017": "10.1098/rspb.2016.2317",
    "hauert_2002_science": "10.1126/science.1070582",
    "szabo_hauert_2002": "10.1103/PhysRevLett.89.118101",
    "garcia_traulsen_2012": "10.1016/j.jtbi.2012.05.011",
    "powers_2016": "10.1098/rstb.2015.0098",
    "powers_2023": "10.1177/17456916231179156",
    "richerson_boyd_2009": "10.1098/rstb.2009.0134",
    "henrich_2004": "10.1016/S0167-2681(03)00094-5",
    "mathew_boyd_2011": "10.1073/pnas.1105604108",
    "choi_bowles_2007": "10.1126/science.1144237",
    "bowles_2009": "10.1126/science.1168112",
    "smaldino_2014": "10.1017/S0140525X13001544",
    "waring_2017": "10.1007/s11625-017-0501-x",
    "michod_1997": "10.1086/286012",
    "queller_strassmann_2009": "10.1098/rstb.2009.0095",
    "west_fisher_gardner_kiers_2015": "10.1073/pnas.1421402112",
    "bourrat_2022": "10.7554/eLife.73715",
    "black_bourrat_rainey_2020": "10.1038/s41559-019-1086-9",
    "hammerschmidt_rainey_2020": "10.1038/s41586-020-2653-6",
    "west_cooper_2021": "10.1038/s41559-020-01384-x",
    "boomsma_gawne_2017": "10.1111/brv.12330",
    "gardner_grafen_2009": "10.1111/j.1420-9101.2008.01681.x",
}


def get_json(url: str, timeout: int = 25):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def oa_pdf_url(doi: str) -> tuple[str | None, str | None, str | None]:
    """Return (pdf_url, host_type, license) from Unpaywall, or (None, status, None)."""
    url = f"https://api.unpaywall.org/v2/{doi}?email={EMAIL}"
    try:
        data = get_json(url)
    except urllib.error.HTTPError as e:
        return None, f"unpaywall_http_{e.code}", None
    except Exception as e:  # noqa: BLE001
        return None, f"unpaywall_err_{type(e).__name__}", None
    if not data.get("is_oa"):
        return None, "closed", None
    locs = []
    if data.get("best_oa_location"):
        locs.append(data["best_oa_location"])
    locs.extend(data.get("oa_locations") or [])
    for loc in locs:
        purl = loc.get("url_for_pdf") or (loc.get("url") if loc.get("url", "").endswith(".pdf") else None)
        if purl:
            return purl, loc.get("host_type"), loc.get("license")
    # OA landing page exists but no direct PDF link
    bl = data.get("best_oa_location") or {}
    return None, "oa_no_pdf_link", bl.get("url_for_landing_page")


def download(url: str, path: str, timeout: int = 60) -> tuple[bool, str]:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/pdf,*/*"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
    except Exception as e:  # noqa: BLE001
        return False, f"download_err_{type(e).__name__}"
    if data[:4] != b"%PDF":
        # some hosts return an HTML interstitial
        return False, f"not_pdf_{len(data)}b"
    with open(path, "wb") as f:
        f.write(data)
    return True, f"ok_{len(data)//1024}kb"


def main() -> None:
    os.makedirs(OUTDIR, exist_ok=True)
    results = []
    for key, doi in DOIS.items():
        path = os.path.join(OUTDIR, f"{key}.pdf")
        if os.path.exists(path):
            results.append((key, doi, "skip_exists", ""))
            continue
        purl, status, extra = oa_pdf_url(doi)
        if purl:
            ok, msg = download(purl, path)
            results.append((key, doi, "OA_DOWNLOADED" if ok else f"OA_FAIL:{msg}", purl))
        else:
            results.append((key, doi, status, extra or ""))
        time.sleep(0.5)

    # manifest
    man = os.path.join(OUTDIR, "_download_manifest.md")
    got = sum(1 for r in results if r[2] == "OA_DOWNLOADED" or r[2] == "skip_exists")
    with open(man, "w") as f:
        f.write("# OA download manifest\n\n")
        f.write(f"Downloaded/present: {got} of {len(results)} DOIs. "
                "Items marked `closed`, `oa_no_pdf_link`, or `OA_FAIL` need institutional access "
                "or manual retrieval.\n\n")
        f.write("| key | doi | status | url/landing |\n|---|---|---|---|\n")
        for key, doi, status, url in sorted(results, key=lambda x: x[2]):
            f.write(f"| {key} | {doi} | {status} | {url} |\n")
    print(f"Done. {got}/{len(results)} available. Manifest: {man}")
    for key, doi, status, _ in results:
        print(f"{status:18s} {key}")


if __name__ == "__main__":
    main()
