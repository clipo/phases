#!/usr/bin/env python3
"""Extract the Word comments of a .docx, with the text each comment anchors to.

Written 2026-10-31 to read Beau DiNapoli's annotated v4 main text
(``MAIN_TEXT_beaueditscomments.docx`` at the repository root) into a
markdown file an auditor can read without Word. The extraction is
mechanical: the comment text is reproduced verbatim from
``word/comments.xml``; the anchor is the document text inside the comment's
range; the paragraph is the one containing the range end. Nothing is
paraphrased. Rule 1 of the verification regime asks that process documents
be produced by a committed script, which is why this is here rather than in
a scratchpad.

Usage:
    python3 scripts/extract_docx_comments.py MAIN_TEXT_beaueditscomments.docx \
        > docs/working/2026-10-31-beau-v4-comments-extracted.md

Only the standard library is used.
"""

from __future__ import annotations

import html
import re
import sys
import zipfile
from pathlib import Path


def _text(xml: str) -> str:
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", xml)).split())


def extract(docx: Path) -> str:
    with zipfile.ZipFile(docx) as z:
        doc = z.read("word/document.xml").decode("utf8")
        try:
            com = z.read("word/comments.xml").decode("utf8")
        except KeyError:
            return f"# {docx.name}\n\nNo comments.xml in this document.\n"

    comments: dict[str, tuple[str, str]] = {}
    for m in re.finditer(
        r'<w:comment [^>]*?w:id="(\d+)"[^>]*?(?:w:author="([^"]*)")?[^>]*>(.*?)</w:comment>',
        com,
        re.S,
    ):
        comments[m.group(1)] = (m.group(2) or "unknown author", _text(m.group(3)))

    paragraphs = re.findall(r"<w:p[ >].*?</w:p>", doc, re.S)
    heading = ""
    rows: list[tuple[int, str, str, str]] = []
    for p in paragraphs:
        style = re.search(r'<w:pStyle w:val="([^"]+)"', p)
        ptext = "".join(html.unescape(t) for t in re.findall(r"<w:t[^>]*>(.*?)</w:t>", p, re.S)).strip()
        if style and style.group(1).lower().startswith("heading"):
            heading = ptext
        tokens = re.findall(
            r'<w:commentRangeStart w:id="(\d+)"/>|<w:commentRangeEnd w:id="(\d+)"/>|<w:t[^>]*>(.*?)</w:t>',
            p,
            re.S,
        )
        open_ranges: dict[str, str] = {}
        for start, end, text in tokens:
            if start:
                open_ranges[start] = ""
            elif end:
                anchor = open_ranges.pop(end, "")
                rows.append((int(end), heading, anchor.strip(), ptext))
            else:
                for k in open_ranges:
                    open_ranges[k] += html.unescape(text)
        for k, v in open_ranges.items():
            rows.append((int(k), heading, v.strip() + " [range continues into later paragraphs]", ptext))

    n_ins = len(re.findall(r"<w:ins ", doc))
    n_del = len(re.findall(r"<w:del ", doc))

    out = [
        f"# Comments in `{docx.name}`, extracted mechanically",
        "",
        f"Produced by `scripts/extract_docx_comments.py` from `{docx.name}`. The anchor is "
        "the document text inside each comment's range; the paragraph is the one containing "
        "the range end. Comment text is verbatim from `word/comments.xml`; nothing is "
        "paraphrased. Numbering follows comment order in the file (id + 1).",
        "",
        f"Tracked changes in the document: {n_ins} insertions, {n_del} deletions. "
        f"Comments: {len(comments)}.",
        "",
    ]
    seen: set[int] = set()
    for cid, h, anchor, para in sorted(rows):
        if cid in seen:
            continue
        seen.add(cid)
        author, text = comments.get(str(cid), ("?", "?"))
        out += [
            f"## [{cid + 1}] Section: {h}",
            "",
            f"**Anchor:** {anchor if anchor else '(empty range)'}",
            "",
            f"**Paragraph:** {para[:600]}{'...' if len(para) > 600 else ''}",
            "",
            f"**Comment ({author}):** {text}",
            "",
        ]
    missing = sorted(int(k) for k in comments if int(k) not in seen)
    if missing:
        out += ["## Comments with no range found in the body", ""]
        out += [f"- [{k + 1}] {comments[str(k)][1]}" for k in missing]
        out.append("")
    return "\n".join(out)


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    sys.stdout.write(extract(Path(argv[1])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
