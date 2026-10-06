"""What each pathway means: descriptions kept with every merged IntPath pathway.

A merged pathway carries the descriptions of the source pathways it was merged
from, so its meaning travels with it:

* Reactome: ``pathway2summation.txt`` (summation text);
* WikiPathways: the GPML ``Comment Source="WikiPathways-description"``;
* KEGG: DESCRIPTION and CLASS of the reference map (``get/mapNNNNN``), shared by
  every organism, so ~550 rate-limited calls cover all KEGG organisms;
* GO: the term definition (``def:`` in go-basic.obo).
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

from . import sources


def reactome(raw: Path) -> dict[str, str]:
    f = sources.fetch("https://reactome.org/download/current/pathway2summation.txt",
                      raw / "reactome" / "pathway2summation.txt")
    out = {}
    with open(f, encoding="utf-8", errors="replace") as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2]:
                out[p[0]] = p[2].strip()
    return out


def wikipathways(raw: Path) -> dict[str, str]:
    zips = sorted((raw / "wikipathways").glob("wikipathways-*-gpml-*.zip"))
    out: dict[str, str] = {}
    if not zips:
        return out
    with zipfile.ZipFile(zips[-1]) as zf:
        for name in zf.namelist():
            m = re.search(r"_(WP\d+)_", name)
            if not m:
                continue
            try:
                root = ET.fromstring(zf.read(name))
            except ET.ParseError:
                continue
            ns = root.tag.split("}")[0] + "}" if root.tag.startswith("{") else ""
            for c in root.findall(f"{ns}Comment"):
                if c.get("Source") == "WikiPathways-description" and (c.text or "").strip():
                    out[m.group(1)] = re.sub(r"\s+", " ", c.text).strip()
                    break
    return out


def kegg(shared: Path, pathway_ids: list[str]) -> dict[str, str]:
    """KEGG organism pathway id (hsa04115) -> 'CLASS. DESCRIPTION' of its reference map."""
    out = {}
    for pid in pathway_ids:
        m = re.fullmatch(r"[a-z]+(\d{5})", pid)
        if not m:
            continue
        ref = f"map{m.group(1)}"
        try:
            text = Path(sources.kegg_get(f"get/{ref}", shared / "kegg_maps" / f"{ref}.txt")).read_text(errors="replace")
        except Exception:
            continue
        desc = re.search(r"^DESCRIPTION\s+(.*?)(?=^\S)", text, re.M | re.S)
        cls = re.search(r"^CLASS\s+(.*)$", text, re.M)
        parts = [cls.group(1).strip() + "." if cls else "", re.sub(r"\s+", " ", desc.group(1)).strip() if desc else ""]
        if any(parts):
            out[pid] = " ".join(p for p in parts if p)
    return out


def go_definitions(obo: Path) -> dict[str, str]:
    out, cur = {}, None
    with open(obo, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("id: GO:"):
                cur = line[4:].strip()
            elif line.startswith("def:") and cur:
                m = re.match(r'def:\s+"(.*)"', line)
                if m:
                    out[cur] = m.group(1)
    return out
