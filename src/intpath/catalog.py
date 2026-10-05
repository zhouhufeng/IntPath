"""Organism catalog for building IntPathV2 across every KEGG organism.

The hand-curated registry (:mod:`intpath.organisms`) covers the key organisms
with every source (Reactome, WikiPathways, MSigDB, GO Consortium, HuRI, ...).
This module generates :class:`Organism` entries for all other KEGG genomes, so
the same builder runs on thousands of organisms with what exists for them:

* KEGG pathway membership (licensed: full tier only);
* GO from UniProt cross-references by taxid;
* STRING v12 physical interactions when STRING has the taxid, plus BioGRID /
  IntAct where those databases have the organism;
* gene identifiers from NCBI gene_info (filtered from the shared group files),
  or, where NCBI Gene has no records (most prokaryotes), from KEGG's own gene list.

KEGG REST is rate-limited (max 3 calls/s, academic use only): every KEGG call
in this module goes through :func:`kegg_get`, which serialises and spaces them.
"""

from __future__ import annotations

import gzip
import json
import re
from dataclasses import asdict
from pathlib import Path

from .organisms import ORGANISMS, Organism
from .sources import fetch, fetch_text, kegg_get

STRING_SPECIES = "https://stringdb-downloads.org/download/species.v12.0.txt"
NCBI_GROUPS = ("Mammalia", "Non-mammalian_vertebrates", "Invertebrates", "Plants", "Fungi", "Protozoa",
               "Archaea_Bacteria")

_VERTEBRATES = {"birds", "reptiles", "amphibians", "fishes", "other vertebrates"}


def _ncbi_group(lineage: str) -> str:
    """KEGG lineage (br08601, e.g. "Eukaryotes;Animals;Birds;...") -> NCBI GENE_INFO group directory."""
    parts = [p.strip().lower() for p in lineage.split(";") if p.strip()]
    if not parts or parts[0] == "prokaryotes":
        return "Archaea_Bacteria"
    if "mammals" in parts:
        return "Mammalia"
    if _VERTEBRATES & set(parts):
        return "Non-mammalian_vertebrates"
    if "animals" in parts:
        return "Invertebrates"
    if "plants" in parts:
        return "Plants"
    if "fungi" in parts:
        return "Fungi"
    return "Protozoa"


def build_catalog(shared: Path) -> list[dict]:
    """All KEGG organisms with taxid, lineage and STRING availability (cached as catalog.json).

    * codes and T numbers: ``list/genome`` ("T01001<TAB>hsa; Homo sapiens (human)");
    * lineage: the KEGG Organism BRITE hierarchy ``br:br08601``;
    * taxid: the TAXONOMY line of each genome entry, fetched 10 entries per call.
    (``list/organism`` used to carry lineage but now answers HTTP 400.)
    """
    out = shared / "catalog" / "catalog.json"
    if out.exists():
        return json.loads(out.read_text())
    cat = shared / "catalog"
    cat.mkdir(parents=True, exist_ok=True)
    genomes = kegg_get("list/genome", cat / "kegg_list_genome.tsv")
    entries: dict[str, dict] = {}
    for line in open(genomes):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or ";" not in f[1]:
            continue
        code, name = (x.strip() for x in f[1].split(";", 1))
        if code and "," not in code:
            entries[code] = {"kegg": code, "t_number": f[0].replace("gn:", ""), "name": name}

    lineage: dict[str, str] = {}
    path: list[str] = []
    for line in open(kegg_get("get/br:br08601", cat / "br08601.keg")):
        if not line[:1].isalpha() or line[:1] in "+!#":
            continue
        level, text = ord(line[0]) - ord("A"), line[1:].strip()
        m = re.match(r"^([a-z]{3,4})\s{2,}(.+)$", text)
        if m and m.group(1) in entries:
            lineage[m.group(1)] = ";".join(path[:level])
        else:
            path[level:] = [re.sub(r"\s*\(\d+\)$", "", text)]

    taxid: dict[str, str] = {}
    tnums = [e["t_number"] for e in entries.values()]
    by_t = {e["t_number"]: e["kegg"] for e in entries.values()}
    from concurrent.futures import ThreadPoolExecutor

    def fetch_chunk(chunk: list[str]) -> str:
        dest = cat / "genome_entries" / f"{chunk[0]}_{len(chunk)}.txt"
        try:
            return Path(kegg_get("get/" + "+".join(f"gn:{t}" for t in chunk), dest)).read_text(errors="replace")
        except Exception:
            return ""  # retried on the next catalog build (the cache file is absent)

    chunks = [tnums[i:i + 10] for i in range(0, len(tnums), 10)]
    with ThreadPoolExecutor(max_workers=3) as pool:  # kegg_get keeps the overall rate within KEGG's limit
        texts = list(pool.map(fetch_chunk, chunks))
    for text in texts:
        for block in text.split("///"):
            t = re.search(r"^ENTRY\s+(T\d+)", block, re.M)
            x = re.search(r"^TAXONOMY\s+TAX:(\d+)", block, re.M)
            if t and x and t.group(1) in by_t:
                taxid[by_t[t.group(1)]] = x.group(1)

    string_taxa = {line.split("\t")[0] for line in fetch_text(STRING_SPECIES).splitlines()[1:]}
    rows = []
    for code, e in entries.items():
        lin = lineage.get(code, "")
        tx = taxid.get(code, "")
        rows.append({**e, "lineage": lin, "taxid": tx, "ncbi_group": _ncbi_group(lin or "Prokaryotes"),
                     "string": tx in string_taxa})
    out.write_text(json.dumps(rows))
    return rows


def organisms_from_catalog(rows: list[dict]) -> list[Organism]:
    """Organism entries for catalog rows not already in the curated registry (matched by KEGG code)."""
    known = {o.kegg for o in ORGANISMS.values() if o.kegg}
    out = []
    for r in rows:
        if r["kegg"] in known or not r["taxid"]:
            continue
        out.append(Organism(
            key=r["kegg"], name=re.sub(r"\s*\(.*\)\s*$", "", r["name"]).strip() or r["name"], taxid=r["taxid"],
            gene_info=f"{r['ncbi_group']}/All_{r['ncbi_group']}", kegg=r["kegg"],
            string_taxid=r["taxid"] if r["string"] else None, aliases=(r["kegg"], r["taxid"]),
        ))
    return out


def register(rows: list[dict]) -> int:
    """Add catalog organisms to the in-process registry (organisms.get then resolves them)."""
    n = 0
    for o in organisms_from_catalog(rows):
        if o.key not in ORGANISMS:
            ORGANISMS[o.key] = o
            n += 1
    return n


def write_registry(rows: list[dict], path: Path) -> None:
    """Organism names/taxids for the web service (one JSON next to the releases)."""
    data = [{**asdict(o), "lineage": r["lineage"]} for o, r in
            ((o, next(x for x in rows if x["kegg"] == o.kegg)) for o in organisms_from_catalog(rows))]
    path.write_text(json.dumps(data))


# --------------------------------------------------------------------------- #
# Shared NCBI gene_info: one pass per group file, split by taxid
# --------------------------------------------------------------------------- #
def gene_info_for_taxid(shared: Path, group: str, taxid: str) -> Path | None:
    """Rows of All_<group>.gene_info.gz for one taxid (split once per group; None if NCBI has none)."""
    src = fetch(f"https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/{group}/All_{group}.gene_info.gz",
                shared / "ncbi" / f"All_{group}.gene_info.gz")
    split = shared / "ncbi" / f"split_{group}"
    done = split / ".complete"
    if not done.exists() or done.stat().st_mtime < src.stat().st_mtime:
        for old in split.glob("*.gene_info.gz"):
            old.unlink()
        split.mkdir(parents=True, exist_ok=True)
        with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh:
            header = fh.readline()
            cur, out = None, None
            for line in fh:
                tax = line.split("\t", 1)[0]
                if tax != cur:
                    if out:
                        out.close()
                    target = split / f"{tax}.gene_info.gz"
                    fresh = not target.exists()
                    out = gzip.open(target, "at", compresslevel=3)
                    if fresh:
                        out.write(header)
                    cur = tax
                out.write(line)
            if out:
                out.close()
        done.touch()
    f = split / f"{taxid}.gene_info.gz"
    return f if f.exists() else None
