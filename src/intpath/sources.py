"""Source download + extraction for any registered organism.

Every extractor returns :class:`SourcePathway` objects with genes (and gene
pairs where the source has topology) already normalised by a GeneMapper.

Licences (check before redistributing any derived file):
  WikiPathways CC0 | Reactome CC0 (CC-BY 4.0 for some files) | GO CC-BY 4.0 |
  NCBI Gene public domain | HGNC CC0 | STRING CC-BY 4.0 | BioGRID MIT |
  IntAct CC-BY 4.0 | KEGG: academic use via REST, redistribution of KEGG-derived
  data requires a licence from Pathway Solutions | BioCyc: subscription licence.
The release builder therefore tags every gene set with its sources, and the
public export can exclude restricted sources (``--public``).
"""

from __future__ import annotations

import re
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from itertools import combinations
from pathlib import Path
from typing import Iterator

from .io import open_text
from .mapping import GeneMapper
from .model import SourcePathway
from .organisms import Organism

RESTRICTED_SOURCES = {"KEGG", "BioCyc"}
UA = {"User-Agent": "IntPath/3 (https://intpath.genohub.org)"}

URLS = {
    "hgnc": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "go_obo": "https://current.geneontology.org/ontology/go-basic.obo",
    "reactome_ncbi": "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt",
    "reactome_lowest": "https://reactome.org/download/current/NCBI2Reactome.txt",
    "reactome_relations": "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
    "wikipathways_gmt_index": "https://data.wikipathways.org/current/gmt/",
    "wikipathways_gpml_index": "https://data.wikipathways.org/current/gpml/",
    "kegg_rest": "https://rest.kegg.jp",
    "string": "https://stringdb-downloads.org/download/{kind}.v12.0/{taxid}.{kind}.v12.0.txt.gz",
    "biogrid": "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip",
}


# --------------------------------------------------------------------------- #
# Download helpers
# --------------------------------------------------------------------------- #
def fetch(url: str, dest: str | Path, *, force: bool = False, retries: int = 3) -> Path:
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r, open(tmp, "wb") as fh:
                while chunk := r.read(1 << 20):
                    fh.write(chunk)
            tmp.rename(dest)
            return dest
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    return dest


def fetch_text(url: str) -> str:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120) as r:
        return r.read().decode("utf-8", "replace")


def latest_listing(index_url: str, pattern: str) -> str:
    names = sorted(set(re.findall(pattern, fetch_text(index_url))))
    if not names:
        raise FileNotFoundError(f"no file matching {pattern!r} at {index_url}")
    return index_url + names[-1]


def fetch_common(raw: Path, org: Organism) -> dict[str, Path]:
    """Files every build needs: NCBI gene_info (+HGNC for human) and GO."""
    out = {"gene_info": fetch(org.gene_info_url, raw / "ncbi" / f"{org.key}.gene_info.gz")}
    if org.hgnc:
        out["hgnc"] = fetch(URLS["hgnc"], raw / "hgnc" / "hgnc_complete_set.txt")
    out["go_obo"] = fetch(URLS["go_obo"], raw / "go" / "go-basic.obo")
    if org.go_gaf_url:
        out["go_gaf"] = fetch(org.go_gaf_url, raw / "go" / f"{org.go_gaf}.gaf.gz")
    return out


# --------------------------------------------------------------------------- #
# Reactome (all organisms with Reactome species; NCBI Gene ids)
# --------------------------------------------------------------------------- #
def reactome(raw: Path, org: Organism, mapper: GeneMapper, lowest_level: bool = False) -> list[SourcePathway]:
    if not org.reactome:
        return []
    key = "reactome_lowest" if lowest_level else "reactome_ncbi"
    f = fetch(URLS[key], raw / "reactome" / Path(URLS[key]).name)
    pws: dict[str, SourcePathway] = {}
    with open_text(f) as fh:
        for line in fh:
            gid, pid, _url, name, ev, species = line.rstrip("\n").split("\t")[:6]
            if species != org.reactome:
                continue
            sym = mapper.map(gid)
            if sym:
                pws.setdefault(pid, SourcePathway("Reactome", name, pid)).genes.add(sym)
    return list(pws.values())


# --------------------------------------------------------------------------- #
# WikiPathways (GMT for membership, GPML for topology)
# --------------------------------------------------------------------------- #
def wikipathways(raw: Path, org: Organism, mapper: GeneMapper, topology: bool = True) -> list[SourcePathway]:
    if not org.wikipathways:
        return []
    url = latest_listing(URLS["wikipathways_gmt_index"], rf"wikipathways-\d+-gmt-{org.wikipathways}\.gmt")
    f = fetch(url, raw / "wikipathways" / Path(url).name)
    pws: dict[str, SourcePathway] = {}
    with open_text(f) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            meta = parts[0].split("%")  # name%WikiPathways_YYYYMMDD%WPnnn%Species
            name, wpid = meta[0], meta[2] if len(meta) > 2 else parts[0]
            p = pws.setdefault(wpid, SourcePathway("WikiPathways", name, wpid))
            p.genes.update(s for s in (mapper.map(g) for g in parts[2:] if g) if s)
    if topology:
        try:
            gurl = latest_listing(URLS["wikipathways_gpml_index"], rf"wikipathways-\d+-gpml-{org.wikipathways}\.zip")
            gz = fetch(gurl, raw / "wikipathways" / Path(gurl).name)
            for wpid, pairs in gpml_pairs(gz, mapper):
                if wpid in pws:
                    for pr, rel in pairs.items():
                        pws[wpid].pairs.setdefault(pr, set()).update(rel)
        except Exception as exc:  # topology is a bonus; membership still usable
            print(f"[wikipathways] GPML topology skipped: {exc}")
    return list(pws.values())


_ARROW_REL = {
    "mim-conversion": "ECrel",
    "mim-transcription-translation": "GErel",
    "TranscriptionTranslation": "GErel",
}


def gpml_pairs(zip_path: Path, mapper: GeneMapper) -> Iterator[tuple[str, dict]]:
    """Port of Wiki.java: Interaction lines -> PPrel/ECrel/GErel ("graphRel"); Groups -> GPrel ("groupRel")."""
    with zipfile.ZipFile(zip_path) as zf:
        for member in zf.namelist():
            if not member.endswith(".gpml"):
                continue
            m = re.search(r"(WP\d+)", member)
            if not m:
                continue
            try:
                root = ET.fromstring(zf.read(member))
            except ET.ParseError:
                continue
            ns = root.tag.split("}")[0] + "}" if root.tag.startswith("{") else ""
            node_genes: dict[str, set[str]] = {}
            group_of: dict[str, str] = {}
            for dn in root.iter(f"{ns}DataNode"):
                if dn.get("Type") not in ("GeneProduct", "Protein", "Rna"):
                    continue
                x = dn.find(f"{ns}Xref")
                cands = [x.get("ID", "") if x is not None else "", dn.get("TextLabel", "")]
                sym = next((s for s in (mapper.map(c) for c in cands if c) if s), None)
                if sym:
                    gid = dn.get("GraphId") or dn.get("elementId") or ""
                    node_genes[gid] = {sym}
                    if dn.get("GroupRef"):
                        group_of[gid] = dn.get("GroupRef")
            groups: dict[str, set[str]] = {}
            for gid, gref in group_of.items():
                groups.setdefault(gref, set()).update(node_genes[gid])
            for grp in root.iter(f"{ns}Group"):  # groups addressable by GraphId in interactions
                gref, ggid = grp.get("GroupId"), grp.get("GraphId")
                if ggid and gref in groups:
                    node_genes[ggid] = groups[gref]
            pairs: dict[tuple[str, str], set[str]] = {}
            for genes in groups.values():
                for a, b in combinations(sorted(genes), 2):
                    pairs.setdefault((a, b), set()).add("GPrel")
            for it in root.iter(f"{ns}Interaction"):
                pts = it.findall(f".//{ns}Point")
                if len(pts) < 2:
                    continue
                src, dst = pts[0].get("GraphRef"), pts[-1].get("GraphRef")
                arrow = pts[-1].get("ArrowHead", "")
                rel = _ARROW_REL.get(arrow, "PPrel")
                for a in node_genes.get(src, ()):
                    for b in node_genes.get(dst, ()):
                        if a != b:
                            pairs.setdefault((a, b), set()).add(rel)
            yield m.group(1), pairs


# --------------------------------------------------------------------------- #
# KEGG (REST + KGML; academic use - see licence note above)
# --------------------------------------------------------------------------- #
def kegg(raw: Path, org: Organism, mapper: GeneMapper, topology: bool = True, delay: float = 0.35) -> list[SourcePathway]:
    if not org.kegg:
        return []
    base = URLS["kegg_rest"]
    d = raw / "kegg" / org.kegg
    lst = fetch(f"{base}/list/pathway/{org.kegg}", d / "list_pathway.txt")
    lnk = fetch(f"{base}/link/{org.kegg}/pathway", d / "link_pathway.txt")
    pws: dict[str, SourcePathway] = {}
    for line in open(lst):
        pid, name = line.rstrip("\n").split("\t")
        pid = pid.replace("path:", "")
        pws[pid] = SourcePathway("KEGG", name, pid)
    kid2sym: dict[str, str | None] = {}

    def sym_of(kid: str) -> str | None:
        if kid not in kid2sym:
            kid2sym[kid] = mapper.map(kid.split(":", 1)[1])
        return kid2sym[kid]

    for line in open(lnk):
        pid, kid = line.rstrip("\n").split("\t")
        pid = pid.replace("path:", "")
        s = sym_of(kid)
        if s and pid in pws:
            pws[pid].genes.add(s)
    if topology:
        for pid, p in pws.items():
            kgml = d / "kgml" / f"{pid}.xml"
            if not kgml.exists():
                try:
                    fetch(f"{base}/get/{pid}/kgml", kgml)
                    time.sleep(delay)  # be polite to the KEGG REST server
                except Exception:
                    continue
            p.pairs = kgml_pairs(kgml, sym_of)
    return list(pws.values())


_KEGG_REL = {"PPrel": "PPrel", "PCrel": "PPrel", "ECrel": "ECrel", "GErel": "GErel"}  # maplink dropped (Normalize.mpRel)


def kgml_pairs(path: Path, sym_of) -> dict[tuple[str, str], set[str]]:
    """Port of KEGG.java: relations between gene entries, and groups -> GPrel."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return {}
    entry_genes: dict[str, set[str]] = {}
    groups: list[str] = []
    for e in root.iter("entry"):
        if e.get("type") == "gene":
            entry_genes[e.get("id")] = {s for s in (sym_of(k) for k in e.get("name", "").split()) if s}
        elif e.get("type") == "group":
            groups.append(e.get("id"))
            entry_genes[e.get("id")] = set()
    for e in root.iter("entry"):
        if e.get("type") == "group":
            for c in e.iter("component"):
                entry_genes[e.get("id")] |= entry_genes.get(c.get("id"), set())
    pairs: dict[tuple[str, str], set[str]] = {}
    for gid in groups:
        for a, b in combinations(sorted(entry_genes[gid]), 2):
            pairs.setdefault((a, b), set()).add("GPrel")
    for r in root.iter("relation"):
        rel = _KEGG_REL.get(r.get("type", ""))
        if not rel:
            continue
        for a in entry_genes.get(r.get("entry1"), ()):
            for b in entry_genes.get(r.get("entry2"), ()):
                if a != b:
                    pairs.setdefault((a, b), set()).add(rel)
    return pairs


# --------------------------------------------------------------------------- #
# BioCyc (licensed): gene membership from a local pathways.col
# --------------------------------------------------------------------------- #
def biocyc_pathways_col(path: str | Path, mapper: GeneMapper) -> list[SourcePathway]:
    pws = []
    with open_text(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = f
                continue
            rec = dict(zip(header, f))
            p = SourcePathway("BioCyc", re.sub(r"<[^>]+>", "", rec.get("NAME", "")), rec.get("UNIQUE-ID", ""))
            for h, v in rec.items():
                if h.startswith("GENE-NAME") and v:
                    s = mapper.map(v)
                    if s:
                        p.genes.add(s)
            if p.genes:
                pws.append(p)
    return pws


# --------------------------------------------------------------------------- #
# PPI sources
# --------------------------------------------------------------------------- #
def fetch_string(raw: Path, org: Organism, kind: str = "protein.physical.links") -> tuple[Path, Path]:
    tx = org.string_taxid or org.taxid
    links = fetch(URLS["string"].format(kind=kind, taxid=tx), raw / "string" / f"{tx}.{kind}.v12.0.txt.gz")
    info = fetch(URLS["string"].format(kind="protein.info", taxid=tx), raw / "string" / f"{tx}.protein.info.v12.0.txt.gz")
    return links, info
