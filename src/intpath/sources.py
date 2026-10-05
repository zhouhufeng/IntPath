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
    "intact": "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip",
    # the certificate covers interactome-atlas.org, not www.interactome-atlas.org
    "huri": "https://interactome-atlas.org/data/HuRI.tsv",
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
    gi = raw / "ncbi" / f"{org.key}.gene_info.gz"
    if org.gene_info_is_shared:  # e.g. All_Archaea_Bacteria: download once, keep only this taxid
        gi = raw.parent / "shared" / "ncbi" / org.gene_info_url.rsplit("/", 1)[1]
    out = {"gene_info": fetch(org.gene_info_url, gi)}
    out["uniprot"] = fetch_uniprot(raw, org)
    if org.hgnc:
        out["hgnc"] = fetch(URLS["hgnc"], raw / "hgnc" / "hgnc_complete_set.txt")
    out["go_obo"] = fetch(URLS["go_obo"], raw / "go" / "go-basic.obo")
    if org.go_gaf_url:
        out["go_gaf"] = fetch(org.go_gaf_url, raw / "go" / f"{org.go_gaf}.gaf.gz")
    else:  # no GO Consortium GAF (e.g. M. tuberculosis): UniProt GO annotations as a minimal GAF
        out["go_gaf"] = uniprot_gaf(raw, org)
    return out


UNIPROT_BY_ORG = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism"
UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=organism_id:{taxid}&fields={fields}"


def fetch_uniprot(raw: Path, org: Organism) -> Path:
    """UniProt accession -> GeneID mapping: the by_organism idmapping file, else the REST stream."""
    if org.uniprot:
        return fetch(f"{UNIPROT_BY_ORG}/{org.uniprot}_idmapping_selected.tab.gz",
                     raw / "uniprot" / f"{org.key}_idmapping_selected.tab.gz")
    return fetch(UNIPROT_REST.format(taxid=org.taxid, fields="accession,id,xref_geneid,gene_oln,go_id"),
                 raw / "uniprot" / f"{org.key}_uniprot_rest.tsv")


def uniprot_gaf(raw: Path, org: Organism) -> Path:
    """Write a minimal GAF 2.2 from UniProt's GO cross-references (evidence code unknown -> 'UniProt')."""
    src = fetch(UNIPROT_REST.format(taxid=org.taxid, fields="accession,id,xref_geneid,gene_oln,go_id"),
                raw / "uniprot" / f"{org.key}_uniprot_rest.tsv")
    out = raw / "go" / f"{org.key}_uniprot.gaf"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(src) as fh, open(out, "w") as w:
        w.write("!gaf-version: 2.2\n! derived from UniProtKB GO cross-references by IntPathV2\n")
        header = fh.readline().rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            acc, go_field = f[col["Entry"]], f[col.get("Gene Ontology IDs", len(f) - 1)]
            for go in (g.strip() for g in go_field.split(";")):
                if go.startswith("GO:"):
                    w.write("\t".join(["UniProtKB", acc, acc, "", go, "", "UniProt", "", "", "", "", "protein",
                                        f"taxon:{org.taxid}", "", "UniProt", "", ""]) + "\n")
    return out


# --------------------------------------------------------------------------- #
# Reactome (all organisms with Reactome species; NCBI Gene ids)
# --------------------------------------------------------------------------- #
REACTOME = "https://reactome.org/download/current"
_REACTOME_REL = {"physical association": "PPrel", "association": "PPrel", "direct interaction": "PPrel"}


def reactome_version() -> str:
    try:
        return fetch_text("https://reactome.org/ContentService/data/database/version").strip()
    except Exception:
        return ""


def reactome(raw: Path, org: Organism, mapper: GeneMapper, topology: bool = True) -> list[SourcePathway]:
    """Reactome pathways (all hierarchy levels) with parent links and, where available, gene pairs.

    Gene pairs come from Reactome's interaction export. Each pair carries a
    context: a complex (-> GPrel, placed in the pathways of Complex_2_Pathway)
    or a reaction (-> PPrel, placed in the lowest-level pathways that contain
    both genes). Pairs are then propagated to every ancestor pathway, matching
    the all-levels membership.
    """
    if not org.reactome:
        return []
    d = raw / "reactome"
    f = fetch(f"{REACTOME}/NCBI2Reactome_All_Levels.txt", d / "NCBI2Reactome_All_Levels.txt")
    pws: dict[str, SourcePathway] = {}
    with open_text(f) as fh:
        for line in fh:
            gid, pid, _url, name, _ev, species = line.rstrip("\n").split("\t")[:6]
            if species != org.reactome:
                continue
            sym = mapper.map(gid)
            if sym:
                pws.setdefault(pid, SourcePathway("Reactome", name, pid)).genes.add(sym)

    rel = fetch(f"{REACTOME}/ReactomePathwaysRelation.txt", d / "ReactomePathwaysRelation.txt")
    children: dict[str, set[str]] = {}
    with open_text(rel) as fh:
        for line in fh:
            parent, child = line.rstrip("\n").split("\t")[:2]
            if child in pws:
                pws[child].parents.add(parent)
            children.setdefault(parent, set()).add(child)
    if not topology:
        return list(pws.values())

    slug = org.reactome.lower().replace(" ", "_")
    try:
        inter = fetch(f"{REACTOME}/interactors/reactome.{slug}.interactions.tab-delimited.txt",
                      d / f"reactome.{slug}.interactions.tab-delimited.txt")
    except Exception as exc:
        print(f"[reactome] no interaction export for {org.name}: {exc}")
        return list(pws.values())
    lowest: dict[str, set[str]] = {}  # gene -> lowest-level pathways
    low = fetch(f"{REACTOME}/NCBI2Reactome.txt", d / "NCBI2Reactome.txt")
    with open_text(low) as fh:
        for line in fh:
            gid, pid, _u, _n, _e, species = line.rstrip("\n").split("\t")[:6]
            if species == org.reactome:
                sym = mapper.map(gid)
                if sym:
                    lowest.setdefault(sym, set()).add(pid)
    complex_pw: dict[str, set[str]] = {}
    if org.key == "sapiens":
        c2p = fetch(f"{REACTOME}/Complex_2_Pathway_human.txt", d / "Complex_2_Pathway_human.txt")
        with open_text(c2p) as fh:
            next(fh)
            for line in fh:
                cx, pid = line.rstrip("\n").split("\t")[:2]
                complex_pw.setdefault(cx, set()).add(pid)

    ancestors_memo: dict[str, set[str]] = {}

    def with_ancestors(pid: str) -> set[str]:
        if pid not in ancestors_memo:
            out = {pid}
            for par in pws[pid].parents if pid in pws else ():
                out |= with_ancestors(par)
            ancestors_memo[pid] = out
        return ancestors_memo[pid]

    def gene_of(field: str) -> str | None:
        for tok in field.split("|"):
            tok = tok.split(":", 1)[-1]
            if tok and tok != "-":
                s = mapper.map(tok)
                if s:
                    return s
        return None

    with open_text(inter) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            a = gene_of(f[2]) or gene_of(f[0]) or gene_of(f[1])
            b = gene_of(f[5]) or gene_of(f[3]) or gene_of(f[4])
            if not a or not b or a == b:
                continue
            ctx = f[7].replace("reactome:", "")
            if ctx in complex_pw:
                targets, kind = complex_pw[ctx], "GPrel"
            else:
                targets = lowest.get(a, set()) & lowest.get(b, set())
                kind = _REACTOME_REL.get(f[6], "PPrel")
            pair = (a, b) if a <= b else (b, a)
            for pid in targets:
                for anc in with_ancestors(pid):
                    if anc in pws:
                        pws[anc].pairs.setdefault(pair, set()).add(kind)
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
def fetch_biogrid(shared: Path) -> Path:
    return fetch(URLS["biogrid"], shared / "biogrid" / "BIOGRID-ORGANISM-LATEST.tab3.zip")


def fetch_intact(shared: Path) -> Path:
    """All IntAct evidence (PSI-MITAB 2.7); rows curated by MINT are labelled MINT by parse_mitab."""
    return fetch(URLS["intact"], shared / "intact" / "intact.zip")


def intact_for_taxid(shared: Path, taxid: str) -> Path:
    """Same-species IntAct rows for one taxid, split once (for every registry organism) from intact.zip.

    intact.txt is ~10 GB; one pass writes small per-taxid files reused by every build.
    """
    from .organisms import ORGANISMS

    split = shared / "intact" / "split"
    target = split / f"{taxid}.mitab.txt"
    src = fetch_intact(shared)
    if target.exists() and target.stat().st_mtime >= src.stat().st_mtime:
        return target
    split.mkdir(parents=True, exist_ok=True)
    taxids = {o.taxid for o in ORGANISMS.values()} | {taxid}
    tax = re.compile(r"taxid:(-?\d+)")
    outs = {t: open(split / f"{t}.mitab.txt.part", "w") for t in taxids}
    try:
        with open_text(src) as fh:
            header = fh.readline()
            for out in outs.values():
                out.write(header)
            for line in fh:
                f = line.split("\t", 11)
                if len(f) < 11:
                    continue
                ta, tb = tax.search(f[9]), tax.search(f[10])
                if ta and tb and ta.group(1) == tb.group(1) and ta.group(1) in outs:
                    outs[ta.group(1)].write(line)
    finally:
        for out in outs.values():
            out.close()
    for t in taxids:
        (split / f"{t}.mitab.txt.part").replace(split / f"{t}.mitab.txt")
    return target


def fetch_huri(shared: Path) -> Path:
    """HuRI, the Human Reference Interactome (Luck et al. 2020; interactome-atlas.org)."""
    return fetch(URLS["huri"], shared / "huri" / "HuRI.tsv")


def _ppi_feeds():
    """name -> callable(shared_dir, organism) -> rows (idA, idB, pmid, method, score[, source])."""
    from . import ppi

    return {
        "BioGRID": lambda shared, org: ppi.parse_biogrid_tab3(fetch_biogrid(shared), taxid=org.taxid, organism=org.name),
        # IntAct rows; evidence curated by MINT is labelled "MINT"
        "IntAct": lambda shared, org: ppi.parse_mitab(intact_for_taxid(shared, org.taxid), taxid=org.taxid),
        "HuRI": lambda shared, org: ppi.parse_pairs_tsv(fetch_huri(shared)) if org.key == "sapiens" else iter(()),
    }


PPI_FEEDS = _ppi_feeds()


def source_versions(raw_root: str | Path, org: Organism, msigdb_version: str | None = None) -> dict[str, str]:
    """Best-effort record of every source release used in a build (stored in stats.json and the DB)."""
    raw_root = Path(raw_root)
    raw = raw_root / org.key
    v: dict[str, str] = {"STRING": "12.0"}
    rv = reactome_version()
    if rv:
        v["Reactome"] = rv
    gmts = sorted((raw / "wikipathways").glob("wikipathways-*-gmt-*.gmt"))
    if gmts:
        v["WikiPathways"] = re.search(r"wikipathways-(\d+)-", gmts[-1].name).group(1)
    obo = raw / "go" / "go-basic.obo"
    if obo.exists():
        with open(obo) as fh:
            for line in fh:
                if line.startswith("data-version:"):
                    v["GO"] = line.split(":", 1)[1].strip()
                    break
    bg = raw_root / "shared" / "biogrid" / "BIOGRID-ORGANISM-LATEST.tab3.zip"
    if bg.exists():
        with zipfile.ZipFile(bg) as zf:
            m = re.search(r"-(\d+\.\d+\.\d+)\.tab3", zf.namelist()[0])
            if m:
                v["BioGRID"] = m.group(1)
    for name, path in (("IntAct/MINT", raw_root / "shared" / "intact" / "intact.zip"),
                       ("HuRI", raw_root / "shared" / "huri" / "HuRI.tsv"),
                       ("NCBI Gene", raw / "ncbi" / f"{org.key}.gene_info.gz"),
                       ("KEGG", raw / "kegg" / (org.kegg or "-") / "list_pathway.txt")):
        if path.exists():
            v[name] = "downloaded " + time.strftime("%Y-%m-%d", time.localtime(path.stat().st_mtime))
    if msigdb_version:
        v["MSigDB"] = msigdb_version
    return v


def fetch_string(raw: Path, org: Organism, kind: str = "protein.physical.links") -> tuple[Path, Path]:
    tx = org.string_taxid or org.taxid
    links = fetch(URLS["string"].format(kind=kind, taxid=tx), raw / "string" / f"{tx}.{kind}.v12.0.txt.gz")
    info = fetch(URLS["string"].format(kind="protein.info", taxid=tx), raw / "string" / f"{tx}.protein.info.v12.0.txt.gz")
    return links, info
