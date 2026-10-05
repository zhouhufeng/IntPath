"""Protein-protein interaction integration (IntPathV2).

Each source is parsed to (symbolA, symbolB, evidence) after HGNC normalisation;
edges are undirected (sorted symbol pair) and fully unified across sources:
an edge keeps the set of supporting sources, the number of distinct
publications, detection methods and the best STRING score. Pathway gene pairs
from the integrated pathways can then be overlaid so every PPI edge knows the
pathways in which it is a curated relation ("pathway-supported" edges).

Supported formats
  * STRING   protein.(physical.)links(.detailed).vX.txt.gz + protein.info (ENSP -> symbol)
  * BioGRID  BIOGRID-ORGANISM-Homo_sapiens-*.tab3.txt (physical only by default)
  * PSI-MITAB 2.5-2.7 (IntAct, MINT, HuRI, ...)
  * any 2-3 column TSV (geneA, geneB[, score]) - e.g. the legacy sapiensSTRING
"""

from __future__ import annotations

import io
import re
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator

from .io import open_text
from .mapping import GeneMapper

Edge = tuple[str, str]


def canon(a: str, b: str) -> Edge:
    return (a, b) if a <= b else (b, a)


EXPERIMENTAL_SOURCES = frozenset({"BioGRID", "IntAct", "MINT", "HuRI"})


@dataclass(slots=True)
class EdgeEvidence:
    sources: set[str] = field(default_factory=set)
    pmids: set[str] = field(default_factory=set)
    methods: set[str] = field(default_factory=set)
    string_score: int = 0
    pathways: set[str] = field(default_factory=set)  # set ids where the pair is a curated pathway relation


class PPINetwork:
    def __init__(self) -> None:
        self.edges: dict[Edge, EdgeEvidence] = {}

    def add(self, a: str, b: str, source: str, *, pmid: str = "", method: str = "", score: int = 0) -> None:
        if not a or not b or a == b:
            return
        ev = self.edges.setdefault(canon(a, b), EdgeEvidence())
        ev.sources.add(source)
        if pmid:
            ev.pmids.add(pmid)
        if method:
            ev.methods.add(method)
        if score > ev.string_score:
            ev.string_score = score

    @staticmethod
    def tier(ev: EdgeEvidence) -> str:
        """Confidence tier (stored, never used to drop edges; users filter at analysis time).

        high:   >= 2 experimental sources, or >= 2 publications, or STRING >= 900
        medium: one experimental source, or STRING >= 700
        low:    everything else
        """
        n_exp = len(ev.sources & EXPERIMENTAL_SOURCES)
        if n_exp >= 2 or len(ev.pmids) >= 2 or ev.string_score >= 900:
            return "high"
        if n_exp >= 1 or ev.string_score >= 700:
            return "medium"
        return "low"

    def neighbours(self) -> dict[str, set[str]]:
        adj: dict[str, set[str]] = {}
        for a, b in self.edges:
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
        return adj

    def filter(self, keep: Callable[[EdgeEvidence], bool]) -> "PPINetwork":
        net = PPINetwork()
        net.edges = {e: ev for e, ev in self.edges.items() if keep(ev)}
        return net

    def overlay_pathways(self, sets) -> int:
        """Mark edges that are also curated gene pairs in integrated pathways."""
        n = 0
        for gs in sets:
            for a, b in gs.pairs:
                ev = self.edges.get(canon(a, b))
                if ev is not None:
                    ev.pathways.add(gs.id)
                    n += 1
        return n

    def write(self, path: str | Path) -> None:
        with open(path, "w") as fh:
            fh.write("gene_a\tgene_b\tsources\tn_sources\tn_pmids\tmethods\tstring_score\ttier\tpathway_sets\n")
            for (a, b), ev in sorted(self.edges.items()):
                fh.write(
                    f"{a}\t{b}\t{','.join(sorted(ev.sources))}\t{len(ev.sources)}\t{len(ev.pmids)}\t"
                    f"{'|'.join(sorted(ev.methods))}\t{ev.string_score}\t{self.tier(ev)}\t{','.join(sorted(ev.pathways))}\n"
                )

    @classmethod
    def read(cls, path: str | Path) -> "PPINetwork":
        net = cls()
        with open_text(path) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                ev = EdgeEvidence(
                    sources=set(f[2].split(",")),
                    methods={m for m in f[5].split("|") if m},
                    string_score=int(f[6] or 0),
                    pathways={p for p in f[8].split(",") if p},
                )
                ev.pmids = {f"n{i}" for i in range(int(f[4] or 0))}  # counts only
                net.edges[(f[0], f[1])] = ev
        return net

    def summary(self) -> dict:
        by_src = Counter(s for ev in self.edges.values() for s in ev.sources)
        by_n = Counter(len(ev.sources) for ev in self.edges.values())
        genes = {g for e in self.edges for g in e}
        tiers = Counter(self.tier(ev) for ev in self.edges.values())
        return {
            "edges": len(self.edges),
            "edges_by_tier": dict(tiers),
            "genes": len(genes),
            "edges_per_source": dict(by_src),
            "edges_by_n_sources": {str(k): v for k, v in sorted(by_n.items())},
            "pathway_supported_edges": sum(1 for ev in self.edges.values() if ev.pathways),
        }


# --------------------------------------------------------------------------- #
# Parsers: each yields (idA, idB, pmid, method, score) in source identifiers
# --------------------------------------------------------------------------- #
def parse_string(links: str | Path, info: str | Path, min_score: int = 700) -> Iterator[tuple[str, str, str, str, int]]:
    ensp2sym: dict[str, str] = {}
    with open_text(info) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            ensp2sym[f[0]] = f[1]
    with open_text(links) as fh:
        header = fh.readline().split()
        si = header.index("combined_score")
        for line in fh:
            f = line.split()
            s = int(f[si])
            if s >= min_score and f[0] in ensp2sym and f[1] in ensp2sym:
                yield ensp2sym[f[0]], ensp2sym[f[1]], "", "", s


def parse_biogrid_tab3(
    path: str | Path, physical_only: bool = True, taxid: str = "9606", organism: str = ""
) -> Iterator[tuple]:
    """BioGRID tab3: one organism file, or the BIOGRID-ORGANISM-LATEST zip with ``organism``
    (binomial, e.g. "Homo sapiens") selecting the member BIOGRID-ORGANISM-Homo_sapiens-*.tab3.txt."""
    if str(path).endswith(".zip"):
        stem = "-" + "_".join(organism.split()[:2])  # e.g. -Saccharomyces_cerevisiae(_S288c)-
        with zipfile.ZipFile(path) as zf:
            members = [m for m in zf.namelist() if stem in m]
            if not members:
                raise FileNotFoundError(f"no BioGRID member for {organism!r} in {path}")
            with io.TextIOWrapper(zf.open(members[0]), encoding="utf-8", errors="replace") as fh:
                yield from _biogrid_lines(fh, physical_only, taxid)
        return
    with open_text(path) as fh:
        yield from _biogrid_lines(fh, physical_only, taxid)


def _biogrid_lines(fh, physical_only: bool, taxid: str) -> Iterator[tuple]:
    header = fh.readline().rstrip("\n").split("\t")
    c = {h: i for i, h in enumerate(header)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[c["Organism ID Interactor A"]] != taxid or f[c["Organism ID Interactor B"]] != taxid:
            continue
        if physical_only and f[c["Experimental System Type"]] != "physical":
            continue
        yield (
            f[c["Entrez Gene Interactor A"]],
            f[c["Entrez Gene Interactor B"]],
            f[c["Publication Source"]],
            f[c["Experimental System"]],
            0,
        )


_MITAB_ID = re.compile(r"(uniprotkb|entrez gene/locuslink|ensembl):([^|()\s]+)", re.I)
_MITAB_TAX = re.compile(r"taxid:(-?\d+)")


_MITAB_SOURCE_DB = re.compile(r'\(([^)]*)\)')


def parse_mitab(path: str | Path, taxid: str = "9606", split_source: bool = True) -> Iterator[tuple]:
    """PSI-MITAB 2.5+: columns 1-2 ids, 3-4 alt ids, 7 method, 9 pmids, 10-11 taxids, 13 source db.

    With ``split_source`` each row also yields its curating database, so MINT
    evidence (curated into IntAct under the IMEx agreement) is kept as its own
    source "MINT"; every other IMEx curator is reported as "IntAct".
    """
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            ta, tb = _MITAB_TAX.findall(f[9])[:1], _MITAB_TAX.findall(f[10])[:1]
            if ta != [taxid] or tb != [taxid]:
                continue
            ida, idb = _MITAB_ID.findall(f[0] + "|" + f[2]), _MITAB_ID.findall(f[1] + "|" + f[3])
            if not ida or not idb:
                continue
            pm = re.findall(r"pubmed:(\d+)", f[8])
            meth = re.findall(r'"MI:\d+"\(([^)]*)\)', f[6])
            row = (ida[0][1], idb[0][1], pm[0] if pm else "", meth[0] if meth else "", 0)
            if split_source:
                db = " ".join(_MITAB_SOURCE_DB.findall(f[12])).lower() if len(f) > 12 else ""
                row += ("MINT" if "mint" in db else "IntAct",)
            yield row


def parse_pairs_tsv(path: str | Path) -> Iterator[tuple]:
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not line.startswith("#"):
                yield f[0], f[1], "", "", int(f[2]) if len(f) > 2 and f[2].isdigit() else 0


def integrate(
    sources: Iterable[tuple[str, Iterable[tuple]]], mapper: GeneMapper | None = None
) -> tuple[PPINetwork, dict[str, dict]]:
    """Map each source to HGNC symbols and merge into one network."""
    net = PPINetwork()
    report: dict[str, dict] = {}
    for name, rows in sources:
        n_in = n_kept = 0
        cache: dict[str, str | None] = {}
        for row in rows:
            a, b, pmid, method, score = row[:5]
            src = row[5] if len(row) > 5 else name
            n_in += 1
            if mapper is not None:
                ma = cache[a] if a in cache else cache.setdefault(a, mapper.map(a))
                mb = cache[b] if b in cache else cache.setdefault(b, mapper.map(b))
            else:
                ma, mb = a, b
            if ma and mb and ma != mb:
                net.add(ma, mb, src, pmid=pmid, method=method, score=score)
                n_kept += 1
        report[name] = {"rows": n_in, "mapped_rows": n_kept, "unmapped_ids": sum(1 for v in cache.values() if v is None)}
    return net, report
