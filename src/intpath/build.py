"""End-to-end release builder for one organism.

    raw/<organism>/...      downloaded source files (cached, not versioned)
    release/<organism>/     intpath_* tables, GMT, PPI network, stats.json

Pipeline: download -> ID normalisation -> pathway name matching + full
unification -> GO sets (propagated, losslessly merged) -> GO<->pathway links
-> PPI integration + pathway overlay -> release files + statistics.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from . import go as golib
from . import ppi as ppilib
from . import msigdb, sources
from .db import write_db
from .io import read_legacy_source, write_release
from .mapping import GeneMapper
from .model import GeneSet, SourcePathway
from .organisms import Organism
from .unify import merge_pathways

DEFAULT_PATHWAY_SOURCES = ("KEGG", "Reactome", "WikiPathways")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def source_stats(pws: list[SourcePathway]) -> dict:
    genes = {g for p in pws for g in p.genes}
    pairs = {pr for p in pws for pr in p.pairs}
    return {
        "pathways": len(pws),
        "genes": len(genes),
        "gene_pairs": len(pairs),
        "mean_genes_per_pathway": round(sum(len(p.genes) for p in pws) / max(len(pws), 1), 2),
        "mean_pairs_per_pathway": round(sum(len(p.pairs) for p in pws) / max(len(pws), 1), 2),
    }


def set_stats(sets: list[GeneSet]) -> dict:
    merged = [s for s in sets if len(s.members) > 1]
    return {
        "sets": len(sets),
        "merged_sets": len(merged),
        "member_pathways_in_merged_sets": sum(len(s.members) for s in merged),
        "genes": len({g for s in sets for g in s.genes}),
        "gene_pairs": len({pr for s in sets for pr in s.pairs}),
        "mean_genes_per_set": round(sum(s.size for s in sets) / max(len(sets), 1), 2),
        "mean_pairs_per_set": round(sum(len(s.pairs) for s in sets) / max(len(sets), 1), 2),
    }


def build(
    org: Organism,
    raw_root: str | Path,
    out_root: str | Path,
    *,
    pathway_sources: tuple[str, ...] = DEFAULT_PATHWAY_SOURCES,
    biocyc_col: str | None = None,
    with_go: bool = True,
    with_ppi: bool = True,
    string_min_score: int = 700,
    extra_ppi: tuple[tuple[str, str, str], ...] = (),  # (label, format[biogrid|mitab|tsv], path)
    topology: bool = True,
    public: bool = False,
    with_msigdb: bool = True,
) -> dict:
    raw, out = Path(raw_root) / org.key, Path(out_root) / org.key
    out.mkdir(parents=True, exist_ok=True)
    stats: dict = {"organism": org.name, "taxid": org.taxid, "built": time.strftime("%Y-%m-%d"), "sources": {}}

    log(f"{org.name}: downloading identifier and ontology files")
    files = sources.fetch_common(raw, org)
    mapper = GeneMapper.for_organism(org, raw)
    log(f"gene mapper: {len(mapper.exact)} exact keys, {len(mapper.alias)} aliases")

    extractors = {"KEGG": sources.kegg, "Reactome": sources.reactome, "WikiPathways": sources.wikipathways}
    pathways: list[SourcePathway] = []
    for src in pathway_sources:
        if public and src in sources.RESTRICTED_SOURCES:
            log(f"skip {src} (restricted licence, public build)")
            continue
        log(f"extracting {src}")
        fn = extractors[src]
        try:
            pws = fn(raw, org, mapper, topology=topology)
        except Exception as exc:  # a source without data for this organism must not stop the build
            log(f"  {src} unavailable for {org.name}: {exc}")
            stats["sources"][src] = {"error": str(exc)[:200]}
            continue
        pws = [p for p in pws if p.genes]
        stats["sources"][src] = source_stats(pws)
        log(f"  {src}: {stats['sources'][src]}")
        pathways += pws
    if biocyc_col and not public:
        pws = sources.biocyc_pathways_col(biocyc_col, mapper)
        stats["sources"]["BioCyc"] = source_stats(pws)
        pathways += pws

    log("merging related pathways (IntPath name alignment + full unification)")
    sets, matches = merge_pathways(pathways, legacy=False, organism=org.key)
    stats["pathways"] = set_stats(sets)
    stats["related_pathway_pairs"] = len(matches)
    with open(out / "related_pathways.tsv", "w") as fh:
        fh.write("source_a\tpathway_a\tsource_b\tpathway_b\talign_score\talign_ratio\tjaccard\n")
        for m in matches:
            ov = "" if m.overlap is None else f"{m.overlap:.3f}"
            fh.write(f"{m.a[0]}\t{m.a[1]}\t{m.b[0]}\t{m.b[1]}\t{m.score}\t{m.ratio:.4f}\t{ov}\n")
    log(f"  {stats['pathways']}")

    if with_go and "go_gaf" in files:
        log("building GO gene sets")
        terms = golib.parse_obo(files["go_obo"])
        ann = golib.parse_gaf(files["go_gaf"], terms, mapper)
        full = golib.propagate(ann, terms)
        go_sets = golib.build_go_sets(terms, full)
        stats["go"] = {
            "annotated_terms_propagated": len(full),
            "go_sets": len(go_sets),
            "go_terms_in_sets": sum(len(s.members) for s in go_sets),
            "by_namespace": {ns: sum(1 for s in go_sets if s.collection == ns) for ns in ("GO:BP", "GO:MF", "GO:CC")},
        }
        stats["go"]["pathway_links"] = golib.link_sets(go_sets, sets)
        log(f"  {stats['go']}")
        sets += go_sets

    msig = None
    if with_msigdb:
        log("MSigDB companion library")
        msig = msigdb.load(raw, org.key, mapper, sets, public=public)
        if msig is not None:
            stats["msigdb"] = {**msig.stats, "sets": len(msig.sets), "equivalents": len(msig.equivalents)}
            log(f"  MSigDB {msig.version}: {len(msig.sets)} sets, {len(msig.equivalents)} equivalence links")
            sets += msig.sets

    net = None
    if with_ppi:
        log("integrating PPIs (STRING, BioGRID, IntAct/MINT, HuRI)")
        shared = Path(raw_root) / "shared"
        feeds = []
        if org.string_taxid:
            links, info = sources.fetch_string(raw, org)
            feeds.append(("STRING", ppilib.parse_string(links, info, string_min_score)))
        for name, feed in sources.PPI_FEEDS.items():
            if name == "HuRI" and org.key != "sapiens":
                continue
            try:
                rows = feed(shared, org)
            except Exception as exc:
                log(f"  {name} unavailable: {exc}")
                continue
            feeds.append((name, rows))
        for label, fmt, path in extra_ppi:
            parser = {"biogrid": ppilib.parse_biogrid_tab3, "mitab": ppilib.parse_mitab, "tsv": ppilib.parse_pairs_tsv}[fmt]
            feeds.append((label, parser(path) if fmt == "tsv" else parser(path, taxid=org.taxid)))
        net, rep = ppilib.integrate(feeds, mapper)
        net.overlay_pathways([s for s in sets if s.collection == "pathway"])
        net.write(out / "intpath_ppi.tsv")
        stats["ppi"] = {"per_source": rep, **net.summary()}
        log(f"  {stats['ppi']}")

    stats["tier"] = "open" if public else "full"
    stats["versions"] = sources.source_versions(raw_root, org, msig.version if msig else None)
    stats["files"] = write_release(sets, out)
    release_genes = {g for s in sets for g in s.genes}
    if net is not None:
        release_genes |= {g for e in net.edges for g in e}
    aliases = mapper.alias_table(release_genes)
    stats["aliases"] = len(aliases)
    stats["mapping"] = dict(mapper.stats)
    log("writing release database")
    meta = {"organism_key": org.key, "organism": org.name, "taxid": org.taxid, "built": stats["built"],
            "tier": stats["tier"], "versions": stats["versions"], "stats": {k: v for k, v in stats.items() if k != "files"}}
    write_db(out / "intpath.sqlite", sets, meta=meta, aliases=aliases, gene_ids=mapper.entrez, ppi=net,
             matches=matches, msigdb_equivalents=msig.equivalents if msig else None)
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    log(f"release written to {out}")
    return stats


def rebuild_legacy(data_root: str | Path, out_root: str | Path, organism: str = "sapiens") -> dict:
    """Re-run the merge on the archived old IntPath normalized files (reproducibility check)."""
    n = Path(data_root) / organism / "normalized"
    order = [("K", "KEGG"), ("C", "BioCyc"), ("W", "WikiPathways")]  # legacy comparison order
    pws = []
    for code, folder in order:
        pws += read_legacy_source(
            n / folder / f"{organism}{folder}NormPthGEN", n / folder / f"{organism}{folder}NormPthGPR", code
        )
    sets, matches = merge_pathways(pws, legacy=True, organism=organism)
    out = Path(out_root) / f"legacy_{organism}"
    files = write_release(sets, out)
    stats = {"related_pathway_pairs": len(matches), **set_stats(sets), "files": files}
    (out / "stats.json").write_text(json.dumps(stats, indent=2))
    return stats
