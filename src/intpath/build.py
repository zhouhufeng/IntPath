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
    legacy_data_root: str | Path | None = "Data",
    exports: bool = True,
    merge_method: str = "intpath",
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
    elif not public and legacy_data_root:
        pws = sources.biocyc_legacy(legacy_data_root, org, mapper)
        if pws:
            log(f"BioCyc from the old IntPath archive: {len(pws)} pathways")
            stats["sources"]["BioCyc"] = source_stats(pws)
            pathways += pws

    log("merging related pathways (IntPath name alignment + full unification)")
    sets, matches = merge_pathways(pathways, legacy=False, organism=org.key, method=merge_method)
    stats["merge_method"] = merge_method
    if merge_method == "intpath":
        from .review import summary, write_log

        stats["merge_review"] = summary(matches)
        write_log(matches, out / "merge_review.tsv")
        log(f"  merge review: {stats['merge_review']}")
    stats["pathways"] = set_stats(sets)
    # the meaning of every member pathway travels with its merged IntPath pathway
    from . import meaning

    desc: dict[str, str] = {}
    for label, fn in (("Reactome", lambda: meaning.reactome(raw)), ("WikiPathways", lambda: meaning.wikipathways(raw)),
                      ("KEGG", lambda: meaning.kegg(Path(raw_root) / "shared",
                                                    [sid for s in sets for src, sid in s.members if src == "KEGG"]))):
        if any(src == label for s in sets for src, _ in s.members):
            try:
                desc.update(fn())
            except Exception as exc:
                log(f"  {label} descriptions unavailable: {exc}")
    for s in sets:
        for _src, sid in s.members:
            if sid in desc:
                s.member_desc[sid] = desc[sid]
    stats["pathways"]["with_description"] = sum(1 for s in sets if s.member_desc)
    stats["related_pathway_pairs"] = sum(1 for m in matches if m.decision == "accept")
    with open(out / "related_pathways.tsv", "w") as fh:
        fh.write("source_a\tpathway_a\tsource_b\tpathway_b\talign_score\talign_ratio\tjaccard\n")
        for m in (m for m in matches if m.decision == "accept"):
            ov = "" if m.overlap is None else f"{m.overlap:.3f}"
            fh.write(f"{m.a[0]}\t{m.a[1]}\t{m.b[0]}\t{m.b[1]}\t{m.score}\t{m.ratio:.4f}\t{ov}\n")
    log(f"  {stats['pathways']}")

    if with_go and "go_gaf" in files:
        log("building GO gene sets")
        terms = golib.parse_obo_cached(files["go_obo"])
        ann = golib.parse_gaf(files["go_gaf"], terms, mapper)
        full = golib.propagate(ann, terms)
        go_sets = golib.build_go_sets(terms, full)
        from .meaning import go_definitions

        defs = go_definitions(files["go_obo"])
        for gs in go_sets:
            gs.member_desc = {t: defs[t] for _src, t in gs.members if t in defs}
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

    net = net_string = None
    if with_ppi:
        # physical PPI: BioGRID + IntAct/MINT + HuRI merged; STRING kept as its own network
        log("integrating physical PPIs (BioGRID, IntAct/MINT, HuRI); STRING separately")
        shared = Path(raw_root) / "shared"
        feeds = []
        if org.string_taxid:
            links, info = sources.fetch_string(raw, org)
            net_string, rep_s = ppilib.integrate([("STRING", ppilib.parse_string(links, info, string_min_score))], mapper)
            net_string.overlay_pathways([s for s in sets if s.collection == "pathway"])
            net_string.write(out / "intpath_string.tsv")
            stats["string"] = {"per_source": rep_s, **net_string.summary()}
            log(f"  STRING: {net_string.summary()['edges']} edges")
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
    if exports:
        stats["files"] = write_release(sets, out)
    else:  # wide builds: the database, GMT and stats only
        from .io import write_gmt

        write_gmt(sets, out / "intpath.gmt")
        stats["files"] = {"gmt": str(out / "intpath.gmt")}
    release_genes = {g for s in sets for g in s.genes}
    for n in (net, net_string):
        if n is not None:
            release_genes |= {g for e in n.edges for g in e}
    aliases = mapper.alias_table(release_genes)
    stats["aliases"] = len(aliases)
    stats["mapping"] = dict(mapper.stats)
    log("writing release database")
    meta = {"organism_key": org.key, "organism": org.name, "taxid": org.taxid, "built": stats["built"],
            "tier": stats["tier"], "versions": stats["versions"], "stats": {k: v for k, v in stats.items() if k != "files"}}
    write_db(out / "intpath.sqlite", sets, meta=meta, aliases=aliases, gene_ids=mapper.entrez, ppi=net, string=net_string,
             matches=matches, msigdb_equivalents=msig.equivalents if msig else None)
    from .diagrams import add_diagrams

    stats["diagrams"] = add_diagrams(out / "intpath.sqlite", raw, org, lambda k: mapper.map(k.split(":", 1)[-1]))
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


# --------------------------------------------------------------------------- #
# Every KEGG organism
# --------------------------------------------------------------------------- #
def _build_one(org: Organism, raw_root: str, out_root: str, public: bool, topology: bool) -> tuple[str, str]:
    """One catalog organism in a worker process; never raises (returns status)."""
    import contextlib
    import io as _io

    out = Path(out_root) / org.key
    buf = _io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            build(org, raw_root, out_root, pathway_sources=() if public else ("KEGG",), with_msigdb=False,
                  topology=topology, public=public, legacy_data_root=None, exports=False)
        status = "ok"
    except Exception as exc:  # one organism must not stop thousands
        status = f"failed: {type(exc).__name__}: {exc}"[:300]
    out.mkdir(parents=True, exist_ok=True)
    (out / "build.log").write_text(buf.getvalue() + "\n" + status + "\n")
    return org.key, status


def build_many(
    raw_root: str | Path,
    out_root: str | Path,
    *,
    public: bool = False,
    workers: int = 8,
    lineage: str | None = None,
    codes: list[str] | None = None,
    limit: int | None = None,
    rebuild: bool = False,
    topology: bool = False,
) -> dict:
    """Build every KEGG organism not in the curated registry (see intpath.catalog).

    Phase 1 (serial): KEGG genome catalog, gene namespaces and, for the full
    tier, each organism's KEGG pathway list and gene links, all through the
    rate-limited KEGG gateway. Phase 2 (parallel processes): the builds, which
    then touch only cached KEGG files.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    from . import catalog

    raw_root, out_root = Path(raw_root), Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    rows = catalog.build_catalog(raw_root / "shared")
    catalog.register(rows)
    catalog.write_registry(rows, out_root / "organisms.json")
    orgs = catalog.organisms_from_catalog(rows)
    if lineage:
        keep = {r["kegg"] for r in rows if r["lineage"].lower().startswith(lineage.lower())}
        orgs = [o for o in orgs if o.kegg in keep]
    if codes:
        orgs = [o for o in orgs if o.kegg in set(codes)]
    if not rebuild:
        orgs = [o for o in orgs if not (out_root / o.key / "intpath.sqlite").exists()]
    if limit:
        orgs = orgs[:limit]
    log(f"{len(orgs)} organisms to build ({'open' if public else 'full'} tier)")

    from concurrent.futures import ThreadPoolExecutor

    def prefetch(org: Organism) -> Organism | None:
        raw = raw_root / org.key
        try:
            sources.gene_namespace(raw, org)
            if not public:
                sources.kegg_get(f"list/pathway/{org.kegg}", raw / "kegg" / org.kegg / "list_pathway.txt")
                sources.kegg_get(f"link/{org.kegg}/pathway", raw / "kegg" / org.kegg / "link_pathway.txt")
            return org
        except Exception as exc:
            log(f"  prefetch failed for {org.key}: {exc}")
            return None

    # the NCBI group files are split by the first organism of each group; do that serially
    for group in sorted({o.gene_info.split("/")[0] for o in orgs if o.gene_info_is_shared}):
        from .catalog import gene_info_for_taxid

        gene_info_for_taxid(raw_root / "shared", group, "0")
    ready: list[Organism] = []
    with ThreadPoolExecutor(max_workers=3) as pool:  # kegg_get keeps KEGG within its rate limit
        for i, org in enumerate(pool.map(prefetch, orgs), 1):
            if org is not None:
                ready.append(org)
            if i % 250 == 0:
                log(f"  prefetched {i}/{len(orgs)}")

    results: dict[str, str] = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(_build_one, o, str(raw_root), str(out_root), public, topology) for o in ready]
        for i, f in enumerate(as_completed(futs), 1):
            key, status = f.result()
            results[key] = status
            if i % 100 == 0 or status != "ok":
                log(f"  [{i}/{len(ready)}] {key}: {status}")
    summary = {"built": sum(1 for s in results.values() if s == "ok"),
               "failed": {k: v for k, v in results.items() if v != "ok"}}
    (out_root / f"build_many_{time.strftime('%Y%m%d_%H%M')}.json").write_text(json.dumps(summary, indent=2))
    log(f"done: {summary['built']} built, {len(summary['failed'])} failed")
    return summary
