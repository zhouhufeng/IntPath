"""Command line interface: ``intpath <command> ...`` (or ``python -m intpath``)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import organisms


def _read_list(path: str) -> list[str]:
    text = sys.stdin.read() if path == "-" else Path(path).read_text()
    return [t for t in text.replace(",", " ").split() if t]


def _read_ranking(path: str) -> dict[str, float]:
    out = {}
    for line in Path(path).read_text().splitlines():
        f = line.replace(",", "\t").split("\t")
        if len(f) >= 2:
            try:
                out[f[0].strip()] = float(f[1])
            except ValueError:
                continue  # header
    return out


def _write(rows: list[dict], path: str | None, limit: int | None = None) -> None:
    rows = rows[:limit] if limit else rows
    fh = open(path, "w") if path else sys.stdout
    if not rows:
        fh.write("# no results\n")
        return
    cols = [c for c in rows[0] if not isinstance(rows[0][c], list)] + [c for c in rows[0] if isinstance(rows[0][c], list)]
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        vals = []
        for c in cols:
            v = r.get(c)
            vals.append(",".join(v) if isinstance(v, list) else (f"{v:.4g}" if isinstance(v, float) else str(v)))
        fh.write("\t".join(vals) + "\n")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="intpath", description="IntPathV2 - integrated pathways, PPIs, GO and enrichment")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("organisms", help="list supported organisms")

    b = sub.add_parser("build", help="download sources and build a release for an organism")
    b.add_argument("organism")
    b.add_argument("--raw", default="Data/intpathv2/raw")
    b.add_argument("--out", default="Data/intpathv2/release")
    b.add_argument("--sources", default="KEGG,Reactome,WikiPathways")
    b.add_argument("--biocyc-col", help="licensed BioCyc pathways.col file")
    b.add_argument("--ppi", action="append", default=[], metavar="LABEL:FORMAT:PATH",
                   help="extra PPI file, FORMAT in biogrid|mitab|tsv (repeatable)")
    b.add_argument("--string-min-score", type=int, default=700)
    b.add_argument("--no-go", action="store_true")
    b.add_argument("--no-ppi", action="store_true")
    b.add_argument("--no-msigdb", action="store_true")
    b.add_argument("--no-topology", action="store_true", help="skip KGML/GPML gene-pair extraction")
    b.add_argument("--public", action="store_true", help="exclude licence-restricted sources (KEGG, BioCyc)")

    bm = sub.add_parser("build-many", help="build every KEGG organism outside the curated registry")
    bm.add_argument("--raw", default="Data/intpathv2/raw")
    bm.add_argument("--out", default=None, help="default: Data/intpathv2/release[-open]")
    bm.add_argument("--public", action="store_true", help="open tier: no KEGG (GO + PPIs only)")
    bm.add_argument("--workers", type=int, default=8)
    bm.add_argument("--lineage", help="KEGG lineage prefix, e.g. 'Eukaryotes' or 'Prokaryotes;Bacteria'")
    bm.add_argument("--codes", help="comma list of KEGG organism codes")
    bm.add_argument("--limit", type=int)
    bm.add_argument("--rebuild", action="store_true")
    bm.add_argument("--topology", action="store_true", help="also download KGML gene pairs (slow: ~150-400 KEGG calls/organism)")

    lg = sub.add_parser("legacy", help="re-run the 2012 merge on archived old IntPath normalized files")
    lg.add_argument("--data", default="Data")
    lg.add_argument("--out", default="Data/intpathv2/release")
    lg.add_argument("--organism", default="sapiens")

    for name, hlp in (("ora", "over-representation analysis"), ("pairs", "gene-pair (network) enrichment"),
                      ("gsea", "preranked GSEA")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("release", help="release directory, e.g. Data/intpathv2/release/sapiens")
        p.add_argument("input", help="gene list (ora/pairs) or 2-column ranking file (gsea); '-' for stdin")
        p.add_argument("--collections", help="comma list: pathway,GO:BP,GO:MF,GO:CC (default all)")
        p.add_argument("--out")
        p.add_argument("--top", type=int)
        p.add_argument("--alpha", type=float, default=0.05)
        if name == "ora":
            p.add_argument("--background")
        if name == "pairs":
            p.add_argument("--with-ppi", action="store_true", help="add merged PPI edges inside each pathway")
        if name == "gsea":
            p.add_argument("--nperm", type=int, default=1000)
            p.add_argument("--seed", type=int, default=0)

    s = sub.add_parser("serve", help="run the web/API server")
    s.add_argument("--release-root", default="Data/intpathv2/release")
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--port", type=int, default=8000)

    a = ap.parse_args(argv)

    if a.cmd == "organisms":
        for o in organisms.ORGANISMS.values():
            flag = " (in old IntPath)" if o.in_old_intpath else ""
            print(f"{o.key:14s}{o.name:36s} taxid={o.taxid:7s} kegg={o.kegg or '-':4s} "
                  f"reactome={'y' if o.reactome else '-'} wp={'y' if o.wikipathways else '-'} go={o.go_gaf or '-'}{flag}")
        return

    if a.cmd == "build":
        from .build import build

        extra = tuple(tuple(x.split(":", 2)) for x in a.ppi)
        stats = build(organisms.get(a.organism), a.raw, a.out, pathway_sources=tuple(a.sources.split(",")),
                      biocyc_col=a.biocyc_col, with_go=not a.no_go, with_ppi=not a.no_ppi,
                      string_min_score=a.string_min_score, extra_ppi=extra, topology=not a.no_topology,
                      public=a.public, with_msigdb=not a.no_msigdb)
        print(json.dumps({k: v for k, v in stats.items() if k != "files"}, indent=2))
        return

    if a.cmd == "build-many":
        from .build import build_many

        out = a.out or ("Data/intpathv2/release-open" if a.public else "Data/intpathv2/release")
        print(json.dumps(build_many(a.raw, out, public=a.public, workers=a.workers, lineage=a.lineage,
                                    codes=a.codes.split(",") if a.codes else None, limit=a.limit,
                                    rebuild=a.rebuild, topology=a.topology), indent=2)[:2000])
        return

    if a.cmd == "legacy":
        from .build import rebuild_legacy

        print(json.dumps(rebuild_legacy(a.data, a.out, a.organism), indent=2))
        return

    if a.cmd == "serve":
        import uvicorn

        from web.app import create_app

        uvicorn.run(create_app(a.release_root), host=a.host, port=a.port)
        return

    from . import enrich
    from .io import read_release

    sets = read_release(a.release)
    cols = a.collections.split(",") if a.collections else None
    if a.cmd == "ora":
        bg = _read_list(a.background) if a.background else None
        rows = enrich.ora(_read_list(a.input), sets, background=bg, collections=cols)
        enrich.cluster(rows, sets, alpha=a.alpha)
    elif a.cmd == "pairs":
        ppi = None
        if a.with_ppi and (Path(a.release) / "intpath_ppi.tsv").exists():
            from .ppi import PPINetwork

            ppi = PPINetwork.read(Path(a.release) / "intpath_ppi.tsv")
        rows = enrich.pair_ora(_read_list(a.input), sets, ppi=ppi, collections=cols or ("pathway",))
    else:
        rows = enrich.gsea(_read_ranking(a.input), sets, collections=cols, nperm=a.nperm, seed=a.seed)
        enrich.cluster(rows, sets, alpha=a.alpha, fdr_key="fdr")
    _write(rows, a.out, a.top)


if __name__ == "__main__":
    main()
