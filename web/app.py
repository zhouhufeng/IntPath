"""IntPathV2 web service (intpath.genohub.org).

    uvicorn web.app:app            # serves Data/intpathv2/release (env INTPATH_RELEASE_ROOT)
    intpath serve --release-root Data/intpathv2/release

REST API (all organisms that have a built release directory):
    GET  /api/organisms
    GET  /api/{org}/stats
    GET  /api/{org}/search?q=...           pathways / GO terms by name or gene
    GET  /api/{org}/set/{set_id}           genes (+source provenance), members, pairs, links
    GET  /api/{org}/gene/{symbol}          sets containing the gene, PPI partners
    POST /api/{org}/enrich/ora             {"genes": [...], "background": [...]?, "collections": [...]?}
    POST /api/{org}/enrich/pairs           {"genes": [...], "with_ppi": true}
    POST /api/{org}/enrich/gsea            {"ranking": {"GENE": score, ...}, "nperm": 1000}
    GET  /download/{org}/{file}            release files
"""

from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from pydantic import BaseModel, Field

from intpath import __version__, enrich
from intpath import organisms as orglib
from intpath.io import read_release
from intpath.ppi import PPINetwork, canon

STATIC = Path(__file__).parent / "static"
MAX_GENES = 20000


class GeneListRequest(BaseModel):
    genes: list[str] = Field(..., max_length=MAX_GENES)
    background: list[str] | None = None
    collections: list[str] | None = None
    min_size: int = 5
    max_size: int = 2000
    alpha: float = 0.05
    with_ppi: bool = False


class RankingRequest(BaseModel):
    ranking: dict[str, float]
    collections: list[str] | None = None
    nperm: int = Field(1000, ge=100, le=10000)
    min_size: int = 15
    max_size: int = 500
    alpha: float = 0.05
    seed: int = 0


def create_app(release_root: str | Path | None = None) -> FastAPI:
    root = Path(release_root or os.environ.get("INTPATH_RELEASE_ROOT", "Data/intpathv2/release")).resolve()
    app = FastAPI(title="IntPath", version=__version__, description="Integrated pathways, PPIs, GO and enrichment")

    def available() -> list[str]:
        return sorted(p.name for p in root.iterdir() if (p / "intpath_genesets.tsv").exists()) if root.exists() else []

    @lru_cache(maxsize=8)
    def load(org: str):
        if org not in available():
            raise HTTPException(404, f"no release for organism {org!r}")
        sets = read_release(root / org)
        by_id = {s.id: s for s in sets}
        by_gene: dict[str, list[str]] = {}
        for s in sets:
            for g in s.genes:
                by_gene.setdefault(g, []).append(s.id)
        ppi_f = root / org / "intpath_ppi.tsv"
        ppi = PPINetwork.read(ppi_f) if ppi_f.exists() else None
        upper: dict[str, str] = {}
        alias_f = root / org / "intpath_gene_aliases.tsv"
        if alias_f.exists():  # ids, previous symbols and unambiguous aliases from the build's GeneMapper
            with open(alias_f) as fh:
                next(fh)
                for line in fh:
                    k, v = line.rstrip("\n").split("\t")
                    upper[k] = v
        upper.update({g.upper(): g for g in by_gene})
        if ppi is not None:
            upper.update({g.upper(): g for e in ppi.edges for g in e})
        return sets, by_id, by_gene, ppi, upper

    @lru_cache(maxsize=8)
    def neighbours(org: str) -> dict[str, set[str]]:
        ppi = load(org)[3]
        return ppi.neighbours() if ppi is not None else {}

    def resolve(org: str, org_genes: list[str]) -> tuple[list[str], list[str]]:
        upper = load(org)[4]
        found, missing = [], []
        for g in org_genes:
            (found.append(upper[g.strip().upper()]) if g.strip().upper() in upper else missing.append(g))
        return found, missing

    def org_key(org: str) -> str:
        try:
            return orglib.get(org).key
        except KeyError:
            return org

    @app.get("/", response_class=HTMLResponse)
    def index():
        return (STATIC / "index.html").read_text()

    @app.get("/api/organisms")
    def organisms():
        avail = set(available())
        return [
            {"key": o.key, "name": o.name, "taxid": o.taxid, "available": o.key in avail, "in_old_intpath": o.in_old_intpath}
            for o in orglib.ORGANISMS.values()
        ]

    @app.get("/api/{org}/stats")
    def stats(org: str):
        f = root / org_key(org) / "stats.json"
        if not f.exists():
            raise HTTPException(404, "no stats")
        return json.loads(f.read_text())

    @app.get("/api/{org}/search")
    def search(org: str, q: str, limit: int = 50):
        sets, _, by_gene, _, upper = load(org_key(org))
        ql = q.strip().lower()
        hits = [s for s in sets if ql in s.name.lower() or ql == s.id.lower()]
        gene = upper.get(q.strip().upper())
        if gene:
            ids = set(by_gene.get(gene, ()))
            hits += [s for s in sets if s.id in ids and s not in hits]
        return [
            {"id": s.id, "name": s.name, "collection": s.collection, "sources": s.sources, "n_genes": s.size}
            for s in hits[:limit]
        ]

    @app.get("/api/{org}/set/{set_id:path}")
    def get_set(org: str, set_id: str):
        _, by_id, *_ = load(org_key(org))
        s = by_id.get(set_id)
        if s is None:
            raise HTTPException(404, "unknown set")
        return {
            "id": s.id,
            "name": s.name,
            "collection": s.collection,
            "members": [{"source": a, "pathway": b} for a, b in s.members],
            "genes": {g: sorted(v) for g, v in sorted(s.genes.items())},
            "pairs": [
                {"a": a, "b": b, "relations": sorted(e["rel"]), "sources": sorted(e["src"])}
                for (a, b), e in sorted(s.pairs.items())
            ],
            "links": [{"id": i, "name": by_id[i].name, "collection": by_id[i].collection} for i in s.links if i in by_id],
        }

    @app.get("/api/{org}/gene/{symbol}")
    def gene(org: str, symbol: str):
        key = org_key(org)
        _, by_id, by_gene, ppi, upper = load(key)
        g = upper.get(symbol.upper())
        if g is None:
            raise HTTPException(404, "gene not in IntPath")
        partners = []
        for h in neighbours(key).get(g, ()):
            ev = ppi.edges[canon(g, h)]
            partners.append({"gene": h, "sources": sorted(ev.sources),
                             "string_score": ev.string_score, "pathway_supported": bool(ev.pathways)})
        return {
            "gene": g,
            "sets": [{"id": i, "name": by_id[i].name, "collection": by_id[i].collection,
                      "sources": sorted(by_id[i].genes[g])} for i in by_gene.get(g, [])],
            "ppi_partners": sorted(partners, key=lambda r: (-len(r["sources"]), -r["string_score"])),
        }

    @app.post("/api/{org}/enrich/ora")
    def run_ora(org: str, req: GeneListRequest):
        key = org_key(org)
        sets = load(key)[0]
        genes, missing = resolve(key, req.genes)
        bg = resolve(key, req.background)[0] if req.background else None
        rows = enrich.ora(genes, sets, background=bg, collections=req.collections, min_size=req.min_size,
                          max_size=req.max_size)
        enrich.cluster(rows, sets, alpha=req.alpha)
        return {"n_input": len(req.genes), "n_mapped": len(genes), "unmapped": missing[:200], "results": rows[:500]}

    @app.post("/api/{org}/enrich/pairs")
    def run_pairs(org: str, req: GeneListRequest):
        key = org_key(org)
        sets, _, _, ppi, _ = load(key)
        genes, missing = resolve(key, req.genes)
        rows = enrich.pair_ora(genes, sets, ppi=ppi if req.with_ppi else None,
                               collections=req.collections or ("pathway",))
        return {"n_input": len(req.genes), "n_mapped": len(genes), "unmapped": missing[:200], "results": rows[:500]}

    @app.post("/api/{org}/enrich/gsea")
    def run_gsea(org: str, req: RankingRequest):
        key = org_key(org)
        sets = load(key)[0]
        upper = load(key)[4]
        ranking = {upper[g.upper()]: v for g, v in req.ranking.items() if g.upper() in upper}
        rows = enrich.gsea(ranking, sets, collections=req.collections, nperm=req.nperm, min_size=req.min_size,
                           max_size=req.max_size, seed=req.seed)
        enrich.cluster(rows, sets, alpha=req.alpha)
        return {"n_input": len(req.ranking), "n_mapped": len(ranking), "results": rows[:500]}

    @app.get("/download/{org}/{name}")
    def download(org: str, name: str):
        f = (root / org_key(org) / name).resolve()
        if f.parent != (root / org_key(org)).resolve() or not f.is_file():
            raise HTTPException(404, "no such file")
        return FileResponse(f, filename=f"{org}_{name}")

    return app


app = create_app()
