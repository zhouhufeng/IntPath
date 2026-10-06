"""IntPathV2 web service (intpath.genohub.org).

    uvicorn web.app:app               # serves $INTPATH_RELEASE_ROOT/<organism>/intpath.sqlite
    intpath serve --release-root Data/intpathv2/release

Each organism release is one read-only SQLite file (intpath.sqlite). Statistics
run on a compact in-memory library (intpath.library); details (provenance, gene
pairs, PPI evidence) are read from SQLite per request.

REST API
    GET  /healthz
    GET  /api/organisms
    GET  /api/{org}/stats                      build statistics and source versions
    GET  /api/{org}/collections                gene set collections with counts
    GET  /api/{org}/search?q=...               sets by name, id or gene
    GET  /api/{org}/set/{set_id}               genes (+sources), members, gene pairs, links, MSigDB equivalents
    GET  /api/{org}/gene/{symbol}              sets containing the gene, PPI partners with evidence
    GET  /api/{org}/map/{set_id}               pathway map: ?format=cyjs|sbml|sbgn &ppi=high|medium|low
                                               &genes=<your list> (highlighted)
    POST /api/{org}/enrich/ora                 {"genes": [...], "background": [...]?, "collections": [...]?}
    POST /api/{org}/enrich/pairs               {"genes": [...], "with_ppi": true, "ppi_tier": "high"}
    POST /api/{org}/enrich/gsea                {"ranking": {"GENE": score}, "nperm": 1000}
    GET  /download/{org}/{file}                release files

Licensed section (/licensed/...): the same API over the full-tier releases
(KEGG, BioCyc, MSigDB BioCarta/KEGG_MEDICUS), for signed-in users only. The
gateway authenticates every /licensed/ request with the shared sign-in gate and
passes X-IGVF-User; requests without it are refused here as well. Downloads are
disabled there: licensed data is offered for analysis, not redistribution.
"""

from __future__ import annotations

import json
import os
import re
import threading
from functools import lru_cache
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from pydantic import BaseModel, Field

from intpath import __version__, diagrams, enrich, maps, sources
from intpath import organisms as orglib
from intpath.library import Library, load_library

STATIC = Path(__file__).parent / "static"
MAX_GENES = 20000
DOWNLOADS = ("intpath.gmt", "intpath_genesets.tsv", "intpath_set_genes.tsv", "intpath_set_genepairs.tsv",
             "intpath_set_members.tsv", "intpath_ppi.tsv", "related_pathways.tsv", "stats.json")


class GeneListRequest(BaseModel):
    genes: list[str] = Field(..., max_length=MAX_GENES)
    background: list[str] | None = Field(None, max_length=60000)
    collections: list[str] | None = None
    min_size: int = Field(5, ge=1)
    max_size: int = Field(2000, le=5000)
    alpha: float = Field(0.05, gt=0, le=1)
    with_ppi: bool = False
    ppi_tier: str = Field("high", pattern="^(high|medium|low|string)$")


class RankingRequest(BaseModel):
    ranking: dict[str, float] = Field(..., max_length=60000)
    collections: list[str] | None = None
    nperm: int = Field(1000, ge=100, le=int(os.environ.get("INTPATH_MAX_NPERM", "2000")))
    min_size: int = Field(15, ge=1)
    max_size: int = Field(500, le=5000)
    alpha: float = Field(0.05, gt=0, le=1)
    seed: int = 0


def _register_generated(root: Path) -> None:
    """Make catalog organisms (every KEGG genome, see intpath.catalog) resolvable by key/code/taxid."""
    f = root / "organisms.json"
    if not f.exists():
        return
    for d in json.loads(f.read_text()):
        if d["key"] in orglib.ORGANISMS:
            continue
        d.pop("lineage", None)
        d["aliases"] = tuple(d.get("aliases", ()))
        orglib.ORGANISMS[d["key"]] = orglib.Organism(**d)


def create_app(release_root: str | Path | None = None, *, licensed: bool = False,
               licensed_root: str | Path | None = None) -> FastAPI:
    root = Path(release_root or os.environ.get("INTPATH_RELEASE_ROOT", "Data/intpathv2/db")).resolve()
    # one unified database (all sources licensed by the operator): KEGG drawings etc. served directly
    unified = os.environ.get("INTPATH_UNIFIED", "1") == "1"
    heavy = threading.BoundedSemaphore(int(os.environ.get("INTPATH_MAX_JOBS", "2")))  # concurrent GSEA runs
    app = FastAPI(title="IntPath (licensed)" if licensed else "IntPath", version=__version__,
                  description="IntPathV2: integrated pathways, PPIs, GO and MSigDB with gene set enrichment")

    if licensed:
        # Local trial only: INTPATH_DEV_USER stands in for the sign-in gate on a workstation.
        # Deploy/docker-compose.yml never sets it; behind the gateway the gate decides.
        dev_user = os.environ.get("INTPATH_DEV_USER", "").strip()

        @app.middleware("http")
        async def require_signed_in(request, call_next):
            if not request.headers.get("x-igvf-user") and not dev_user:
                return JSONResponse({"error": "sign in at https://intpath.genohub.org/licensed/"}, status_code=401)
            return await call_next(request)

    _register_generated(root)

    avail_cache: dict = {"t": 0.0, "v": []}

    def available() -> list[str]:
        # thousands of release directories: rescan at most once a minute
        import time as _t

        if _t.time() - avail_cache["t"] > 60:
            avail_cache["v"] = sorted(p.name for p in root.iterdir() if (p / "intpath.sqlite").exists()) \
                if root.exists() else []
            avail_cache["t"] = _t.time()
        return avail_cache["v"]

    def org_key(org: str) -> str:
        try:
            return orglib.get(org).key
        except KeyError:
            return org

    load_lock = threading.Lock()

    @lru_cache(maxsize=12)
    def _load_cached(key: str) -> Library:
        return load_library(root / key / "intpath.sqlite", key)

    def _load(key: str) -> Library:
        with load_lock:  # the preload thread and a first request must not load twice
            return _load_cached(key)

    def lib(org: str) -> Library:
        key = org_key(org)
        if key not in available():
            raise HTTPException(404, f"no release for organism {org!r}")
        return _load(key)

    def query(org: str, sql: str, args=()) -> list[dict]:
        con = lib(org).connect()
        try:
            return [dict(r) for r in con.execute(sql, args)]
        finally:
            con.close()

    def read_stats(key: str) -> dict:
        f = root / key / "stats.json"
        return json.loads(f.read_text()) if f.exists() else {}

    @app.get("/", response_class=HTMLResponse)
    def index():
        # the page changes with every release: never let the browser reuse an old copy
        return HTMLResponse((STATIC / "index.html").read_text(), headers={"Cache-Control": "no-store"})

    @app.get("/static/{name}")
    def static_file(name: str):
        if name not in ("logo.png", "icon.png"):
            raise HTTPException(404, "not found")
        return FileResponse(STATIC / name, media_type="image/png", headers={"Cache-Control": "public, max-age=86400"})

    @app.get("/healthz")
    def healthz():
        return {"ok": True, "organisms": available(), "version": __version__, "licensed": licensed}

    @app.get("/api/organisms")
    def organisms(all: bool = False):  # noqa: A002
        """Organisms with a release (all=true also lists registered organisms without one)."""
        avail = set(available())
        out = []
        for o in orglib.ORGANISMS.values():
            if o.key not in avail and not all:
                continue
            out.append({"key": o.key, "name": o.name, "taxid": o.taxid, "kegg": o.kegg,
                        "available": o.key in avail, "in_old_intpath": o.in_old_intpath,
                        "curated": o.key in orglib.CURATED})
        out.sort(key=lambda r: (not r["curated"], r["name"].lower()))
        return out

    @app.get("/api/{org}/stats")
    def stats(org: str):
        st = read_stats(org_key(org))
        if not st:
            raise HTTPException(404, "no stats")
        st.pop("files", None)
        return st

    @app.get("/api/{org}/collections")
    def collections(org: str):
        return query(org, "SELECT collection, COUNT(*) AS n FROM gset GROUP BY collection ORDER BY collection")

    @app.get("/api/{org}/search")
    def search(org: str, q: str = "", limit: int = 50):
        term = q.strip()
        if not term:
            return []
        L = lib(org)
        limit = min(max(limit, 1), 200)
        rows = query(org, "SELECT set_id AS id, name, collection, sources, n_genes FROM gset "
                          "WHERE name LIKE ? OR set_id = ? "
                          "OR set_id IN (SELECT set_id FROM set_member WHERE name LIKE ?) ORDER BY "
                          "CASE WHEN lower(name) = lower(?) THEN 0 ELSE 1 END, "
                          "CASE WHEN collection = 'pathway' THEN 0 WHEN collection LIKE 'GO:%' THEN 1 "
                          "WHEN collection = 'msigdb:H' THEN 2 ELSE 3 END, length(name) LIMIT ?",
                     (f"%{term}%", term, f"%{term}%", term, limit))
        gene = L.aliases.get(term.upper())
        if gene:
            seen = {r["id"] for r in rows}
            for sid in L.by_gene.get(gene, [])[: max(0, limit - len(rows))]:
                if sid not in seen:
                    s = L.by_id[sid]
                    rows.append({"id": s.id, "name": s.name, "collection": s.collection,
                                 "sources": ",".join(s.sources), "n_genes": s.size})
        for r in rows:
            r["sources"] = [x for x in (r["sources"] or "").split(",") if x]
        return rows

    @app.get("/api/{org}/set/{set_id:path}")
    def get_set(org: str, set_id: str):
        head = query(org, "SELECT set_id AS id, name, collection, sources FROM gset WHERE set_id = ?", (set_id,))
        if not head:
            raise HTTPException(404, "unknown set")
        s = head[0]
        s["sources"] = [x for x in (s["sources"] or "").split(",") if x]
        s["members"] = query(org, "SELECT *, source_set_id AS pathway FROM set_member WHERE set_id = ?", (set_id,))
        s["genes"] = {r["symbol"]: r["sources"].split(",") for r in
                      query(org, "SELECT symbol, sources FROM set_gene WHERE set_id = ? ORDER BY symbol", (set_id,))}
        s["pairs"] = [{"a": r["gene_a"], "b": r["gene_b"], "relations": r["relations"].split(","),
                       "sources": r["sources"].split(",")}
                      for r in query(org, "SELECT * FROM set_pair WHERE set_id = ? LIMIT 2000", (set_id,))]
        s["links"] = query(org, "SELECT g.set_id AS id, g.name, g.collection FROM set_link l JOIN gset g "
                                "ON g.set_id = l.set_b WHERE l.set_a = ?", (set_id,))
        s["msigdb_equivalents"] = query(org, "SELECT msigdb_name, systematic_name, collection, jaccard "
                                             "FROM msigdb_equivalent WHERE set_id = ?", (set_id,))
        return s

    def _with_hits(d: dict, hits: set[str], focus: str = "") -> dict:
        for n in d["nodes"]:
            n["hit"] = any(g in hits for g in n.get("genes", ()))
        keep = set(d.get("subpathways", {}).get(focus, ())) if focus else set()
        if keep:  # a sub-pathway drawn inside an ancestor's diagram: flag its part
            for n in d["nodes"]:
                n["focus"] = n["id"] in keep or n["kind"] in ("shape", "label", "title")
        d.pop("subpathways", None)
        return d

    @app.get("/api/{org}/diagrams/{set_id:path}")
    def diagrams_of(org: str, set_id: str):
        """Drawings available for a gene set: KEGG / WikiPathways (positioned), Reactome (embedded viewer)."""
        members = query(org, "SELECT * FROM set_member WHERE set_id = ?", (set_id,))
        try:
            stored = {(r["source"], r["source_set_id"]): r["name"] for r in
                      query(org, "SELECT source, source_set_id, name FROM diagram")}
            drawn_in = {r["source_set_id"]: r["diagram_id"] for r in
                        query(org, "SELECT source_set_id, diagram_id FROM diagram_of WHERE source = 'Reactome'")}
        except Exception:  # release built before diagrams existed
            stored, drawn_in = {}, {}
        out = []
        for m in members:
            key = (m["source"], m["source_set_id"])
            if m["source"] == "Reactome":
                diag = drawn_in.get(key[1])
                if diag and ("Reactome", diag) in stored:
                    out.append({"source": "Reactome", "id": diag, "focus": key[1] if diag != key[1] else "",
                                "name": stored[("Reactome", diag)], "member_name": m.get("name"), "kind": "drawing"})
            elif key in stored:
                out.append({"source": key[0], "id": key[1], "name": stored[key], "member_name": m.get("name"),
                            "kind": "drawing"})
            elif m["source"] == "KEGG" and (licensed or unified):
                out.append({"source": "KEGG", "id": key[1], "name": key[1], "member_name": m.get("name"),
                            "kind": "drawing"})  # fetched on demand
        order = {"KEGG": 0, "WikiPathways": 1, "Reactome": 2}
        return sorted(out, key=lambda d: order.get(d["source"], 9))

    @app.get("/api/{org}/diagram/{source}/{sid}")
    def diagram(org: str, source: str, sid: str, genes: str = "", focus: str = ""):
        """One positioned pathway drawing, with the given genes flagged ``hit``."""
        L = lib(org)
        hits = set(L.resolve([g for g in genes.replace(",", " ").split() if g][:20000])[0]) if genes else set()
        con = L.connect()
        try:
            row = con.execute("SELECT data FROM diagram WHERE source = ? AND source_set_id = ?", (source, sid)).fetchone()
        except Exception:
            row = None
        finally:
            con.close()
        if row is not None:
            return _with_hits(diagrams.unpack(row[0]), hits, focus)
        if source == "KEGG" and (licensed or unified) and re.fullmatch(r"[a-z]{2,4}\d{5}", sid):
            # organisms built without KGML: fetch the drawing once through the rate-limited KEGG gateway
            cache = root / L.organism / "kgml" / f"{sid}.xml"
            try:
                sources.kegg_get(f"get/{sid}/kgml", cache)
            except Exception:
                raise HTTPException(404, "KEGG drawing unavailable") from None
            d = diagrams.kgml_diagram(cache.read_bytes(), lambda k: L.aliases.get(k.split(":", 1)[-1].upper()))
            if d is None:
                raise HTTPException(404, "no drawing for this KEGG map")
            return _with_hits(d, hits)
        raise HTTPException(404, "no drawing")

    @app.get("/api/{org}/map/{set_id:path}")
    def pathway_map(org: str, set_id: str, format: str = "cyjs", ppi: str = "",  # noqa: A002
                    genes: str = "", max_nodes: int = Query(maps.MAX_NODES, alias="max")):
        """Map of one gene set: Cytoscape.js JSON (default), SBML-qual (format=sbml) or SBGN-ML (format=sbgn).

        ``genes``: the user's gene list (comma/space separated); those genes are flagged ``hit``.
        ``ppi``: high | medium | low adds merged PPI edges among the set's genes.
        """
        if ppi not in ("", "high", "medium", "low", "string"):
            raise HTTPException(400, "ppi must be high, medium, low (physical PPI) or string")
        if format not in ("cyjs", "sbml", "sbgn"):
            raise HTTPException(400, "format must be cyjs, sbml or sbgn")
        L = lib(org)
        wanted = [g for g in genes.replace(",", " ").split() if g][:20000]
        hits = set(L.resolve(wanted)[0]) if wanted else set()
        con = L.connect()
        try:
            graph = maps.set_graph(con, set_id, ppi_tier=ppi or None, highlight=hits,
                                   adjacency=L.adjacency.get(ppi) if ppi else None,
                                   max_nodes=min(max(max_nodes, 10), maps.MAX_NODES) if format == "cyjs" else maps.MAX_NODES)
        except KeyError:
            raise HTTPException(404, "unknown set") from None
        finally:
            con.close()
        stem = "".join(c if c.isalnum() else "_" for c in graph["name"])[:60] or "pathway"
        if format == "sbml":
            return Response(maps.to_sbml_qual(graph), media_type="application/xml",
                            headers={"Content-Disposition": f'attachment; filename="{stem}.sbml.xml"'})
        if format == "sbgn":
            return Response(maps.to_sbgn_af(graph), media_type="application/xml",
                            headers={"Content-Disposition": f'attachment; filename="{stem}.sbgn"'})
        return maps.to_cytoscape(graph)

    @app.get("/api/{org}/gene/{symbol}")
    def gene(org: str, symbol: str, limit: int = 300):
        L = lib(org)
        g = L.aliases.get(symbol.strip().upper())
        if g is None:
            raise HTTPException(404, "gene not in IntPath")
        sets = query(org, "SELECT g.set_id AS id, g.name, g.collection, sg.sources FROM set_gene sg JOIN gset g "
                          "ON g.set_id = sg.set_id WHERE sg.symbol = ? ORDER BY g.collection, g.name", (g,))
        lim = min(max(limit, 1), 2000)
        partners = query(org, "SELECT CASE WHEN gene_a = ? THEN gene_b ELSE gene_a END AS gene, sources, n_pmids, "
                              "string_score, tier, pathway_sets != '' AS pathway_supported FROM ppi "
                              "WHERE gene_a = ? OR gene_b = ? ORDER BY n_sources DESC, n_pmids DESC LIMIT ?",
                         (g, g, g, lim))
        try:
            string = query(org, "SELECT CASE WHEN gene_a = ? THEN gene_b ELSE gene_a END AS gene, string_score "
                                "FROM string_ppi WHERE gene_a = ? OR gene_b = ? ORDER BY string_score DESC LIMIT ?",
                           (g, g, g, lim))
        except Exception:  # releases built before STRING was separate
            string = []
        for r in sets + partners:
            r["sources"] = r["sources"].split(",")
        return {"gene": g, "sets": sets, "ppi_partners": partners, "string_partners": string}

    @app.post("/api/{org}/enrich/ora")
    def run_ora(org: str, req: GeneListRequest):
        L = lib(org)
        genes, missing = L.resolve(req.genes)
        bg = L.resolve(req.background)[0] if req.background else None
        rows = enrich.ora(genes, L.sets, background=bg, collections=req.collections, min_size=req.min_size,
                          max_size=req.max_size, index=L.set_index(), universe=L.universe(req.collections))
        enrich.cluster(rows, L.by_id, alpha=req.alpha)
        return {"n_input": len(req.genes), "n_mapped": len(genes), "unmapped": missing[:200], "results": rows[:500]}

    @app.post("/api/{org}/enrich/pairs")
    def run_pairs(org: str, req: GeneListRequest):
        L = lib(org)
        genes, missing = L.resolve(req.genes)
        rows = enrich.pair_ora(genes, L.sets, ppi=L.adjacency[req.ppi_tier] if req.with_ppi else None,
                               collections=req.collections or ("pathway",))
        return {"n_input": len(req.genes), "n_mapped": len(genes), "unmapped": missing[:200], "results": rows[:500]}

    @app.post("/api/{org}/enrich/gsea")
    def run_gsea(org: str, req: RankingRequest):
        L = lib(org)
        ranking: dict[str, float] = {}
        for g, v in req.ranking.items():
            s = L.aliases.get(g.strip().upper())
            if s and s not in ranking:
                ranking[s] = v
        if not heavy.acquire(timeout=120):
            raise HTTPException(503, "server busy, retry shortly")
        try:
            rows = enrich.gsea(ranking, L.sets, collections=req.collections, nperm=req.nperm,
                               min_size=req.min_size, max_size=req.max_size, seed=req.seed)
        finally:
            heavy.release()
        enrich.cluster(rows, L.by_id, alpha=req.alpha)
        return {"n_input": len(req.ranking), "n_mapped": len(ranking), "results": rows[:500]}

    @app.get("/download/{org}/{name}")
    def download(org: str, name: str):
        if licensed or name not in DOWNLOADS:
            raise HTTPException(404, "no such file")
        f = root / org_key(org) / name
        if not f.is_file():
            raise HTTPException(404, "no such file")
        return FileResponse(f, filename=f"intpath_{org_key(org)}_{name}")

    if os.environ.get("INTPATH_PRELOAD", "1") == "1":
        # load released organisms in the background so the first request is fast
        def preload():
            for o in available():
                try:
                    _load(o)
                except Exception as exc:  # pragma: no cover
                    print(f"[intpath] failed to load {o}: {exc}")

        threading.Thread(target=preload, daemon=True).start()

    if not licensed:
        lic = Path(licensed_root or os.environ.get("INTPATH_LICENSED_ROOT", "")) if (
            licensed_root or os.environ.get("INTPATH_LICENSED_ROOT")) else None
        if lic is not None and lic.exists():
            app.mount("/licensed", create_app(lic, licensed=True))

    return app


app = create_app()
