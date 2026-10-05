"""IntPathV2 enrichment engine.

Methods
-------
ora           Over-representation analysis: one-sided hypergeometric test per
              gene set (the 2012 "Identify Pathways" tool), BH-FDR, fold
              enrichment, plus *source consensus*: for each overlapping gene,
              how many source databases put it in the integrated pathway.
pair_ora      Gene-pair (network) enrichment, unique to IntPath: tests whether
              the curated relations among the query genes (pathway gene pairs,
              optionally plus merged PPI edges) are over-represented in a
              pathway. Detects pathways whose *wiring* - not just membership -
              is hit by the query.
gsea          Preranked GSEA (Subramanian et al. 2005; weighted KS statistic,
              p = 1) with gene-set-size-matched permutation nulls shared
              across sets (fgsea-simple style), NES, nominal p, GSEA FDR, BH.
cluster       Redundancy-aware reporting: significant hits are grouped into
              themes using IntPath merge groups, GO<->pathway links and gene
              overlap, so one biological signal is shown once.

Only the standard library and numpy are required.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np

from .model import GeneSet
from .ppi import PPINetwork, canon


def _members(s) -> frozenset[str]:
    """Gene membership of a GeneSet (build/test objects) or LibSet (compact serving library)."""
    m = getattr(s, "members_set", None)
    return m if m is not None else s.members


# --------------------------------------------------------------------------- #
# Statistics helpers
# --------------------------------------------------------------------------- #
def _lchoose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(population N, K successes, n draws).

    Sums the upper tail with the term ratio
    P(i+1)/P(i) = (K-i)(n-i) / ((i+1)(N-K-n+i+1)), starting from the exact log
    term at the lower bound and stopping once terms no longer change the sum.
    """
    lo, hi = max(k, 0, n - (N - K)), min(K, n)
    if lo > hi:
        return 0.0 if k > hi else 1.0
    if lo == 0:
        return 1.0
    log_first = _lchoose(K, lo) + _lchoose(N - K, n - lo) - _lchoose(N, n)
    term, total = 1.0, 1.0  # relative to the first term
    for i in range(lo, hi):
        term *= (K - i) * (n - i) / ((i + 1) * (N - K - n + i + 1))
        total += term
        if term < total * 1e-17:
            break
    return min(1.0, math.exp(log_first + math.log(total)))


def bh(pvals: Sequence[float]) -> list[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    out = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * n / rank)
        out[i] = prev
    return out


def _select(sets: Iterable, collections: Iterable[str] | None) -> list:
    """Collections match exactly or as a prefix: "msigdb" selects every msigdb:* set,
    "msigdb:C2" selects msigdb:C2:CGP and msigdb:C2:CP:PID, "GO" selects GO:BP/MF/CC."""
    if not collections:
        return list(sets)
    cols = tuple(collections)
    return [s for s in sets if any(s.collection == c or s.collection.startswith(c + ":") for c in cols)]


# --------------------------------------------------------------------------- #
# Over-representation
# --------------------------------------------------------------------------- #
def ora(
    query: Iterable[str],
    sets: Iterable[GeneSet],
    *,
    background: Iterable[str] | None = None,
    collections: Iterable[str] | None = None,
    min_size: int = 5,
    max_size: int = 2000,
    min_overlap: int = 2,
    index: dict[str, list] | None = None,
    universe: set[str] | frozenset[str] | None = None,
) -> list[dict]:
    """Hypergeometric ORA.

    ``index`` (gene -> sets containing it) and ``universe`` (all genes of the
    selected collections) are optional precomputations a server passes in: only
    sets that contain a query gene can reach ``min_overlap``, so with an index
    the others are never scanned. Results are identical either way.
    """
    lib = _select(sets, collections)
    if universe is None:
        universe = {g for s in lib for g in _members(s)}
    if background is not None:
        universe = set(universe) & set(background)
    q = set(query) & universe
    N, n = len(universe), len(q)
    if index is not None:
        selected = {id(s) for s in lib}
        seen: dict[int, object] = {}
        for g in q:
            for s in index.get(g, ()):
                if id(s) in selected:
                    seen[id(s)] = s
        candidates = list(seen.values())
    else:
        candidates = lib
    restricted = background is not None
    rows = []
    for s in candidates:
        mem = _members(s)
        members = universe & mem if restricted else mem
        K = len(members)
        if not (min_size <= K <= max_size):
            continue
        hits = q & members
        k = len(hits)
        if k < min_overlap:
            continue
        support = [s.n_support(g) for g in hits]
        rows.append(
            {
                "id": s.id,
                "name": s.name,
                "collection": s.collection,
                "sources": ",".join(s.sources),
                "set_size": K,
                "overlap": k,
                "expected": round(n * K / N, 3) if N else 0.0,
                "fold": round((k / n) / (K / N), 3) if n and K else 0.0,
                "p": hypergeom_sf(k, N, K, n),
                "genes": sorted(hits),
                # fraction of hit genes supported by >= 2 source databases
                "consensus": round(sum(1 for x in support if x >= 2) / k, 3) if s.collection == "pathway" else None,
            }
        )
    for r, f in zip(rows, bh([r["p"] for r in rows])):
        r["fdr"] = f
    rows.sort(key=lambda r: (r["p"], -r["overlap"], r["id"]))
    return rows


def pair_ora(
    query: Iterable[str],
    sets: Iterable[GeneSet],
    *,
    ppi: PPINetwork | dict | None = None,
    collections: Iterable[str] | None = ("pathway",),
    min_pairs: int = 5,
    min_overlap: int = 1,
) -> list[dict]:
    """Gene-pair enrichment over curated relations (+ optional PPI edges within each set).

    ``ppi`` is a PPINetwork or a ready adjacency dict (gene -> neighbours), e.g.
    one confidence tier of a served Library.
    """
    lib = _select(sets, collections)
    q = set(query)
    set_pairs: dict[str, set] = {}
    adj = ppi.neighbours() if isinstance(ppi, PPINetwork) else (ppi or {})
    for s in lib:
        pairs = {canon(a, b) for a, b in s.pairs if a != b}
        if adj:
            genes = _members(s)
            pairs |= {canon(g, h) for g in genes for h in adj.get(g, ()) if h in genes}
        if len(pairs) >= min_pairs:
            set_pairs[s.id] = pairs
    universe = set().union(*set_pairs.values()) if set_pairs else set()
    qpairs = {e for e in universe if e[0] in q and e[1] in q}
    N, n = len(universe), len(qpairs)
    rows = []
    by_id = {s.id: s for s in lib}
    for sid, pairs in set_pairs.items():
        hit = pairs & qpairs
        if len(hit) < min_overlap:
            continue
        s = by_id[sid]
        rows.append(
            {
                "id": sid,
                "name": s.name,
                "collection": s.collection,
                "sources": ",".join(s.sources),
                "set_pairs": len(pairs),
                "overlap_pairs": len(hit),
                "expected": round(n * len(pairs) / N, 3) if N else 0.0,
                "p": hypergeom_sf(len(hit), N, len(pairs), n),
                "pairs": sorted(f"{a}-{b}" for a, b in hit),
            }
        )
    for r, f in zip(rows, bh([r["p"] for r in rows])):
        r["fdr"] = f
    rows.sort(key=lambda r: r["p"])
    return rows


# --------------------------------------------------------------------------- #
# Preranked GSEA
# --------------------------------------------------------------------------- #
def _es_batch(pos: np.ndarray, w: np.ndarray, n_genes: int) -> tuple[np.ndarray, np.ndarray]:
    """Enrichment scores for a batch of hit-position rows (sorted ascending).

    pos: (B, k) 0-based ranks of the hits; w: (B, k) hit weights.
    Returns (ES, argpeak) where argpeak indexes the hit at the extremum.
    """
    k = pos.shape[1]
    nr = w.sum(axis=1, keepdims=True)
    nr[nr == 0] = 1.0
    cw = np.cumsum(w, axis=1) / nr
    miss_before = (pos - np.arange(k)) / max(n_genes - k, 1)
    top = cw - miss_before  # running sum right after each hit
    bottom = (cw - w / nr) - miss_before  # just before each hit
    imax, imin = top.argmax(axis=1), bottom.argmin(axis=1)
    rows = np.arange(pos.shape[0])
    vmax, vmin = top[rows, imax], bottom[rows, imin]
    es = np.where(np.abs(vmax) >= np.abs(vmin), vmax, vmin)
    return es, np.where(np.abs(vmax) >= np.abs(vmin), imax, imin)


def gsea(
    ranking: dict[str, float],
    sets: Iterable[GeneSet],
    *,
    collections: Iterable[str] | None = None,
    min_size: int = 15,
    max_size: int = 500,
    nperm: int = 1000,
    weight: float = 1.0,
    seed: int = 0,
) -> list[dict]:
    genes = sorted(ranking, key=lambda g: -ranking[g])
    scores = np.array([ranking[g] for g in genes], dtype=float)
    absw = np.abs(scores) ** weight
    rank_of = {g: i for i, g in enumerate(genes)}
    N = len(genes)
    lib = []
    for s in _select(sets, collections):
        idx = np.array(sorted(rank_of[g] for g in _members(s) if g in rank_of), dtype=int)
        if min_size <= len(idx) <= max_size:
            lib.append((s, idx))
    if not lib:
        return []

    # shared null: each row is a random ordering prefix -> first k columns are a uniform k-subset
    rng = np.random.default_rng(seed)
    kmax = max(len(i) for _, i in lib)
    perm = np.empty((nperm, kmax), dtype=int)
    for start in range(0, nperm, 64):
        stop = min(start + 64, nperm)
        keys = rng.random((stop - start, N))
        part = np.argpartition(keys, kmax - 1, axis=1)[:, :kmax] if kmax < N else np.argsort(keys, axis=1)
        order = np.take_along_axis(keys, part, axis=1).argsort(axis=1)
        perm[start:stop] = np.take_along_axis(part, order, axis=1)
    null_by_k: dict[int, np.ndarray] = {}

    rows, null_nes = [], []
    for s, idx in lib:
        k = len(idx)
        es, peak = _es_batch(idx[None, :], absw[idx][None, :], N)
        es, peak = float(es[0]), int(peak[0])
        if k not in null_by_k:
            pos = np.sort(perm[:, :k], axis=1)
            null_by_k[k] = _es_batch(pos, absw[pos], N)[0]
        null = null_by_k[k]
        pos_null, neg_null = null[null >= 0], null[null < 0]
        if es >= 0:
            p = (np.sum(pos_null >= es) + 1) / (len(pos_null) + 1)
            denom = pos_null.mean() if len(pos_null) else 1.0
            lead = idx[: peak + 1]
        else:
            p = (np.sum(neg_null <= es) + 1) / (len(neg_null) + 1)
            denom = abs(neg_null.mean()) if len(neg_null) else 1.0
            lead = idx[peak:]
        nes = es / denom if denom else 0.0
        nn = np.where(null >= 0, null / (pos_null.mean() if len(pos_null) else 1), null / (abs(neg_null.mean()) if len(neg_null) else 1))
        null_nes.append(nn)
        rows.append(
            {
                "id": s.id,
                "name": s.name,
                "collection": s.collection,
                "sources": ",".join(s.sources),
                "set_size": k,
                "es": round(es, 4),
                "nes": round(float(nes), 4),
                "p": float(p),
                "leading_edge": [genes[i] for i in lead],
            }
        )
    # GSEA FDR (Subramanian 2005): compare NES with the pooled normalised null, per sign
    pool = np.concatenate(null_nes)
    pool_pos, pool_neg = np.sort(pool[pool >= 0]), np.sort(pool[pool < 0])
    obs = np.array([r["nes"] for r in rows])
    obs_pos, obs_neg = np.sort(obs[obs >= 0]), np.sort(obs[obs < 0])
    for r in rows:
        x = r["nes"]
        if x >= 0:
            f_null = (len(pool_pos) - np.searchsorted(pool_pos, x, "left")) / max(len(pool_pos), 1)
            f_obs = (len(obs_pos) - np.searchsorted(obs_pos, x, "left")) / max(len(obs_pos), 1)
        else:
            f_null = np.searchsorted(pool_neg, x, "right") / max(len(pool_neg), 1)
            f_obs = np.searchsorted(obs_neg, x, "right") / max(len(obs_neg), 1)
        r["fdr"] = float(min(1.0, f_null / f_obs)) if f_obs else 1.0
    for r, f in zip(rows, bh([r["p"] for r in rows])):
        r["padj"] = f
    rows.sort(key=lambda r: (r["p"], -abs(r["nes"])))
    return rows


# --------------------------------------------------------------------------- #
# Redundancy-aware clustering of results
# --------------------------------------------------------------------------- #
def cluster(results: list[dict], sets, *, alpha: float = 0.05, min_jaccard: float = 0.5,
            fdr_key: str = "fdr", max_hits: int = 500) -> list[dict]:
    """Annotate significant results with ``theme`` / ``theme_leader`` (in place, returns them)."""
    by_id = sets if isinstance(sets, dict) else {s.id: s for s in sets}
    leaders: list[tuple[dict, frozenset]] = []
    n_done = 0
    for r in results:
        if r.get(fdr_key, 1.0) > alpha or n_done >= max_hits:
            continue
        n_done += 1
        s = by_id.get(r["id"])
        a = _members(s) if s is not None else frozenset()
        for lead, b in leaders:
            t = by_id.get(lead["id"])
            if s is None or t is None:
                continue
            linked = t.id in s.links or s.id in t.links
            # Jaccard >= j needs min/max size >= j: skip the intersection when impossible
            if not linked and (not a or not b or min(len(a), len(b)) < min_jaccard * max(len(a), len(b))):
                continue
            if linked or len(a & b) / len(a | b) >= min_jaccard:
                r["theme"], r["theme_leader"] = lead["theme"], lead["name"]
                break
        else:
            r["theme"], r["theme_leader"] = len(leaders) + 1, r["name"]
            leaders.append((r, a))
    return results


__all__ = ["ora", "pair_ora", "gsea", "cluster", "hypergeom_sf", "bh"]
