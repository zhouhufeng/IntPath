"""IntPath V3 enrichment engine.

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


# --------------------------------------------------------------------------- #
# Statistics helpers
# --------------------------------------------------------------------------- #
def _lchoose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(population N, K successes, n draws)."""
    lo, hi = max(k, 0, n - (N - K)), min(K, n)
    if lo > hi:
        return 0.0 if k > hi else 1.0
    base = _lchoose(N, n)
    logs = [_lchoose(K, i) + _lchoose(N - K, n - i) - base for i in range(lo, hi + 1)]
    m = max(logs)
    return min(1.0, math.exp(m) * sum(math.exp(x - m) for x in logs))


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


def _select(sets: Iterable[GeneSet], collections: Iterable[str] | None) -> list[GeneSet]:
    cols = set(collections) if collections else None
    return [s for s in sets if cols is None or s.collection in cols or s.collection.split(":")[0] in cols]


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
) -> list[dict]:
    lib = _select(sets, collections)
    universe = {g for s in lib for g in s.genes}
    if background is not None:
        universe &= set(background)
    q = set(query) & universe
    N, n = len(universe), len(q)
    rows = []
    for s in lib:
        members = universe & s.genes.keys()
        K = len(members)
        if not (min_size <= K <= max_size):
            continue
        hits = q & members
        k = len(hits)
        if k < min_overlap:
            continue
        support = [len(s.genes[g]) for g in hits]
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
    rows.sort(key=lambda r: (r["p"], -r["overlap"]))
    return rows


def pair_ora(
    query: Iterable[str],
    sets: Iterable[GeneSet],
    *,
    ppi: PPINetwork | None = None,
    collections: Iterable[str] | None = ("pathway",),
    min_pairs: int = 5,
    min_overlap: int = 1,
) -> list[dict]:
    """Gene-pair enrichment over curated relations (+ optional PPI edges within each set)."""
    lib = _select(sets, collections)
    q = set(query)
    set_pairs: dict[str, set] = {}
    adj = ppi.neighbours() if ppi is not None else {}
    for s in lib:
        pairs = {canon(a, b) for a, b in s.pairs if a != b}
        if adj:
            genes = s.genes.keys()
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
        idx = np.array(sorted(rank_of[g] for g in s.genes if g in rank_of), dtype=int)
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
def cluster(results: list[dict], sets: Iterable[GeneSet], *, alpha: float = 0.05, min_jaccard: float = 0.5,
            fdr_key: str = "fdr") -> list[dict]:
    """Annotate significant results with ``theme`` / ``theme_leader`` (in place, returns them)."""
    by_id = {s.id: s for s in sets}
    leaders: list[dict] = []
    for r in results:
        if r.get(fdr_key, 1.0) > alpha:
            continue
        s = by_id.get(r["id"])
        for lead in leaders:
            t = by_id.get(lead["id"])
            if s is None or t is None:
                continue
            a, b = s.genes.keys(), t.genes.keys()
            inter = len(a & b)
            jac = inter / (len(a) + len(b) - inter) if a or b else 0.0
            if jac >= min_jaccard or t.id in s.links or s.id in t.links:
                r["theme"], r["theme_leader"] = lead["theme"], lead["name"]
                break
        else:
            r["theme"], r["theme_leader"] = len(leaders) + 1, r["name"]
            leaders.append(r)
    return results


__all__ = ["ora", "pair_ora", "gsea", "cluster", "hypergeom_sf", "bh"]
