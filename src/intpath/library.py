"""Compact in-memory gene set library for serving enrichment from a release database.

Holding every set as ``GeneSet`` (gene -> set of sources) costs ~0.7 KB per
membership; with GO, MSigDB and the merged PPIs a human release would need
several GB. Enrichment only needs membership, pathway source support (for
consensus), pathway gene pairs and links, so a :class:`LibSet` keeps exactly
those, and everything else is queried from SQLite on demand.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class LibSet:
    id: str
    name: str
    collection: str
    sources: list[str]
    members: frozenset[str]
    support: dict[str, int] | None = None  # pathway sets: gene -> number of supporting sources
    pairs: frozenset[tuple[str, str]] = frozenset()
    links: list[str] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.members)

    def n_support(self, gene: str) -> int:
        return self.support.get(gene, 1) if self.support else 1


@dataclass
class Library:
    organism: str
    sets: list[LibSet]
    by_id: dict[str, LibSet]
    by_gene: dict[str, list[str]]
    aliases: dict[str, str]  # UPPER-CASE alias/id/symbol -> symbol
    adjacency: dict[str, dict[str, frozenset[str]]]  # PPI tier -> gene -> neighbours (cumulative tiers)
    meta: dict[str, str]
    db_path: Path

    def set_index(self) -> dict[str, list[LibSet]]:
        """gene -> LibSet objects containing it (built once)."""
        idx = getattr(self, "_set_index", None)
        if idx is None:
            idx = {g: [self.by_id[sid] for sid in sids] for g, sids in self.by_gene.items()}
            self._set_index = idx
        return idx

    def universe(self, collections) -> frozenset[str]:
        """All genes of the selected collections (cached per selection)."""
        key = tuple(sorted(collections)) if collections else ()
        cache = self.__dict__.setdefault("_universe", {})
        if key not in cache:
            from .enrich import _select

            cache[key] = frozenset(g for s in _select(self.sets, collections) for g in s.members)
        return cache[key]

    def resolve(self, ids) -> tuple[list[str], list[str]]:
        found, missing = [], []
        for x in ids:
            s = self.aliases.get(str(x).strip().upper())
            (found if s else missing).append(s or x)
        return list(dict.fromkeys(found)), missing

    def connect(self) -> sqlite3.Connection:
        con = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True, check_same_thread=False)
        con.row_factory = sqlite3.Row
        return con


def load_library(db_path: str | Path, organism: str = "") -> Library:
    db_path = Path(db_path)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    intern: dict[str, str] = {}

    def I(s: str) -> str:  # noqa: E743 - share one string object per gene symbol
        return intern.setdefault(s, s)

    meta = dict(con.execute("SELECT key, value FROM meta"))
    info = {sid: (name, col, src) for sid, name, col, src in con.execute("SELECT set_id, name, collection, sources FROM gset")}
    members: dict[str, list[str]] = {}
    support: dict[str, dict[str, int]] = {}
    for sid, sym, src in con.execute("SELECT set_id, symbol, sources FROM set_gene ORDER BY set_id"):
        sym = I(sym)
        members.setdefault(sid, []).append(sym)
        if info[sid][1] == "pathway":
            n = src.count(",") + 1
            if n > 1:
                support.setdefault(sid, {})[sym] = n
    pairs: dict[str, list[tuple[str, str]]] = {}
    for sid, a, b in con.execute("SELECT set_id, gene_a, gene_b FROM set_pair"):
        pairs.setdefault(sid, []).append((I(a), I(b)) if a <= b else (I(b), I(a)))
    links: dict[str, list[str]] = {}
    for a, b in con.execute("SELECT set_a, set_b FROM set_link"):
        links.setdefault(a, []).append(b)

    sets = []
    by_gene: dict[str, list[str]] = {}
    for sid, (name, col, src) in info.items():
        mem = frozenset(members.get(sid, ()))
        s = LibSet(sid, name, col, [x for x in src.split(",") if x], mem,
                   support.get(sid) if col == "pathway" else None,
                   frozenset(pairs.get(sid, ())), links.get(sid, []))
        sets.append(s)
        for g in mem:
            by_gene.setdefault(g, []).append(sid)

    aliases = {k: I(v) for k, v in con.execute("SELECT alias, symbol FROM alias")}
    for (sym,) in con.execute("SELECT symbol FROM gene"):
        aliases[sym.upper()] = I(sym)

    order = {"high": 0, "medium": 1, "low": 2}
    raw_adj: dict[str, dict[str, set[str]]] = {"high": {}, "medium": {}, "low": {}}
    for a, b, tier in con.execute("SELECT gene_a, gene_b, tier FROM ppi"):
        a, b = I(a), I(b)
        for t, rank in order.items():  # cumulative: "medium" includes high-confidence edges
            if order.get(tier, 2) <= rank:
                raw_adj[t].setdefault(a, set()).add(b)
                raw_adj[t].setdefault(b, set()).add(a)
    adjacency = {t: {g: frozenset(n) for g, n in adj.items()} for t, adj in raw_adj.items()}
    # STRING is a separate network ("string" key); older releases merged it into ppi
    try:
        s_adj: dict[str, set[str]] = {}
        for a, b in con.execute("SELECT gene_a, gene_b FROM string_ppi"):
            a, b = I(a), I(b)
            s_adj.setdefault(a, set()).add(b)
            s_adj.setdefault(b, set()).add(a)
        adjacency["string"] = {g: frozenset(n) for g, n in s_adj.items()}
    except sqlite3.OperationalError:
        adjacency["string"] = {}
    con.close()
    return Library(organism or meta.get("organism_key", ""), sets, {s.id: s for s in sets}, by_gene, aliases,
                   adjacency, meta, db_path)
