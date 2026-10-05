"""Release database: one SQLite file per organism release (``intpath.sqlite``).

SQLite (standard library, single file, read-only at serve time) replaces the
DuckDB choice in the plan (D7): the web service only needs indexed lookups,
and the statistics run on the compact in-memory library (intpath.library).

Tables
  meta             key, value                      build date, organism, source versions, tier
  gene             symbol, gene_id                 genes present in the release
  alias            alias, symbol                   unambiguous ids / previous / alias symbols -> symbol
  gset             set_id, name, collection, sources, n_genes, n_pairs
  set_member       set_id, source, source_set_id   original pathways / GO terms merged into a set
  set_gene         set_id, symbol, sources         membership with provenance
  set_pair         set_id, gene_a, gene_b, relations, sources
  set_link         set_a, set_b                    hierarchy, GO<->pathway and other cross-links
  ppi              gene_a, gene_b, sources, n_sources, n_pmids, methods, string_score, tier, pathway_sets
  msigdb_equivalent msigdb_name, systematic_name, collection, set_id, jaccard
  match_log        source_a, pathway_a, source_b, pathway_b, align_score, align_ratio, jaccard
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from .model import GeneSet
from .names import Match
from .ppi import PPINetwork

SCHEMA = """
CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE gene (symbol TEXT PRIMARY KEY, gene_id TEXT);
CREATE TABLE alias (alias TEXT PRIMARY KEY, symbol TEXT NOT NULL);
CREATE TABLE gset (set_id TEXT PRIMARY KEY, name TEXT, collection TEXT, sources TEXT, n_genes INT, n_pairs INT);
CREATE TABLE set_member (set_id TEXT, source TEXT, source_set_id TEXT);
CREATE TABLE set_gene (set_id TEXT, symbol TEXT, sources TEXT);
CREATE TABLE set_pair (set_id TEXT, gene_a TEXT, gene_b TEXT, relations TEXT, sources TEXT);
CREATE TABLE set_link (set_a TEXT, set_b TEXT);
CREATE TABLE ppi (gene_a TEXT, gene_b TEXT, sources TEXT, n_sources INT, n_pmids INT, methods TEXT,
                  string_score INT, tier TEXT, pathway_sets TEXT);
CREATE TABLE msigdb_equivalent (msigdb_name TEXT, systematic_name TEXT, collection TEXT, set_id TEXT, jaccard REAL);
CREATE TABLE match_log (source_a TEXT, pathway_a TEXT, source_b TEXT, pathway_b TEXT,
                        align_score INT, align_ratio REAL, jaccard REAL);
"""
INDEXES = """
CREATE INDEX set_gene_set ON set_gene(set_id);
CREATE INDEX set_gene_sym ON set_gene(symbol);
CREATE INDEX set_pair_set ON set_pair(set_id);
CREATE INDEX set_member_set ON set_member(set_id);
CREATE INDEX set_member_src ON set_member(source_set_id);
CREATE INDEX set_link_a ON set_link(set_a);
CREATE INDEX ppi_a ON ppi(gene_a);
CREATE INDEX ppi_b ON ppi(gene_b);
CREATE INDEX gset_collection ON gset(collection);
CREATE INDEX msig_set ON msigdb_equivalent(set_id);
"""


def write_db(
    path: str | Path,
    sets: list[GeneSet],
    *,
    meta: dict,
    aliases: dict[str, str] | None = None,
    gene_ids: dict[str, str] | None = None,
    ppi: PPINetwork | None = None,
    matches: list[Match] | None = None,
    msigdb_equivalents: list[tuple] | None = None,
) -> Path:
    path = Path(path)
    tmp = path.with_suffix(".sqlite.part")
    tmp.unlink(missing_ok=True)
    con = sqlite3.connect(tmp)
    con.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;" + SCHEMA)
    con.executemany("INSERT INTO meta VALUES (?,?)", [(k, json.dumps(v) if not isinstance(v, str) else v)
                                                       for k, v in meta.items()])
    genes = {g for s in sets for g in s.genes}
    if ppi is not None:
        genes |= {g for e in ppi.edges for g in e}
    gene_ids = gene_ids or {}
    con.executemany("INSERT INTO gene VALUES (?,?)", ((g, gene_ids.get(g, "")) for g in sorted(genes)))
    if aliases:
        con.executemany("INSERT OR REPLACE INTO alias VALUES (?,?)", aliases.items())
    con.executemany(
        "INSERT INTO gset VALUES (?,?,?,?,?,?)",
        ((s.id, s.name, s.collection, ",".join(s.sources), s.size, len(s.pairs)) for s in sets),
    )
    con.executemany("INSERT INTO set_member VALUES (?,?,?)", ((s.id, a, b) for s in sets for a, b in s.members))
    con.executemany(
        "INSERT INTO set_gene VALUES (?,?,?)",
        ((s.id, g, ",".join(sorted(src))) for s in sets for g, src in s.genes.items()),
    )
    con.executemany(
        "INSERT INTO set_pair VALUES (?,?,?,?,?)",
        ((s.id, a, b, ",".join(sorted(e["rel"])), ",".join(sorted(e["src"]))) for s in sets for (a, b), e in s.pairs.items()),
    )
    con.executemany("INSERT INTO set_link VALUES (?,?)", ((s.id, l) for s in sets for l in s.links))
    if ppi is not None:
        con.executemany(
            "INSERT INTO ppi VALUES (?,?,?,?,?,?,?,?,?)",
            (
                (a, b, ",".join(sorted(ev.sources)), len(ev.sources), len(ev.pmids), "|".join(sorted(ev.methods)),
                 ev.string_score, PPINetwork.tier(ev), ",".join(sorted(ev.pathways)))
                for (a, b), ev in ppi.edges.items()
            ),
        )
    if msigdb_equivalents:
        con.executemany("INSERT INTO msigdb_equivalent VALUES (?,?,?,?,?)", msigdb_equivalents)
    if matches:
        con.executemany(
            "INSERT INTO match_log VALUES (?,?,?,?,?,?,?)",
            ((m.a[0], m.a[1], m.b[0], m.b[1], m.score, round(m.ratio, 4), m.overlap) for m in matches),
        )
    con.executescript(INDEXES + "ANALYZE;")
    con.commit()
    con.close()
    tmp.replace(path)
    return path
