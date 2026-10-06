"""Pathway maps: graphs of integrated pathways and their standard exports.

An integrated pathway merges several source diagrams, so it has no hand-drawn
layout; its map is the graph of its genes and unified gene relationships,
laid out automatically in the browser (Cytoscape.js). Exports:

* Cytoscape.js JSON (``elements``), for the web view and Cytoscape desktop;
* SBML Level 3 with the Qualitative Models package (SBML-qual): genes as
  qualitative species, directed relations as transitions. This is the SBML
  flavour for gene-relationship networks; core SBML describes quantitative
  reaction models;
* SBGN-ML, Activity Flow language: genes as biological activities, directed
  relations as influence arcs, for drawing tools such as Newt or SBGN-ED.
"""

from __future__ import annotations

import sqlite3
from xml.sax.saxutils import escape, quoteattr

# Unified relation -> (directed?, SBML-qual input sign, SBGN Activity Flow arc class)
RELATION_STYLE = {
    "ECrel": (True, "unknown", "unknown influence"),  # successive enzymatic steps
    "PPrel": (True, "unknown", "unknown influence"),  # binding / modification / control
    "GErel": (True, "unknown", "unknown influence"),  # transcription factor -> target
    "GPrel": (False, None, None),  # same complex: undirected, not a transition
    "PPI": (False, None, None),  # physical interaction from interaction databases
}
MAX_NODES = 600
_TIERS = {"high": ("high",), "medium": ("high", "medium"), "low": ("high", "medium", "low")}


def set_graph(con: sqlite3.Connection, set_id: str, *, ppi_tier: str | None = None,
              highlight: set[str] | None = None, max_nodes: int = MAX_NODES,
              adjacency: dict[str, frozenset[str]] | None = None) -> dict:
    """Nodes and edges of one gene set from a release database.

    Curated edges are the set's unified gene pairs; with ``ppi_tier``, merged PPI
    edges among the set's genes are added (kind "ppi"). Very large sets keep the
    ``max_nodes`` best-connected genes so the browser stays responsive.
    """
    con.row_factory = sqlite3.Row
    head = con.execute("SELECT set_id, name, collection FROM gset WHERE set_id = ?", (set_id,)).fetchone()
    if head is None:
        raise KeyError(set_id)
    genes = {r["symbol"]: r["sources"].split(",") for r in
             con.execute("SELECT symbol, sources FROM set_gene WHERE set_id = ?", (set_id,))}
    edges = []
    for r in con.execute("SELECT gene_a, gene_b, relations, sources FROM set_pair WHERE set_id = ?", (set_id,)):
        if r["gene_a"] in genes and r["gene_b"] in genes and r["gene_a"] != r["gene_b"]:
            edges.append({"source": r["gene_a"], "target": r["gene_b"], "relations": r["relations"].split(","),
                          "sources": r["sources"].split(","), "kind": "curated"})
    if ppi_tier and adjacency is not None:  # in-memory PPI network of the served library (fast)
        curated = {frozenset((e["source"], e["target"])) for e in edges}
        for g in genes:
            for h in adjacency.get(g, ()):
                if g < h and h in genes and frozenset((g, h)) not in curated:
                    edges.append({"source": g, "target": h, "relations": ["PPI"], "sources": ["PPI"],
                                  "kind": "ppi", "tier": ppi_tier})
    elif ppi_tier:
        allowed = _TIERS[ppi_tier]
        curated = {frozenset((e["source"], e["target"])) for e in edges}
        names = list(genes)
        for i in range(0, len(names), 400):  # SQLite parameter limit; gene_a < gene_b, so gene_a suffices
            chunk = names[i:i + 400]
            q = f"SELECT gene_a, gene_b, sources, tier FROM ppi WHERE gene_a IN ({','.join('?' * len(chunk))})"
            for r in con.execute(q, chunk):
                pair = frozenset((r["gene_a"], r["gene_b"]))
                if r["gene_b"] in genes and r["tier"] in allowed and pair not in curated:
                    edges.append({"source": r["gene_a"], "target": r["gene_b"], "relations": ["PPI"],
                                  "sources": r["sources"].split(","), "kind": "ppi", "tier": r["tier"]})
    degree: dict[str, int] = {}
    for e in edges:
        degree[e["source"]] = degree.get(e["source"], 0) + 1
        degree[e["target"]] = degree.get(e["target"], 0) + 1
    shown = sorted(genes, key=lambda g: (-degree.get(g, 0), g))[:max_nodes]
    keep = set(shown)
    hl = highlight or set()
    return {
        "id": head["set_id"], "name": head["name"], "collection": head["collection"],
        "n_genes": len(genes), "truncated": len(genes) > max_nodes,
        "nodes": [{"id": g, "sources": genes[g], "degree": degree.get(g, 0), "hit": g in hl} for g in shown],
        "edges": [e for e in edges if e["source"] in keep and e["target"] in keep],
    }


def _directed(e: dict) -> bool:
    rel = e["relations"][0] if e["relations"] else "PPrel"
    return e["kind"] == "curated" and RELATION_STYLE.get(rel, (True,))[0]


def to_cytoscape(graph: dict) -> dict:
    els = [{"data": {"id": n["id"], "label": n["id"], "sources": n["sources"], "n_sources": len(n["sources"]),
                     "degree": n["degree"], "hit": n["hit"]}} for n in graph["nodes"]]
    for i, e in enumerate(graph["edges"]):
        els.append({"data": {"id": f"e{i}", "source": e["source"], "target": e["target"],
                             "relation": e["relations"][0] if e["relations"] else "PPrel",
                             "relations": e["relations"], "sources": e["sources"], "kind": e["kind"],
                             "directed": _directed(e)}})
    return {"data": {"id": graph["id"], "name": graph["name"], "collection": graph["collection"]},
            "elements": els, "n_genes": graph["n_genes"], "truncated": graph["truncated"]}


def _sid(text: str) -> str:
    """A valid SBML/SBGN SId from a gene symbol or set id."""
    return "g_" + "".join(c if c.isalnum() or c == "_" else "_" for c in text)


def to_sbml_qual(graph: dict) -> str:
    """SBML Level 3 Version 1 core + qual Version 1 (qualitative network)."""
    out = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1"',
        '      xmlns:qual="http://www.sbml.org/sbml/level3/version1/qual/version1" qual:required="true">',
        f'  <model id={quoteattr(_sid(graph["id"]))} name={quoteattr(graph["name"])}>',
        '    <notes><body xmlns="http://www.w3.org/1999/xhtml"><p>IntPathV2 gene set '
        f'{escape(graph["id"])} ({escape(graph["collection"])}): genes are qualitative species; directed '
        'unified gene relationships (ECrel, PPrel, GErel) are transitions. Undirected relations (GPrel, PPI) '
        'are not transitions and are omitted.</p></body></notes>',
        '    <listOfCompartments>',
        '      <compartment id="cell" constant="true"/>',
        '    </listOfCompartments>',
        '    <qual:listOfQualitativeSpecies>',
    ]
    for n in graph["nodes"]:
        out.append(f'      <qual:qualitativeSpecies qual:id={quoteattr(_sid(n["id"]))} qual:name={quoteattr(n["id"])}'
                   ' qual:compartment="cell" qual:constant="false" qual:maxLevel="1"/>')
    out.append('    </qual:listOfQualitativeSpecies>')
    transitions = []
    for e in graph["edges"]:
        if not _directed(e):
            continue
        k = len(transitions) + 1
        rel = e["relations"][0]
        a, b = _sid(e["source"]), _sid(e["target"])
        label = f"{rel}: {e['source']} -> {e['target']}"
        transitions += [
            f'      <qual:transition qual:id="t{k}" qual:name={quoteattr(label)}>',
            '        <qual:listOfInputs>',
            f'          <qual:input qual:id="t{k}_in" qual:qualitativeSpecies="{a}" qual:transitionEffect="none"'
            f' qual:sign="{RELATION_STYLE.get(rel, (True, "unknown"))[1]}"/>',
            '        </qual:listOfInputs>',
            '        <qual:listOfOutputs>',
            f'          <qual:output qual:id="t{k}_out" qual:qualitativeSpecies="{b}" qual:transitionEffect="assignmentLevel"/>',
            '        </qual:listOfOutputs>',
            '        <qual:listOfFunctionTerms>',
            '          <qual:defaultTerm qual:resultLevel="0"/>',
            '          <qual:functionTerm qual:resultLevel="1">',
            '            <math xmlns="http://www.w3.org/1998/Math/MathML"><apply><geq/>'
            f'<ci>{a}</ci><cn type="integer">1</cn></apply></math>',
            '          </qual:functionTerm>',
            '        </qual:listOfFunctionTerms>',
            '      </qual:transition>',
        ]
    if transitions:
        out += ['    <qual:listOfTransitions>', *transitions, '    </qual:listOfTransitions>']
    out += ['  </model>', '</sbml>']
    return "\n".join(out) + "\n"


def to_sbgn_af(graph: dict, positions: dict[str, tuple[float, float]] | None = None) -> str:
    """SBGN-ML 0.3, Activity Flow. Genes without browser positions are placed on a grid."""
    pos = dict(positions or {})
    cols = max(1, int(len(graph["nodes"]) ** 0.5))
    for i, n in enumerate(graph["nodes"]):
        pos.setdefault(n["id"], ((i % cols) * 140.0, (i // cols) * 70.0))
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<sbgn xmlns="http://sbgn.org/libsbgn/0.3">',
           f'  <map id={quoteattr(_sid(graph["id"]))} language="activity flow">',
           f'    <notes><body xmlns="http://www.w3.org/1999/xhtml"><p>{escape(graph["name"])} (IntPathV2 '
           f'{escape(graph["id"])})</p></body></notes>']
    for n in graph["nodes"]:
        x, y = pos[n["id"]]
        out += [f'    <glyph id={quoteattr(_sid(n["id"]))} class="biological activity">',
                f'      <label text={quoteattr(n["id"])}/>',
                f'      <bbox x="{x:.1f}" y="{y:.1f}" w="100" h="40"/>',
                '    </glyph>']
    k = 0
    for e in graph["edges"]:
        if not _directed(e):
            continue
        k += 1
        arc = RELATION_STYLE.get(e["relations"][0], (True, None, "unknown influence"))[2]
        (x1, y1), (x2, y2) = pos[e["source"]], pos[e["target"]]
        out += [f'    <arc id="a{k}" class="{arc}" source={quoteattr(_sid(e["source"]))} target={quoteattr(_sid(e["target"]))}>',
                f'      <start x="{x1 + 50:.1f}" y="{y1 + 20:.1f}"/>',
                f'      <end x="{x2 + 50:.1f}" y="{y2 + 20:.1f}"/>',
                '    </arc>']
    out += ['  </map>', '</sbgn>']
    return "\n".join(out) + "\n"
