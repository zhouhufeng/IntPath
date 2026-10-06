"""Source pathway diagrams (KEGG KGML, WikiPathways GPML) as positioned drawings.

An IntPathV2 map should look like the pathway diagram a biologist knows, not a
hairball. KEGG and WikiPathways publish their hand-drawn layouts: every gene box,
compound, compartment and arrow with coordinates and an interaction type. This
module converts both into one drawing format the web UI renders at the original
positions (Cytoscape.js, preset layout), with the user's genes highlighted.

Drawing format (JSON):
    {"source", "id", "name", "width", "height",
     "nodes": [{"id", "kind", "label", "x", "y", "w", "h", "shape", "genes", "parent", "fill", "z"}],
     "edges": [{"id", "source", "target", "arrow", "tail", "line", "label", "kind", "color"}]}
kind (nodes): gene | compound | map | group | shape | label | title | point
arrow/tail: Cytoscape arrow shapes (triangle, tee, circle, triangle-cross, none)
line: solid | dashed | dotted

Reactome diagrams come from Reactome's published layout files
(download/current/diagram/<stId>.json + .graph.json): positioned entities in
compartments, reactions with inputs, outputs, catalysts and regulators. A
sub-pathway without its own diagram is shown in its nearest ancestor's diagram,
with its reactions in focus.

Licences: KEGG drawings are KEGG's (licensed tier only). WikiPathways is CC0;
Reactome is CC BY 4.0 (attribution shown in the UI).
"""

from __future__ import annotations

import json
import re
import sqlite3
import xml.etree.ElementTree as ET
import zipfile
import zlib
from pathlib import Path
from typing import Callable

# --------------------------------------------------------------------------- #
# KEGG KGML
# --------------------------------------------------------------------------- #
# relation subtype -> (arrow, line, label); the first subtype with an arrow wins
_KGML_SUBTYPE = {
    "activation": ("triangle", "solid", ""),
    "expression": ("triangle", "solid", "e"),
    "inhibition": ("tee", "solid", ""),
    "repression": ("tee", "solid", "e"),
    "indirect effect": ("triangle", "dashed", ""),
    "state change": ("triangle", "dotted", ""),
    "binding/association": ("none", "solid", ""),
    "dissociation": ("triangle-cross", "solid", ""),
    "missing interaction": ("triangle", "dotted", ""),
    "phosphorylation": (None, None, "+p"),
    "dephosphorylation": (None, None, "-p"),
    "glycosylation": (None, None, "+g"),
    "ubiquitination": (None, None, "+u"),
    "methylation": (None, None, "+m"),
}


def kgml_diagram(xml: str | bytes, sym_of: Callable[[str], str | None]) -> dict | None:
    """KGML -> drawing. ``sym_of`` maps a KEGG gene id ("hsa:7157") to a symbol.

    Returns None for global maps drawn only with lines (e.g. map01100).
    """
    root = ET.fromstring(xml)
    nodes: dict[str, dict] = {}
    groups: dict[str, list[str]] = {}
    for e in root.iter("entry"):
        eid, etype = "k" + e.get("id"), e.get("type", "")
        g = e.find("graphics")
        if etype == "group":
            groups[eid] = ["k" + c.get("id") for c in e.iter("component")]
            nodes[eid] = {"id": eid, "kind": "group", "label": "", "genes": [], "z": 0}
            continue
        if g is None or g.get("type") == "line" or g.get("x") is None:
            continue
        label = (g.get("name") or e.get("name", "")).split(",")[0].strip().rstrip(".")
        kind = {"gene": "gene", "ortholog": "gene", "enzyme": "gene", "compound": "compound",
                "map": "map", "group": "group"}.get(etype, "gene")
        if kind == "map" and label.startswith("TITLE:"):
            kind, label = "title", label[6:]
        genes = []
        if kind == "gene":
            genes = sorted({s for s in (sym_of(x) for x in e.get("name", "").split()) if s})
        nodes[eid] = {
            "id": eid, "kind": kind, "label": label,
            "x": float(g.get("x")), "y": float(g.get("y")),
            "w": float(g.get("width") or 46), "h": float(g.get("height") or 17),
            "shape": {"circle": "ellipse", "roundrectangle": "round-rectangle"}.get(g.get("type"), "rectangle"),
            "genes": genes, "fill": g.get("bgcolor") if kind == "gene" else None,
            "z": 2 if kind in ("gene", "compound") else 1,
        }
    if not any(n.get("x") is not None for n in nodes.values()):
        return None
    for gid, members in groups.items():
        for m in members:
            if m in nodes:
                nodes[m]["parent"] = gid
    for gid in list(groups):  # groups whose members were all dropped
        if not any(nodes.get(m, {}).get("parent") == gid for m in groups[gid]):
            nodes.pop(gid, None)

    edges: list[dict] = []

    def add(a: str, b: str, arrow: str = "triangle", line: str = "solid", label: str = "", kind: str = "relation",
            tail: str = "none") -> None:
        if a in nodes and b in nodes and a != b:
            edges.append({"id": f"e{len(edges)}", "source": a, "target": b, "arrow": arrow, "tail": tail,
                          "line": line, "label": label, "kind": kind})

    for r in root.iter("relation"):
        a, b, rtype = "k" + r.get("entry1"), "k" + r.get("entry2"), r.get("type")
        arrow, line, labels, via = None, "solid", [], None
        for st in r.iter("subtype"):
            name = st.get("name", "")
            if name == "compound":
                via = "k" + st.get("value", "")
                continue
            ar, ln, lb = _KGML_SUBTYPE.get(name, (None, None, ""))
            if ar and arrow is None:
                arrow, line = ar, ln
            if lb:
                labels.append(lb)
        if rtype == "maplink":
            add(a, b, "none", "dashed", "", "maplink")
            continue
        arrow = arrow or ("triangle" if rtype in ("ECrel", "GErel", "PPrel", "PCrel") else "none")
        if via in nodes:  # drawn through the shared compound, as on the KEGG map
            add(a, via, "none", line, "", "relation")
            add(via, b, arrow, line, " ".join(labels), "relation")
        else:
            add(a, b, arrow, line, " ".join(labels), "relation")
    for rx in root.iter("reaction"):
        enzymes = ["k" + x for x in rx.get("id", "").split()]
        rev = rx.get("type") == "reversible"
        subs = ["k" + s.get("id") for s in rx.iter("substrate")]
        prods = ["k" + p.get("id") for p in rx.iter("product")]
        enz = next((x for x in enzymes if x in nodes), None)
        for s in subs:
            for p in prods:
                if enz:
                    add(s, enz, "none", "solid", "", "reaction", "triangle" if rev else "none")
                    add(enz, p, "triangle", "solid", "", "reaction")
                else:
                    add(s, p, "triangle", "solid", "", "reaction", "triangle" if rev else "none")
    xs = [n["x"] + n["w"] / 2 for n in nodes.values() if "x" in n]
    ys = [n["y"] + n["h"] / 2 for n in nodes.values() if "y" in n]
    return {"source": "KEGG", "id": root.get("name", "").replace("path:", ""), "name": root.get("title", ""),
            "link": root.get("link", ""), "width": max(xs, default=0), "height": max(ys, default=0),
            "nodes": list(nodes.values()), "edges": edges}


# --------------------------------------------------------------------------- #
# WikiPathways GPML (2013a)
# --------------------------------------------------------------------------- #
_GPML_ARROW = {
    "Arrow": "triangle", "mim-conversion": "triangle", "mim-stimulation": "triangle",
    "mim-necessary-stimulation": "triangle", "mim-modification": "triangle", "mim-translocation": "triangle",
    "mim-transcription-translation": "triangle", "mim-gap": "triangle", "mim-catalysis": "circle",
    "mim-inhibition": "tee", "TBar": "tee", "mim-cleavage": "triangle-cross", "mim-binding": "none",
    "mim-covalent-bond": "none", "Line": "none",
}
_GENE_TYPES = {"GeneProduct", "Protein", "Rna"}


def gpml_diagram(xml: str | bytes, sym_of: Callable[[str], str | None]) -> dict:
    """GPML -> drawing. ``sym_of`` maps an Xref id or a text label to a gene symbol."""
    root = ET.fromstring(xml)
    ns = root.tag.split("}")[0] + "}" if root.tag.startswith("{") else ""
    nodes: dict[str, dict] = {}
    group_node: dict[str, str] = {}  # GroupId -> node id
    for grp in root.iter(f"{ns}Group"):
        gid = grp.get("GraphId") or "grp_" + grp.get("GroupId", "")
        group_node[grp.get("GroupId", "")] = gid
        nodes[gid] = {"id": gid, "kind": "group", "label": "", "genes": [], "z": 0,
                      "style": grp.get("Style", "Group"), "_ref": grp.get("GroupRef")}

    def place(el, kind: str, label: str) -> dict | None:
        g = el.find(f"{ns}Graphics")
        if g is None or g.get("CenterX") is None:
            return None
        return {"id": el.get("GraphId") or f"n{len(nodes)}", "kind": kind, "label": label,
                "x": float(g.get("CenterX")), "y": float(g.get("CenterY")),
                "w": float(g.get("Width") or 60), "h": float(g.get("Height") or 20),
                "shape": {"Oval": "ellipse", "RoundedRectangle": "round-rectangle", "Mitochondria": "round-rectangle",
                          "Nucleus": "ellipse", "Cell": "round-rectangle", "Organelle": "round-rectangle",
                          "Membrane": "round-rectangle"}.get(g.get("ShapeType", ""), "rectangle"),
                "genes": [], "fill": ("#" + g.get("FillColor")) if g.get("FillColor") not in (None, "Transparent") else None,
                "color": ("#" + g.get("Color")) if g.get("Color") else None, "z": int(float(g.get("ZOrder") or 0))}

    for dn in root.iter(f"{ns}DataNode"):
        dtype = dn.get("Type", "")
        kind = "gene" if dtype in _GENE_TYPES else {"Metabolite": "compound", "Pathway": "map"}.get(dtype, "gene")
        n = place(dn, kind, (dn.get("TextLabel") or "").replace("\n", " "))
        if n is None:
            continue
        if kind == "gene":
            x = dn.find(f"{ns}Xref")
            cands = [x.get("ID", "") if x is not None else "", dn.get("TextLabel", "")]
            sym = next((s for s in (sym_of(c) for c in cands if c) if s), None)
            n["genes"] = [sym] if sym else []
        if dn.get("GroupRef") in group_node:
            n["parent"] = group_node[dn.get("GroupRef")]
        nodes[n["id"]] = n
    for sh in root.iter(f"{ns}Shape"):
        n = place(sh, "shape", (sh.get("TextLabel") or "").replace("\n", " "))
        if n:
            n["z"] = min(n["z"], -1)
            nodes[n["id"]] = n
    for lb in root.iter(f"{ns}Label"):
        n = place(lb, "label", (lb.get("TextLabel") or "").replace("\n", " "))
        if n:
            nodes[n["id"]] = n
    for gid, n in list(nodes.items()):  # nested groups
        ref = n.pop("_ref", None)
        if ref in group_node and group_node[ref] != gid:
            n["parent"] = group_node[ref]
    for gid in [k for k, v in nodes.items() if v["kind"] == "group"]:
        if not any(v.get("parent") == gid for v in nodes.values()):
            del nodes[gid]  # empty or unreferenced group

    edges: list[dict] = []
    for it in root.iter(f"{ns}Interaction"):
        g = it.find(f"{ns}Graphics")
        if g is None:
            continue
        pts = g.findall(f"{ns}Point")
        if len(pts) < 2:
            continue
        first, last = pts[0], pts[-1]

        def end(p) -> str:
            ref = p.get("GraphRef")
            if ref in nodes:
                return ref
            if ref in group_node.values() and ref in nodes:
                return ref
            pid = f"p{len(nodes)}"
            nodes[pid] = {"id": pid, "kind": "point", "label": "", "x": float(p.get("X")), "y": float(p.get("Y")),
                          "w": 1, "h": 1, "genes": [], "z": 0}
            return pid

        # anchors on this interaction (other interactions may point at them)
        x1, y1, x2, y2 = (float(first.get("X")), float(first.get("Y")), float(last.get("X")), float(last.get("Y")))
        for an in g.findall(f"{ns}Anchor"):
            pos = float(an.get("Position") or 0.5)
            aid = an.get("GraphId")
            if aid:
                nodes[aid] = {"id": aid, "kind": "point", "label": "", "x": x1 + (x2 - x1) * pos,
                              "y": y1 + (y2 - y1) * pos, "w": 1, "h": 1, "genes": [], "z": 0}
        a, b = end(first), end(last)
        head = last.get("ArrowHead") or ""
        tail = first.get("ArrowHead") or ""
        if not head and tail:
            a, b, head, tail = b, a, tail, ""
        edges.append({"id": it.get("GraphId") or f"e{len(edges)}", "source": a, "target": b,
                      "arrow": _GPML_ARROW.get(head, "none"), "tail": _GPML_ARROW.get(tail, "none"),
                      "line": "dashed" if g.get("LineStyle") == "Broken" else "solid", "label": "",
                      "kind": "interaction", "color": ("#" + g.get("Color")) if g.get("Color") else None})
    board = root.find(f"{ns}Graphics")
    return {"source": "WikiPathways", "id": (re.search(r"WP\d+", root.get("Version", "")) or [""])[0],
            "name": root.get("Name", ""), "width": float(board.get("BoardWidth") or 0) if board is not None else 0,
            "height": float(board.get("BoardHeight") or 0) if board is not None else 0,
            "nodes": list(nodes.values()), "edges": [e for e in edges if e["source"] != e["target"]]}


# --------------------------------------------------------------------------- #
# Reactome diagram layouts
# --------------------------------------------------------------------------- #
# renderableClass -> (kind, shape, fill)
_REACTOME_STYLE = {
    "Protein": ("gene", "round-rectangle", "#a8e4b0"),
    "ProteinDrug": ("gene", "round-rectangle", "#a8e4b0"),
    "Gene": ("gene", "rectangle", "#c9e9c9"),
    "RNA": ("gene", "rhomboid", "#c9e9c9"),
    "Complex": ("gene", "cut-rectangle", "#a9cfe8"),
    "ComplexDrug": ("gene", "cut-rectangle", "#a9cfe8"),
    "EntitySet": ("gene", "round-rectangle", "#a9cfe8"),
    "EntitySetDrug": ("gene", "round-rectangle", "#a9cfe8"),
    "Chemical": ("chemical", "ellipse", "#c7dbe9"),
    "ChemicalDrug": ("chemical", "ellipse", "#c7dbe9"),
    "Entity": ("gene", "rectangle", "#dddddd"),
    "ProcessNode": ("map", "round-rectangle", "#e6f4d7"),
    "EncapsulatedNode": ("map", "round-rectangle", "#e6f4d7"),
}
_CONNECTOR = {"INPUT": ("in", "none"), "OUTPUT": ("out", "triangle"), "CATALYST": ("in", "circle"),
              "ACTIVATOR": ("in", "triangle"), "INHIBITOR": ("in", "tee")}


def reactome_diagram(layout: dict, graph: dict, sym_of: Callable[[str], str | None]) -> dict:
    """Reactome layout (+ graph) JSON -> drawing; ``subpathways`` maps stId -> node ids in that sub-pathway."""
    gnodes = {n["dbId"]: n for n in graph.get("nodes", [])}
    memo: dict[int, set[str]] = {}

    def genes_of(db_id: int, depth: int = 0) -> set[str]:
        if db_id in memo:
            return memo[db_id]
        memo[db_id] = set()
        n = gnodes.get(db_id)
        out: set[str] = set()
        if n is not None and depth < 30:
            for name in n.get("geneNames") or []:
                sym = sym_of(name)
                if sym:
                    out.add(sym)
                    break
            for child in n.get("children") or []:
                out |= genes_of(child, depth + 1)
        memo[db_id] = out
        return out

    nodes, edges = [], []
    for c in layout.get("compartments", []):
        pr = c["prop"]
        nodes.append({"id": f"c{c['id']}", "kind": "shape", "label": c.get("displayName", ""),
                      "x": pr["x"] + pr["width"] / 2, "y": pr["y"] + pr["height"] / 2, "w": pr["width"],
                      "h": pr["height"], "shape": "round-rectangle", "genes": [], "z": -1})
    reaction_of: dict[int, str] = {}
    for e in layout.get("edges", []):
        rid = f"r{e['id']}"
        reaction_of[e["id"]] = rid
        nodes.append({"id": rid, "kind": "reaction", "label": "", "x": e["position"]["x"], "y": e["position"]["y"],
                      "w": 8, "h": 8, "shape": "rectangle", "genes": [], "z": 3, "reactome": e.get("reactomeId")})
    for n in layout.get("nodes", []):
        kind, shape, fill = _REACTOME_STYLE.get(n.get("renderableClass", ""), ("gene", "rectangle", "#dddddd"))
        pr = n["prop"]
        nid = f"n{n['id']}"
        nodes.append({"id": nid, "kind": kind, "label": n.get("displayName", ""), "x": pr["x"] + pr["width"] / 2,
                      "y": pr["y"] + pr["height"] / 2, "w": pr["width"], "h": pr["height"], "shape": shape,
                      "fill": fill, "genes": sorted(genes_of(n.get("reactomeId"))) if kind == "gene" else [],
                      "z": 2, "stId": n.get("stId") if kind == "map" else None})
        for con in n.get("connectors", []):
            r = reaction_of.get(con.get("edgeId"))
            direction, arrow = _CONNECTOR.get(con.get("type", ""), ("in", "none"))
            if r is None:
                continue
            a, b = (nid, r) if direction == "in" else (r, nid)
            edges.append({"id": f"e{len(edges)}", "source": a, "target": b, "arrow": arrow, "tail": "none",
                          "line": "solid", "label": "", "kind": con.get("type", "").lower()})
    # sub-pathway -> nodes taking part in its reactions (to focus a sub-pathway inside an ancestor diagram)
    by_reactome = {nd["reactome"]: nd["id"] for nd in nodes if nd.get("reactome")}
    subpathways = {}
    for sp in graph.get("subpathways", []):
        rx = {by_reactome[ev] for ev in sp.get("events", []) if ev in by_reactome}
        members = set(rx)
        for e in edges:
            if e["source"] in rx or e["target"] in rx:
                members.update((e["source"], e["target"]))
        if members:
            subpathways[sp["stId"]] = sorted(members)
    return {"source": "Reactome", "id": layout.get("stableId", ""), "name": layout.get("displayName", ""),
            "width": layout.get("maxX", 0), "height": layout.get("maxY", 0), "nodes": nodes, "edges": edges,
            "subpathways": subpathways}


REACTOME_DIAGRAMS = "https://reactome.org/download/current/diagram"


def reactome_diagram_ids(raw: Path, member_ids: list[str], parents: dict[str, set[str]]) -> dict[str, str]:
    """Reactome pathway -> stId of the diagram that draws it (itself, else the nearest ancestor with one)."""
    from .sources import fetch

    cache = raw / "reactome" / "diagram"
    cache.mkdir(parents=True, exist_ok=True)
    has: dict[str, bool] = {}

    def has_diagram(st: str) -> bool:
        if st not in has:
            miss = cache / f"{st}.none"
            if miss.exists():
                has[st] = False
            else:
                try:
                    fetch(f"{REACTOME_DIAGRAMS}/{st}.json", cache / f"{st}.json", retries=1)
                    has[st] = True
                except Exception:
                    miss.touch()
                    has[st] = False
        return has[st]

    out: dict[str, str] = {}
    for st in member_ids:
        frontier, seen = [st], set()
        while frontier:
            cur = frontier.pop(0)
            if cur in seen:
                continue
            seen.add(cur)
            if has_diagram(cur):
                out[st] = cur
                break
            frontier.extend(sorted(parents.get(cur, ())))
    return out


# --------------------------------------------------------------------------- #
# Storing drawings in a release database
# --------------------------------------------------------------------------- #
SCHEMA = "CREATE TABLE IF NOT EXISTS diagram (source TEXT, source_set_id TEXT, name TEXT, data BLOB, " \
         "PRIMARY KEY (source, source_set_id))"
# pathway -> the diagram that draws it (Reactome sub-pathways are drawn in an ancestor's diagram)
SCHEMA_OF = "CREATE TABLE IF NOT EXISTS diagram_of (source TEXT, source_set_id TEXT, diagram_id TEXT, " \
            "PRIMARY KEY (source, source_set_id))"


def pack(d: dict) -> bytes:
    return zlib.compress(json.dumps(d, separators=(",", ":")).encode(), 6)


def unpack(blob: bytes) -> dict:
    return json.loads(zlib.decompress(blob))


def add_diagrams(db_path: str | Path, raw: str | Path, org, sym_of: Callable[[str], str | None]) -> dict:
    """Store KGML/GPML drawings for every KEGG / WikiPathways member pathway of a release."""
    raw = Path(raw)
    con = sqlite3.connect(db_path)
    con.execute(SCHEMA)
    members = con.execute("SELECT DISTINCT source, source_set_id FROM set_member "
                          "WHERE source IN ('KEGG', 'WikiPathways')").fetchall()
    stats = {"KEGG": 0, "WikiPathways": 0, "missing": 0}
    kgml_dir = raw / "kegg" / (org.kegg or "-") / "kgml"
    zips = sorted((raw / "wikipathways").glob("wikipathways-*-gpml-*.zip"))
    zf = zipfile.ZipFile(zips[-1]) if zips else None
    wp_member = {}
    if zf is not None:
        for name in zf.namelist():
            m = re.search(r"_(WP\d+)_", name)
            if m:
                wp_member[m.group(1)] = name
    rows = []
    for source, sid in members:
        try:
            if source == "KEGG" and (kgml_dir / f"{sid}.xml").exists():
                d = kgml_diagram((kgml_dir / f"{sid}.xml").read_bytes(), sym_of)
            elif source == "WikiPathways" and sid in wp_member:
                d = gpml_diagram(zf.read(wp_member[sid]), sym_of)
            else:
                d = None
        except ET.ParseError:
            d = None
        if d is None:
            stats["missing"] += 1
            continue
        rows.append((source, sid, d["name"], pack(d)))
        stats[source] += 1
    con.executemany("INSERT OR REPLACE INTO diagram VALUES (?,?,?,?)", rows)

    # Reactome: every member pathway -> its diagram (own or nearest ancestor's), diagrams stored once
    con.execute(SCHEMA_OF)
    reactome = [r[0] for r in con.execute("SELECT DISTINCT source_set_id FROM set_member WHERE source = 'Reactome'")]
    stats["Reactome"] = 0
    if reactome:
        from .sources import fetch

        parents: dict[str, set[str]] = {}
        rel = raw / "reactome" / "ReactomePathwaysRelation.txt"
        if rel.exists():
            for line in open(rel):
                parent, child = line.rstrip("\n").split("\t")[:2]
                parents.setdefault(child, set()).add(parent)
        of = reactome_diagram_ids(raw, reactome, parents)
        cache = raw / "reactome" / "diagram"
        for diag in sorted(set(of.values())):
            try:
                layout = json.loads((cache / f"{diag}.json").read_text())
                graph = json.loads(Path(fetch(f"{REACTOME_DIAGRAMS}/{diag}.graph.json",
                                              cache / f"{diag}.graph.json")).read_text())
                d = reactome_diagram(layout, graph, sym_of)
            except Exception:
                continue
            con.execute("INSERT OR REPLACE INTO diagram VALUES (?,?,?,?)", ("Reactome", diag, d["name"], pack(d)))
            stats["Reactome"] += 1
        con.executemany("INSERT OR REPLACE INTO diagram_of VALUES (?,?,?)",
                        [("Reactome", st, diag) for st, diag in of.items()])
    con.commit()
    con.close()
    return stats
