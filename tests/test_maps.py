import sqlite3
import xml.etree.ElementTree as ET

from intpath import maps
from intpath.db import write_db
from intpath.model import GeneSet
from intpath.ppi import PPINetwork


def make_db(tmp_path):
    gs = GeneSet("IP1", "Toy signalling", "pathway", {g: {"KEGG", "Reactome"} for g in "ABCDE"},
                 members=[("KEGG", "hsa0"), ("Reactome", "R-1")])
    gs.pairs = {("A", "B"): {"rel": {"PPrel"}, "src": {"KEGG"}}, ("B", "C"): {"rel": {"GErel"}, "src": {"KEGG"}},
                ("C", "D"): {"rel": {"GPrel"}, "src": {"Reactome"}}}
    net = PPINetwork()
    net.add("D", "E", "BioGRID", pmid="1")
    net.add("D", "E", "IntAct", pmid="2")
    net.add("A", "Z", "BioGRID", pmid="3")  # partner outside the set
    path = tmp_path / "intpath.sqlite"
    write_db(path, [gs], meta={}, ppi=net)
    return sqlite3.connect(path)


def test_graph_highlights_and_ppi(tmp_path):
    con = make_db(tmp_path)
    g = maps.set_graph(con, "IP1", highlight={"A", "C"})
    assert {n["id"] for n in g["nodes"] if n["hit"]} == {"A", "C"}
    assert len(g["edges"]) == 3
    g = maps.set_graph(con, "IP1", ppi_tier="high")
    ppi = [e for e in g["edges"] if e["kind"] == "ppi"]
    assert [(e["source"], e["target"]) for e in ppi] == [("D", "E")]  # A-Z left out: Z is not in the set
    cy = maps.to_cytoscape(g)
    directed = {(e["data"]["source"], e["data"]["target"]) for e in cy["elements"] if e["data"].get("directed")}
    assert directed == {("A", "B"), ("B", "C")}  # GPrel and PPI are undirected


def test_sbml_qual_and_sbgn_are_well_formed(tmp_path):
    g = maps.set_graph(make_db(tmp_path), "IP1")
    sbml = ET.fromstring(maps.to_sbml_qual(g))
    q = "{http://www.sbml.org/sbml/level3/version1/qual/version1}"
    assert len(sbml.findall(f".//{q}qualitativeSpecies")) == 5
    assert len(sbml.findall(f".//{q}transition")) == 2
    sbgn = ET.fromstring(maps.to_sbgn_af(g))
    ns = "{http://sbgn.org/libsbgn/0.3}"
    assert len(sbgn.findall(f".//{ns}glyph")) == 5 and len(sbgn.findall(f".//{ns}arc")) == 2
