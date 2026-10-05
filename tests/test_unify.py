from intpath.model import SourcePathway
from intpath.unify import merge_pathways


def sp(src, name, genes, pairs=()):
    p = SourcePathway(src, name, genes=set(genes))
    for a, b, r in pairs:
        p.pairs[(a, b)] = {r}
    return p


def test_full_unification_keeps_everything_with_provenance():
    pws = [
        sp("KEGG", "Focal adhesion", "ABC", [("A", "B", "PPrel")]),
        sp("WikiPathways", "Focal Adhesion", "BCD", [("A", "B", "GPrel"), ("C", "D", "PPrel")]),
        sp("Reactome", "Apoptosis", "XY"),
    ]
    sets, matches = merge_pathways(pws, legacy=False)
    assert len(matches) == 1 and len(sets) == 2
    fa = next(s for s in sets if len(s.members) == 2)
    assert set(fa.genes) == set("ABCD")
    assert fa.genes["B"] == {"KEGG", "WikiPathways"} and fa.genes["D"] == {"WikiPathways"}
    assert fa.pairs[("A", "B")] == {"rel": {"PPrel", "GPrel"}, "src": {"KEGG", "WikiPathways"}}


def test_overlap_guard_blocks_lookalike_names_without_shared_genes():
    pws = [sp("A", "Wnt signaling pathway", "ABCDE"), sp("B", "Wnt signaling pathways", "VWXYZ")]
    sets, matches = merge_pathways(pws, legacy=False, min_jaccard=0.1)
    assert not matches and len(sets) == 2
    sets, matches = merge_pathways(pws, legacy=True)
    assert len(matches) == 1
