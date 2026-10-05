import random
from pathlib import Path

import pytest

from intpath.names import align, find_related_pairs, group_related, lcs_length, numbered_entity_mismatch

DATA = Path(__file__).resolve().parents[1] / "Data" / "sapiens"


def lcs_dp(a, b):
    prev = [0] * (len(b) + 1)
    for ca in a:
        cur = [0]
        for j, cb in enumerate(b):
            cur.append(prev[j] + 1 if ca == cb else max(prev[j + 1], cur[j]))
        prev = cur
    return prev[-1]


def test_lcs_matches_dynamic_programming():
    rng = random.Random(1)
    for _ in range(300):
        a = "".join(rng.choice("abcde ") for _ in range(rng.randint(0, 40)))
        b = "".join(rng.choice("abcde ") for _ in range(rng.randint(0, 40)))
        assert lcs_length(a, b) == lcs_dp(a, b)


@pytest.mark.parametrize(
    "a,b,score,ratio",
    [  # values printed by the 2012/2021 Java run (RelPthNamsGEN)
        ("Focal adhesion", "Focal Adhesion", 14, 1.0),
        ("Non-homologous end-joining", "Non-homologous end joining", 25, 0.9615384615384616),
        ("TGF-beta signaling pathway", "TGF Beta Signaling Pathway", 25, 0.9615384615384616),
        ("Proteasome", "Proteasome Degradation", 10, 0.625),
    ],
)
def test_alignment_reproduces_legacy_values(a, b, score, ratio):
    aln = align(a, b)
    assert aln.score == score and aln.ratio == pytest.approx(ratio)


@pytest.mark.parametrize(
    "a,b,expected",
    [
        ("Signaling by FGFR1", "Signaling by FGFR2", True),
        ("IL-1 signaling", "IL-17 signaling", True),
        ("RHOC GTPase cycle", "RHOG GTPase cycle", True),
        ("Familial hyperlipidemia type 1", "Familial hyperlipidemia type 2", True),
        ("Type I diabetes mellitus", "Type II diabetes mellitus", True),
        ("Type II diabetes mellitus", "Type 2 diabetes", False),
        ("heme biosynthesis I", "heme biosynthesis II", False),  # BioCyc variants: merged since 2012
        ("TGF-beta signaling pathway", "TGF Beta Signaling Pathway", False),
        ("PI3K-Akt signaling pathway", "PI3K-AKT Signaling Pathway", False),
    ],
)
def test_entity_guard(a, b, expected):
    assert numbered_entity_mismatch(a, b) is expected


V2_ORGANISMS = ["sapiens", "musculus", "cerevisiae", "tuberculosis"]


@pytest.mark.parametrize("organism", V2_ORGANISMS)
def test_reproduces_legacy_v2_release(organism):
    """The port reproduces each V2 organism's related pathways and integrated pathway-gene table exactly."""
    root = DATA.parent / organism
    if not (root / "normalized").exists():
        pytest.skip("V2 data not present")
    from intpath.io import read_legacy_source
    from intpath.unify import merge_pathways

    pws = []
    for code, folder in (("K", "KEGG"), ("C", "BioCyc"), ("W", "WikiPathways")):
        pws += read_legacy_source(root / "normalized" / folder / f"{organism}{folder}NormPthGEN", None, code)
    sets, matches = merge_pathways(pws, legacy=True, organism=organism)

    legacy_pairs = set()
    for line in open(root / "integrated" / "ReltedPathNames" / "RelPthNamsGEN", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) == 6:
            legacy_pairs.add(frozenset([(f[1], f[0]), (f[3], f[2])]))
    assert {frozenset([m.a, m.b]) for m in matches} == legacy_pairs

    legacy_rows = set()
    for line in open(root / "integrated" / "IntPathData" / f"{organism}IntPathGenes", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 1:
            legacy_rows.add((f[0].strip(), f[1]))
    assert {(s.name, g) for s in sets for g in s.genes} == legacy_rows


def test_human_groups():
    if not (DATA / "normalized").exists():
        pytest.skip("V2 data not present")
    from intpath.io import read_legacy_source

    n = DATA / "normalized"
    names = {}
    for code, folder in (("K", "KEGG"), ("C", "BioCyc"), ("W", "WikiPathways")):
        names[code] = [p.name for p in read_legacy_source(n / folder / f"sapiens{folder}NormPthGEN", None, code)]
    groups = group_related(find_related_pairs(names))
    assert len(groups) == 57 and sum(len(g.members) for g in groups) == 136
