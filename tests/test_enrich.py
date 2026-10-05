import numpy as np
import pytest

from intpath.enrich import bh, cluster, gsea, hypergeom_sf, ora, pair_ora
from intpath.model import GeneSet


def test_hypergeom_against_scipy():
    stats = pytest.importorskip("scipy.stats")
    for k, N, K, n in [(5, 20000, 100, 200), (0, 100, 10, 10), (10, 500, 40, 30), (3, 50, 5, 5)]:
        assert hypergeom_sf(k, N, K, n) == pytest.approx(stats.hypergeom.sf(k - 1, N, K, n), rel=1e-9, abs=1e-300)


def test_bh_monotone():
    # sorted p*n/rank = 0.04, 0.06, 0.0533, 0.5 -> step-up minimum
    assert bh([0.01, 0.04, 0.03, 0.5]) == pytest.approx([0.04, 0.16 / 3, 0.16 / 3, 0.5])


def lib():
    genes = [f"G{i}" for i in range(1000)]
    sets = [GeneSet(f"S{j}", f"set {j}", "pathway", {g: {"A"} for g in genes[j * 50:(j + 1) * 50]}) for j in range(20)]
    sets[0].genes.update({g: {"A", "B"} for g in genes[:10]})
    sets[0].pairs = {(genes[i], genes[i + 1]): {"rel": {"PPrel"}, "src": {"A"}} for i in range(49)}
    sets[1].pairs = {(genes[50 + i], genes[51 + i]): {"rel": {"PPrel"}, "src": {"A"}} for i in range(49)}
    return genes, sets


def test_ora_finds_planted_signal_and_consensus():
    genes, sets = lib()
    rows = ora(genes[:30] + genes[500:505], sets)
    assert rows[0]["id"] == "S0" and rows[0]["fdr"] < 1e-10
    assert rows[0]["consensus"] == pytest.approx(10 / 30, abs=1e-3)


def test_pair_ora():
    genes, sets = lib()
    rows = pair_ora(genes[:20], sets)
    assert rows[0]["id"] == "S0" and rows[0]["overlap_pairs"] == 19


def test_gsea_planted_up_and_down():
    genes, sets = lib()
    rng = np.random.default_rng(0)
    ranking = {g: float(rng.normal()) for g in genes}
    for g in genes[:50]:
        ranking[g] += 3
    for g in genes[50:100]:
        ranking[g] -= 3
    rows = {r["id"]: r for r in gsea(ranking, sets, nperm=500)}
    assert rows["S0"]["nes"] > 2 and rows["S0"]["fdr"] < 0.01
    assert rows["S1"]["nes"] < -2 and rows["S1"]["fdr"] < 0.01
    assert rows["S5"]["fdr"] > 0.05
    assert set(rows["S0"]["leading_edge"]) <= set(genes[:50])


def test_cluster_groups_redundant_hits():
    genes, sets = lib()
    dup = GeneSet("S0b", "set 0 copy", "GO:BP", dict(sets[0].genes))
    rows = ora(genes[:30], sets + [dup])
    cluster(rows, sets + [dup])
    themes = {r["id"]: r.get("theme") for r in rows}
    assert themes["S0"] == themes["S0b"]
