import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")
from fastapi.testclient import TestClient  # noqa: E402

from intpath.db import write_db  # noqa: E402
from intpath.model import GeneSet  # noqa: E402
from intpath.ppi import PPINetwork  # noqa: E402


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("INTPATH_PRELOAD", "0")
    genes = [f"G{i}" for i in range(300)]
    sets = [
        GeneSet(f"IP{j}", f"pathway {j}", "pathway", {g: {"KEGG", "Reactome"} for g in genes[j * 30:(j + 1) * 30]},
                members=[("KEGG", f"p{j}"), ("Reactome", f"r{j}")])
        for j in range(10)
    ]
    sets[0].pairs = {(genes[i], genes[i + 1]): {"rel": {"PPrel"}, "src": {"KEGG"}} for i in range(29)}
    sets[0].links = ["IP1"]
    sets.append(GeneSet("M1", "HALLMARK_TEST", "msigdb:H", {g: {"MSigDB"} for g in genes[:25]}, members=[("MSigDB", "HALLMARK_TEST")]))
    net = PPINetwork()
    for i in range(30, 50):
        net.add(genes[i], genes[i + 1], "BioGRID", pmid=str(i))
        net.add(genes[i], genes[i + 1], "HuRI")
    (tmp_path / "sapiens").mkdir()
    write_db(tmp_path / "sapiens" / "intpath.sqlite", sets, meta={"organism_key": "sapiens"},
             aliases={"ALIAS0": "G0"}, ppi=net)
    (tmp_path / "sapiens" / "stats.json").write_text('{"built": "2026-10-05", "tier": "open"}')
    from web.app import create_app

    return TestClient(create_app(tmp_path))


def test_api_roundtrip(client):
    orgs = client.get("/api/organisms").json()
    assert any(o["key"] == "sapiens" and o["available"] and o["curated"] for o in orgs)
    assert client.get("/api/sapiens/stats").json()["tier"] == "open"
    r = client.post("/api/human/enrich/ora", json={"genes": [f"g{i}" for i in range(1, 20)] + ["alias0", "NOPE"]}).json()
    assert r["n_mapped"] == 20 and r["unmapped"] == ["NOPE"] and r["results"][0]["id"] in ("IP0", "M1")
    assert r["results"][0]["theme"] == 1
    hall = client.post("/api/sapiens/enrich/ora", json={"genes": [f"G{i}" for i in range(20)], "collections": ["msigdb"]}).json()
    assert {x["collection"] for x in hall["results"]} == {"msigdb:H"}
    r = client.post("/api/sapiens/enrich/pairs", json={"genes": [f"G{i}" for i in range(20)]}).json()
    assert r["results"][0]["overlap_pairs"] == 19
    r = client.post("/api/sapiens/enrich/pairs", json={"genes": [f"G{i}" for i in range(30, 60)], "with_ppi": True,
                                                       "ppi_tier": "high"}).json()
    assert r["results"][0]["id"] == "IP1" and r["results"][0]["overlap_pairs"] == 20  # G30-G31 ... G49-G50
    s = client.get("/api/sapiens/set/IP0").json()
    assert len(s["genes"]) == 30 and len(s["members"]) == 2 and s["links"][0]["id"] == "IP1"
    assert client.get("/api/sapiens/search", params={"q": "G5"}).json()[0]["id"] == "IP0"
    g = client.get("/api/sapiens/gene/G30").json()
    assert g["ppi_partners"][0]["gene"] == "G31" and g["ppi_partners"][0]["tier"] == "high"
    assert client.get("/api/sapiens/collections").json()
    assert client.get("/download/sapiens/..%2Fsecret").status_code == 404
    assert client.get("/healthz").json()["ok"]
    assert client.get("/").status_code == 200


def test_licensed_section_requires_sign_in_and_blocks_downloads(tmp_path, monkeypatch):
    monkeypatch.setenv("INTPATH_PRELOAD", "0")
    for tier in ("open", "full"):
        d = tmp_path / tier / "sapiens"
        d.mkdir(parents=True)
        sets = [GeneSet("IP0", "Apoptosis", "pathway", {f"G{i}": {"Reactome"} for i in range(20)}, members=[("Reactome", "R1")])]
        if tier == "full":
            sets.append(GeneSet("IPK", "Apoptosis - KEGG", "pathway", {f"G{i}": {"KEGG"} for i in range(20)}, members=[("KEGG", "hsa04210")]))
        write_db(d / "intpath.sqlite", sets, meta={})
        (d / "stats.json").write_text('{"tier": "%s"}' % tier)
        (d / "intpath.gmt").write_text("x")
    from web.app import create_app

    c = TestClient(create_app(tmp_path / "open", licensed_root=tmp_path / "full"))
    body = {"genes": [f"G{i}" for i in range(10)]}
    assert {r["id"] for r in c.post("/api/sapiens/enrich/ora", json=body).json()["results"]} == {"IP0"}
    assert c.post("/licensed/api/sapiens/enrich/ora", json=body).status_code == 401
    signed = {"X-IGVF-User": "alice"}
    res = c.post("/licensed/api/sapiens/enrich/ora", json=body, headers=signed).json()["results"]
    assert {r["id"] for r in res} == {"IP0", "IPK"} and any("KEGG" in r["sources"] for r in res)
    assert c.get("/download/sapiens/intpath.gmt").status_code == 200
    assert c.get("/licensed/download/sapiens/intpath.gmt", headers=signed).status_code == 404
