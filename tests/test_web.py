import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")
from fastapi.testclient import TestClient  # noqa: E402

from intpath.io import write_release  # noqa: E402
from intpath.model import GeneSet  # noqa: E402


@pytest.fixture()
def client(tmp_path):
    genes = [f"G{i}" for i in range(300)]
    sets = [
        GeneSet(f"IP{j}", f"pathway {j}", "pathway", {g: {"KEGG", "Reactome"} for g in genes[j * 30:(j + 1) * 30]},
                members=[("KEGG", f"p{j}"), ("Reactome", f"r{j}")])
        for j in range(10)
    ]
    sets[0].pairs = {(genes[i], genes[i + 1]): {"rel": {"PPrel"}, "src": {"KEGG"}} for i in range(29)}
    write_release(sets, tmp_path / "sapiens")
    from web.app import create_app

    return TestClient(create_app(tmp_path))


def test_api_roundtrip(client):
    orgs = client.get("/api/organisms").json()
    assert any(o["key"] == "sapiens" and o["available"] for o in orgs)
    r = client.post("/api/human/enrich/ora", json={"genes": [f"g{i}" for i in range(20)] + ["NOPE"]}).json()
    assert r["n_mapped"] == 20 and r["unmapped"] == ["NOPE"] and r["results"][0]["id"] == "IP0"
    r = client.post("/api/sapiens/enrich/pairs", json={"genes": [f"G{i}" for i in range(20)]}).json()
    assert r["results"][0]["overlap_pairs"] == 19
    s = client.get("/api/sapiens/set/IP0").json()
    assert len(s["genes"]) == 30 and len(s["members"]) == 2
    assert client.get("/api/sapiens/search", params={"q": "G5"}).json()[0]["id"] == "IP0"
    assert client.get("/download/sapiens/..%2Fsecret").status_code == 404
    assert client.get("/").status_code == 200
