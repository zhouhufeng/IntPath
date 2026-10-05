"""MSigDB as a companion reference library (Docs/Plans/IntPathV2Plans.md §4.5).

MSigDB's canonical-pathway and GO collections are copies of Reactome,
WikiPathways, KEGG and GO, which IntPathV2 already integrates from the live
sources. Merging them again would count the same evidence twice, so:

* role "set":  collections IntPathV2 has no other source for (Hallmarks,
  positional, CGP, PID, TF/miRNA targets, cancer modules, HPO, oncogenic,
  immunologic, cell type, ...) become gene sets with collection
  ``msigdb:<MSigDB collection>`` (e.g. ``msigdb:H``, ``msigdb:C7:IMMUNESIGDB``);
* role "link": MSigDB copies of sources IntPathV2 integrates are not stored as
  sets; each is mapped through ``exactSource`` (R-HSA-..., WP..., hsa..., GO:...)
  to the IntPathV2 set holding that source pathway, with the Jaccard between
  the two versions (version-drift QC);
* tier "full": BioCarta, KEGG_LEGACY and KEGG_MEDICUS carry extra licence terms
  and are only built without ``--public``.

Licence: MSigDB v2022.1+ is CC BY 4.0 except the collections above. Cite
Subramanian et al. 2005 PNAS and Liberzon et al. 2015 Cell Systems.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

from .mapping import GeneMapper
from .model import GeneSet
from .sources import fetch, fetch_text

RELEASES = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"


@dataclass(frozen=True)
class Collection:
    stem: str  # file stem before ".v<version>.<Hs|Mm>.json"
    role: str  # "set" | "link"
    tier: str = "open"  # "open" | "full"


HUMAN = (
    Collection("h.all", "set"),
    Collection("c1.all", "set"),
    Collection("c2.cgp", "set"),
    Collection("c2.cp.pid", "set"),
    Collection("c2.cp.biocarta", "set", "full"),
    Collection("c2.cp.kegg_medicus", "set", "full"),
    Collection("c2.cp.kegg_legacy", "link", "full"),
    Collection("c2.cp.reactome", "link"),
    Collection("c2.cp.wikipathways", "link"),
    Collection("c3.all", "set"),
    Collection("c4.all", "set"),
    Collection("c5.go", "link"),
    Collection("c5.hpo", "set"),
    Collection("c6.all", "set"),
    Collection("c7.all", "set"),
    Collection("c8.all", "set"),
    Collection("c9.all", "set"),
)
MOUSE = (
    Collection("mh.all", "set"),
    Collection("m1.all", "set"),
    Collection("m2.cgp", "set"),
    Collection("m2.cp.biocarta", "set", "full"),
    Collection("m2.cp.reactome", "link"),
    Collection("m2.cp.wikipathways", "link"),
    Collection("m3.all", "set"),
    Collection("m5.go", "link"),
    Collection("m5.mpt", "set"),
    Collection("m7.all", "set"),
    Collection("m8.all", "set"),
)
SPECIES = {"sapiens": ("Hs", HUMAN), "musculus": ("Mm", MOUSE)}


def latest_version(species: str) -> str:
    """Newest release folder for 'Hs' or 'Mm', e.g. '2026.1'."""
    found = re.findall(rf'href="(\d{{4}}\.\d+)\.{species}/"', fetch_text(RELEASES + "/"))
    if not found:
        raise RuntimeError("could not list MSigDB releases")
    return max(found, key=lambda v: tuple(int(x) for x in v.split(".")))


@dataclass
class MSigDBResult:
    version: str
    sets: list[GeneSet]
    equivalents: list[tuple[str, str, str, str, float]]  # (name, systematic, collection, intpath set id, jaccard)
    stats: dict


def load(
    raw: Path,
    organism: str,
    mapper: GeneMapper,
    integrated: list[GeneSet],
    *,
    public: bool = True,
    version: str | None = None,
    min_size: int = 5,
    max_size: int = 2000,
) -> MSigDBResult | None:
    if organism not in SPECIES:
        return None  # MSigDB is curated for human and mouse only (projection: roadmap)
    sp, collections = SPECIES[organism]
    version = version or latest_version(sp)
    # source pathway id -> IntPathV2 set (for "link" collections)
    owner: dict[str, GeneSet] = {}
    for gs in integrated:
        for _src, sid in gs.members:
            owner[sid] = gs
    sets: list[GeneSet] = []
    equivalents = []
    stats: dict = {"version": f"{version}.{sp}", "collections": {}}
    for col in collections:
        if public and col.tier != "open":
            continue
        name = f"{col.stem}.v{version}.{sp}.json"
        try:
            path = fetch(f"{RELEASES}/{version}.{sp}/{name}", raw / "msigdb" / f"{version}.{sp}" / name)
        except Exception as exc:  # a collection missing from a release is not fatal
            stats["collections"][col.stem] = f"unavailable: {exc}"
            continue
        data = json.loads(Path(path).read_text())
        n_sets = n_links = 0
        for set_name, rec in data.items():
            genes = {g for g in (mapper.map(s) for s in rec.get("geneSymbols", [])) if g}
            if col.role == "link":
                target = owner.get(rec.get("exactSource", ""))
                if target is not None:
                    inter = len(genes & target.genes.keys())
                    union = len(genes) + target.size - inter
                    equivalents.append((set_name, rec.get("systematicName", ""), rec.get("collection", ""),
                                        target.id, round(inter / union, 4) if union else 0.0))
                    n_links += 1
                continue
            if not (min_size <= len(genes) <= max_size):
                continue
            sets.append(
                GeneSet(
                    id=rec.get("systematicName") or set_name,
                    name=set_name,
                    collection="msigdb:" + rec.get("collection", col.stem.upper()),
                    genes={g: {"MSigDB"} for g in genes},
                    members=[("MSigDB", set_name)],
                )
            )
            n_sets += 1
        stats["collections"][col.stem] = {"role": col.role, "sets": n_sets, "links": n_links, "tier": col.tier}
    return MSigDBResult(f"{version}.{sp}", sets, equivalents, stats)
