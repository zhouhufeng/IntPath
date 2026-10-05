"""Core data structures shared by the builders and the enrichment engine."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

# Unified gene-relationship vocabulary (Table 1 of the 2012 paper) plus the
# relationship classes IntPathV2 adds when PPIs and complexes are merged in.
RELATIONS = {
    "ECrel": "enzyme-enzyme relation (successive reaction steps)",
    "PPrel": "protein-protein interaction (binding, modification, shared control)",
    "GErel": "gene expression interaction (transcription factor -> target)",
    "GPrel": "proteins in the same complex/group, not necessarily direct",
    "PPI": "physical protein-protein interaction from an interaction database (IntPathV2)",
}


@dataclass
class SourcePathway:
    """A pathway as extracted and ID-normalised from one source database."""

    source: str
    name: str
    source_id: str = ""
    genes: set[str] = field(default_factory=set)
    # (geneA, geneB) -> set of unified relation types, as given by the source
    pairs: dict[tuple[str, str], set[str]] = field(default_factory=dict)
    # source ids of parent pathways in the source's own hierarchy (Reactome)
    parents: set[str] = field(default_factory=set)

    @property
    def key(self) -> tuple[str, str]:
        return (self.source, self.name)


@dataclass
class GeneSet:
    """An IntPath gene set: an integrated pathway, a GO term, or any library entry.

    ``genes`` maps gene symbol -> set of sources supporting the membership, so
    the "full unification" provenance survives into enrichment results.
    """

    id: str
    name: str
    collection: str  # "pathway", "GO:BP", "GO:MF", "GO:CC", ...
    genes: dict[str, set[str]] = field(default_factory=dict)
    members: list[tuple[str, str]] = field(default_factory=list)  # (source, original name or id)
    pairs: dict[tuple[str, str], dict[str, set[str]]] = field(default_factory=dict)  # pair -> {"rel":..., "src":...}
    links: list[str] = field(default_factory=list)  # related set ids (e.g. pathway <-> GO cross-links)

    @property
    def sources(self) -> list[str]:
        return sorted({s for s, _ in self.members})

    @property
    def size(self) -> int:
        return len(self.genes)

    @property
    def members_set(self) -> frozenset[str]:
        return frozenset(self.genes)

    def n_support(self, gene: str) -> int:
        return len(self.genes.get(gene, ())) or 1


def stable_id(prefix: str, members: list[tuple[str, str]]) -> str:
    """Deterministic ID from the sorted member list (stable across rebuilds)."""
    h = hashlib.sha1("|".join(f"{s}:{n}" for s, n in sorted(members)).encode()).hexdigest()[:10].upper()
    return f"{prefix}{h}"
