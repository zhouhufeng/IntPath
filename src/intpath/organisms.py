"""Organism registry: IntPathV2 is multi-organism, like old IntPath.

Old IntPath (2012) shipped H. sapiens, M. musculus, S. cerevisiae and M. tuberculosis
H37Rv. IntPathV2 keeps those and adds the common model organisms. Every builder is
parameterised by an :class:`Organism`; adding a species means adding one entry
here (or a YAML/JSON file loaded with :func:`load_custom`), not new code.

Gene identifiers are normalised to the NCBI Gene official symbol of the
species (NCBI ``gene_info``), which exists for every organism; for human the
HGNC table is layered on top for richer alias/previous-symbol resolution.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

NCBI_GENE_INFO = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/{group}/{name}.gene_info.gz"
GO_GAF = "https://current.geneontology.org/annotations/{gaf}.gaf.gz"


@dataclass(frozen=True)
class Organism:
    key: str  # short name used in paths/URLs, e.g. "sapiens" (old IntPath convention)
    name: str  # binomial
    taxid: str
    gene_info: str  # "<group>/<file stem>" under NCBI GENE_INFO
    kegg: str | None = None  # KEGG organism code
    reactome: str | None = None  # species name in Reactome mapping files
    wikipathways: str | None = None  # species token in WikiPathways GMT file names
    string_taxid: str | None = None
    go_gaf: str | None = None  # GO Consortium GAF stem
    biocyc: str | None = None  # BioCyc PGDB (subscription licence since 2017)
    hgnc: bool = False  # human: layer HGNC on top of NCBI gene_info
    # UniProt by_organism idmapping file code (e.g. "HUMAN_9606"); None -> UniProt REST stream by taxid
    uniprot: str | None = None
    in_old_intpath: bool = False
    aliases: tuple[str, ...] = field(default=())

    @property
    def gene_info_is_shared(self) -> bool:
        """True when gene_info is a multi-organism file (All_*) that must be filtered by taxid."""
        return self.gene_info.split("/")[-1].startswith("All_")

    @property
    def gene_info_url(self) -> str:
        group, name = self.gene_info.split("/")
        return NCBI_GENE_INFO.format(group=group, name=name)

    @property
    def go_gaf_url(self) -> str | None:
        return GO_GAF.format(gaf=self.go_gaf) if self.go_gaf else None


ORGANISMS: dict[str, Organism] = {
    o.key: o
    for o in [
        Organism("sapiens", "Homo sapiens", "9606", "Mammalia/Homo_sapiens", "hsa", "Homo sapiens", "Homo_sapiens",
                 "9606", "goa_human", "HUMAN", hgnc=True, in_old_intpath=True, uniprot="HUMAN_9606", aliases=("human", "hsa", "9606")),
        Organism("musculus", "Mus musculus", "10090", "Mammalia/Mus_musculus", "mmu", "Mus musculus", "Mus_musculus",
                 "10090", "mgi", "MOUSE", in_old_intpath=True, uniprot="MOUSE_10090", aliases=("mouse", "mmu", "10090")),
        Organism("cerevisiae", "Saccharomyces cerevisiae", "559292", "Fungi/Saccharomyces_cerevisiae", "sce",
                 "Saccharomyces cerevisiae", "Saccharomyces_cerevisiae", "4932", "sgd", "YEAST", in_old_intpath=True,
                 uniprot="YEAST_559292", aliases=("yeast", "sce", "4932", "559292")),
        Organism("tuberculosis", "Mycobacterium tuberculosis H37Rv", "83332",
                 "Archaea_Bacteria/All_Archaea_Bacteria", "mtu", None, "Mycobacterium_tuberculosis",
                 "83332", None, "MTBRV", in_old_intpath=True, aliases=("mtb", "mtu", "83332", "1773")),
        Organism("norvegicus", "Rattus norvegicus", "10116", "Mammalia/Rattus_norvegicus", "rno", "Rattus norvegicus",
                 "Rattus_norvegicus", "10116", "rgd", "RAT", uniprot="RAT_10116", aliases=("rat", "rno", "10116")),
        Organism("rerio", "Danio rerio", "7955", "Non-mammalian_vertebrates/Danio_rerio", "dre", "Danio rerio",
                 "Danio_rerio", "7955", "zfin", None, uniprot="DANRE_7955", aliases=("zebrafish", "dre", "7955")),
        Organism("melanogaster", "Drosophila melanogaster", "7227", "Invertebrates/Drosophila_melanogaster", "dme",
                 "Drosophila melanogaster", "Drosophila_melanogaster", "7227", "fb", "FLY", uniprot="DROME_7227", aliases=("fly", "dme", "7227")),
        Organism("elegans", "Caenorhabditis elegans", "6239", "Invertebrates/Caenorhabditis_elegans", "cel",
                 "Caenorhabditis elegans", "Caenorhabditis_elegans", "6239", "wb", "WORM", uniprot="CAEEL_6239", aliases=("worm", "cel", "6239")),
        Organism("thaliana", "Arabidopsis thaliana", "3702", "Plants/Arabidopsis_thaliana", "ath", None,
                 "Arabidopsis_thaliana", "3702", "tair", "ARA", uniprot="ARATH_3702", aliases=("arabidopsis", "ath", "3702")),
        Organism("coli", "Escherichia coli K-12 MG1655", "511145",
                 "Archaea_Bacteria/Escherichia_coli_str._K-12_substr._MG1655", "eco", None, None, "511145", "ecocyc",
                 "ECOLI", uniprot="ECOLI_83333", aliases=("ecoli", "eco", "511145", "83333")),
    ]
}


CURATED = frozenset(ORGANISMS)  # hand-curated entries (all sources); others come from intpath.catalog


def get(name: str) -> Organism:
    """Look up by key, binomial, alias, KEGG code or taxid (case-insensitive)."""
    q = name.strip().lower()
    for o in ORGANISMS.values():
        if q in (o.key, o.name.lower(), o.taxid, o.kegg or "", *o.aliases):
            return o
    raise KeyError(f"unknown organism {name!r}; known: {', '.join(ORGANISMS)}")


def load_custom(path: str | Path) -> Organism:
    """Register an extra organism from a JSON file with the Organism fields."""
    data = json.loads(Path(path).read_text())
    data["aliases"] = tuple(data.get("aliases", ()))
    o = Organism(**data)
    ORGANISMS[o.key] = o
    return o


def table() -> list[dict]:
    return [asdict(o) for o in ORGANISMS.values()]
