# IntPath

**IntPath** is an integrated pathway and gene-relationship database for model organisms and
important pathogens. It also provides gene set enrichment analysis on the integrated data.
This repository holds the public scripts, documentation and summary statistics of **IntPathV2**,
which is in development. It also keeps the IntPath v1 code for reference.

**Website: https://intpath.genohub.org.** Sign in with a
[Genohub Community](https://discussion.genohub.org) account. Access is approved: request
membership of the [IntPath Users](https://discussion.genohub.org/g/intpath) group, and use the
[IntPath category](https://discussion.genohub.org/c/intpath) for questions and feedback.

There are two versions:

- **IntPath v1**: the database and methods published in 2012, covering human, mouse, yeast and
  *M. tuberculosis*. Everything described in the paper below refers to this version. Its Java
  code is in `legacy/java/`, and [Docs/INTPATH_V1.md](Docs/INTPATH_V1.md) summarises it.
- **IntPathV2**: the new version, in development. It is the Python package, web service and
  documentation in this repository.

Please cite:

> Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L. IntPath: an integrated pathway gene
> relationship database for model organisms and important pathogens. *BMC Systems Biology*
> 2012, 6(Suppl 2):S2. https://doi.org/10.1186/1752-0509-6-S2-S2

## What IntPath does

Pathway databases describe the same biology with different formats, gene identifiers,
relationship vocabularies and pathway names. IntPath extracts the pathways, normalises gene IDs
and relationships, finds the pathways that describe the same biological process under related
names, and **fully unifies** them: every gene and gene pair of the merged pathways is kept ("no
deletion, no introduced noise"). Each IntPath pathway lists the original pathway names it was
merged from, with what each one means.

IntPathV2 keeps this method and extends it:

| | IntPath v1 (2012) | IntPathV2 |
|---|---|---|
| organisms | human, mouse, yeast, *M. tuberculosis* | the same 4; more organisms planned |
| pathways | merged by name alignment, checked by hand | many more pathways; merged by name alignment plus a review of all pathway names by meaning, every decision recorded with its reason |
| gene pairs | unified gene-gene relations per pathway | the same, used to draw pathway maps and networks |
| PPIs | used for distances | merged physical PPI network with confidence tiers, kept separate from the pathways |
| GO | – | GO BP/MF/CC gene sets, kept separate |
| MSigDB | – | MSigDB gene sets (human, mouse), kept separate |
| enrichment | hypergeometric "Identify Pathways" | over-representation, gene-pair (network) enrichment, preranked GSEA, results grouped into themes |
| maps | – | pathway maps drawn in the browser with your genes highlighted, plus a network view of the merged gene pairs |
| code | Java + MySQL | Python package + REST API + web UI |

The Python port reproduces the IntPath v1 database exactly for all four of its organisms: every
related pathway pair, integrated pathway name and pathway-gene row. See
[Docs/METHODS.md](Docs/METHODS.md) and [Docs/INTPATH_V1.md](Docs/INTPATH_V1.md).

## IntPath in numbers: v1 vs V2

v1 numbers are from the archived builds, reproduced exactly by the Python port. V2 numbers are
from the build of 6 October 2026.

| | | *H. sapiens* | *M. musculus* | *S. cerevisiae* | *M. tuberculosis* H37Rv |
|---|---|---|---|---|---|
| source pathways | v1 | 661 | 678 | 406 | 351 |
| | V2 | 4,198 | 2,425 | 1,090 | 143 |
| related pathway pairs (merged) | v1 | 87 | 129 | 136 | 53 |
| | V2 | 301 | 170 | 74 | 1 |
| integrated pathways | v1 | 582 | 555 | 285 | 299 |
| | V2 | 3,932 | 2,278 | 1,024 | 142 |
| genes | v1 | 7,134 | 8,013 | 1,833 | 1,146 |
| | V2 | 14,849 | 14,790 | 3,199 | 1,174 |
| pathway-gene rows | v1 | 23,873 | 24,878 | 5,285 | 4,778 |
| | V2 | 209,697 | 136,817 | 20,228 | 3,892 |
| gene pairs | v1 | 50,852 | 64,029 | 3,956 | 6,227 |
| | V2 | 193,265 | 137,639 | 5,510 | 3,153 |

- **Related pathway pairs** are pairs of pathways judged to describe the same biological
  process and merged into one IntPath pathway. In human, the 301 pairs merge 464 source
  pathways into 208 IntPath pathways.
- **Gene pairs** are the unified gene-gene relations of all pathways. **Genes** per pathway are
  used for enrichment.
- For human, V2 has 2.1× the genes, 3.8× the gene pairs and 6.8× the integrated pathways of v1.
- *M. tuberculosis* has fewer pathways and gene pairs than in v1. Most of its v1 pathways came
  from a source that is not in this build, and current pathway releases have few pathways for it.

Kept alongside the pathways (not merged into them):

| | *H. sapiens* | *M. musculus* | *S. cerevisiae* | *M. tuberculosis* H37Rv |
|---|---|---|---|---|
| physical PPI edges (merged) | 1,151,307 | 85,239 | 236,630 | 93 |
| GO gene sets | 10,724 | 11,718 | 3,451 | 1,063 |
| MSigDB gene sets | 21,902 | 4,786 | – | – |

## Quick start

```bash
pip install -e '.[web,test]'
intpath organisms                         # supported organisms
intpath build sapiens                     # download the data and build Data/intpathv2/db/sapiens
intpath ora   Data/intpathv2/db/sapiens genes.txt --top 20
intpath pairs Data/intpathv2/db/sapiens genes.txt --with-ppi
intpath gsea  Data/intpathv2/db/sapiens ranking.tsv --collections pathway,GO:BP
intpath serve                             # web UI + REST API at http://127.0.0.1:8000 (docs at /docs)
pytest
```

## Repository layout

```
src/intpath/   names.py (name alignment + union-find)  unify.py (full unification)
               review.py + data/merge_curation.tsv (recorded merge decisions)
               mapping.py (gene IDs)  organisms.py  sources.py  meaning.py (pathway descriptions)
               ppi.py  go.py  enrich.py (ORA, pairs, GSEA, themes)  diagrams.py  maps.py  build.py  cli.py
web/           FastAPI service and single-page web UI
Docs/          METHODS, DATA_FORMATS, INTPATH_V1
legacy/java/   IntPath v1 Java sources: human pipeline, tools/ (alignment, DB writers), organisms/<org>/
stats/         summary counts of the current builds
```

## Licence

The code is public. The data that IntPath builds are subject to the terms of the databases they
are built from; check those terms before redistributing a built database.
