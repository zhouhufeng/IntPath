# IntPath

**IntPath** is an integrated pathway and gene-relationship database for model organisms and
important pathogens. It also provides gene set enrichment analysis on the integrated data.
This repository holds the public scripts, documentation and summary statistics for
**IntPath V3** (in development, to be served at **https://intpath.genohub.org**). It also
keeps the V2 code for reference.

> Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L. IntPath: an integrated pathway gene
> relationship database for model organisms and important pathogens. *BMC Systems Biology*
> 2012, 6(Suppl 2):S2. https://doi.org/10.1186/1752-0509-6-S2-S2

## What IntPath does

Pathway databases describe the same biology with different formats, gene identifiers,
relationship vocabularies and pathway names. IntPath extracts the pathways, normalises gene IDs
and relationships, finds the pathways that different databases (and each database internally)
describe under related names, and **fully unifies** them. Every gene and gene pair is kept,
with the sources that support it ("no deletion, no introduced noise").

V3 keeps this method and extends it in three ways:

| | V2 (2012–2021) | V3 |
|---|---|---|
| organisms | human, mouse, yeast, *M. tuberculosis* | the same 4 plus rat, zebrafish, fly, worm, *Arabidopsis*, *E. coli*; extensible registry |
| pathways | KEGG, WikiPathways, BioCyc | KEGG, **Reactome**, WikiPathways, BioCyc (licensed) |
| PPIs | STRING (used for distances) | **merged** STRING + BioGRID + IntAct/MITAB with evidence, overlaid on pathways |
| GO | – | **GO BP/MF/CC**, propagated, lossless merge of redundant terms, linked to pathways |
| enrichment | hypergeometric "Identify Pathways" | ORA with **source consensus**, **gene-pair (network) enrichment**, **preranked GSEA**, redundancy-aware **themes** |
| code | Java + MySQL | Python package + REST API + web UI |

The Python port reproduces all four V2 releases exactly (human, mouse, yeast,
*M. tuberculosis*): every related-pathway pair, integrated pathway name and pathway-gene row.
The first V3 human build integrates KEGG, Reactome and WikiPathways into 4,003 pathway sets
(14,849 genes, 173,645 gene pairs). It adds 10,721 GO sets and 85,998 STRING physical
interactions. See [docs/METHODS.md](docs/METHODS.md) and
[docs/LEGACY_V2.md](docs/LEGACY_V2.md).

## Quick start

```bash
pip install -e '.[web,test]'
intpath organisms                         # supported organisms
intpath build sapiens                     # download sources and build a release in Data/v3/release/sapiens
intpath build sapiens --public            # only openly redistributable sources (no KEGG/BioCyc)
intpath ora   Data/v3/release/sapiens genes.txt --top 20
intpath pairs Data/v3/release/sapiens genes.txt --with-ppi
intpath gsea  Data/v3/release/sapiens ranking.tsv --collections pathway,GO:BP
intpath serve                             # web UI + REST API at http://127.0.0.1:8000 (docs at /docs)
pytest
```

## Repository layout

```
src/intpath/   names.py (name alignment + union-find)  unify.py (full unification)
               mapping.py (gene IDs)  organisms.py  sources.py (KEGG/Reactome/WikiPathways/BioCyc)
               ppi.py  go.py  enrich.py (ORA, pairs, GSEA, themes)  build.py  cli.py
web/           FastAPI service and single-page UI for intpath.genohub.org
docs/          METHODS, DATA_FORMATS, ROADMAP, DEPLOY, LEGACY_V2
legacy/java/   IntPath V2 Java sources: human pipeline, tools/ (alignment, DB writers), organisms/<org>/
stats/         summary statistics of builds
```

## Data licensing

The code is public. Built data inherit their sources' licences. WikiPathways, Reactome, GO,
NCBI Gene, HGNC, STRING, BioGRID and IntAct are open (CC0, CC-BY or MIT). **KEGG** and
**BioCyc** content may not be redistributed without a licence, so public downloads are built
with `--public`.
