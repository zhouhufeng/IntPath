# IntPath

**IntPath** is an integrated pathway and gene-relationship database for model organisms and
important pathogens. It also provides gene set enrichment analysis on the integrated data.
This repository holds the public scripts, documentation and summary statistics for
**IntPathV2** (in development; a first public release is live at **https://intpath.genohub.org**). It also
keeps the old IntPath code for reference.

There are two versions:

- **Old IntPath**: the database and methods published in 2012, covering human, mouse, yeast and
  *M. tuberculosis*. Everything described in the paper below refers to this version. Its Java
  code is in `legacy/java/`, and [Docs/OLD_INTPATH.md](Docs/OLD_INTPATH.md) summarises it.
- **IntPathV2**: the new version, in development. It is the Python package, web service and
  documentation in this repository.

> Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L. IntPath: an integrated pathway gene
> relationship database for model organisms and important pathogens. *BMC Systems Biology*
> 2012, 6(Suppl 2):S2. https://doi.org/10.1186/1752-0509-6-S2-S2

## What IntPath does

Pathway databases describe the same biology with different formats, gene identifiers,
relationship vocabularies and pathway names. IntPath extracts the pathways, normalises gene IDs
and relationships, finds the pathways that different databases (and each database internally)
describe under related names, and **fully unifies** them. Every gene and gene pair is kept,
with the sources that support it ("no deletion, no introduced noise").

IntPathV2 keeps this method and extends it in three ways:

| | Old IntPath (2012) | IntPathV2 |
|---|---|---|
| organisms | human, mouse, yeast, *M. tuberculosis* | the same 4 plus rat, zebrafish, fly, worm, *Arabidopsis*, *E. coli*; extensible registry |
| pathways | KEGG, WikiPathways, BioCyc | KEGG, **Reactome** (with topology and hierarchy), WikiPathways, BioCyc (licensed) |
| PPIs | STRING (used for distances) | **merged** STRING + BioGRID + IntAct + MINT + HuRI, with evidence and confidence tiers, overlaid on pathways |
| MSigDB | – | companion library (Hallmarks, C1–C9, mouse M-collections), copies of Reactome/WikiPathways/GO linked rather than double-counted |
| GO | – | **GO BP/MF/CC**, propagated, lossless merge of redundant terms, linked to pathways |
| enrichment | hypergeometric "Identify Pathways" | ORA with **source consensus**, **gene-pair (network) enrichment**, **preranked GSEA**, redundancy-aware **themes** |
| code | Java + MySQL | Python package + REST API + web UI |

The Python port reproduces the old IntPath database exactly for all four of its organisms (human, mouse, yeast,
*M. tuberculosis*): every related-pathway pair, integrated pathway name and pathway-gene row.
The first public IntPathV2 release (open tier, October 2026) covers human, mouse, yeast and
*M. tuberculosis*. For human it has 3,747 integrated Reactome/WikiPathways pathways, 10,724 GO
sets, 20,952 MSigDB sets and 1.2 million merged PPI edges. See [Docs/METHODS.md](Docs/METHODS.md) and
[Docs/OLD_INTPATH.md](Docs/OLD_INTPATH.md).

## Quick start

```bash
pip install -e '.[web,test]'
intpath organisms                         # supported organisms
intpath build sapiens                     # download sources and build a release in Data/intpathv2/release/sapiens
intpath build sapiens --public            # only openly redistributable sources (no KEGG/BioCyc)
intpath ora   Data/intpathv2/release/sapiens genes.txt --top 20
intpath pairs Data/intpathv2/release/sapiens genes.txt --with-ppi
intpath gsea  Data/intpathv2/release/sapiens ranking.tsv --collections pathway,GO:BP
intpath serve                             # web UI + REST API at http://127.0.0.1:8000 (docs at /docs)
pytest
```

## Repository layout

```
src/intpath/   names.py (name alignment + union-find)  unify.py (full unification)
               mapping.py (gene IDs)  organisms.py  sources.py (KEGG/Reactome/WikiPathways/BioCyc)
               ppi.py  go.py  enrich.py (ORA, pairs, GSEA, themes)  build.py  cli.py
web/           FastAPI service and single-page UI for intpath.genohub.org
Docs/          METHODS, DATA_FORMATS, ROADMAP, DEPLOY, OLD_INTPATH
legacy/java/   old IntPath Java sources: human pipeline, tools/ (alignment, DB writers), organisms/<org>/
stats/         summary statistics of builds
```

## Data licensing

The code is public. Built data inherit their sources' licences. WikiPathways, Reactome, GO,
NCBI Gene, HGNC, STRING, BioGRID and IntAct are open (CC0, CC-BY or MIT). **KEGG** and
**BioCyc** content may not be redistributed without a licence, so public downloads are built
with `--public`.
