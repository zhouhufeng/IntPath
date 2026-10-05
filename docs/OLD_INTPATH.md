# Old IntPath: archived builds

Old IntPath covered four organisms. Every one is reproduced exactly by the IntPathV2 Python port in
legacy mode (see METHODS.md):

| organism | KEGG / BioCyc / WikiPathways pathways | related name pairs | integrated pathways | genes | pathway-gene rows | gene pairs |
|---|---|---|---|---|---|---|
| *H. sapiens* | 237 / 289 / 135 | 87 | 582 | 7,134 | 23,873 | 50,852 |
| *M. musculus* | 218 / 322 / 138 | 129 | 555 | 8,013 | 24,878 | 64,029 |
| *S. cerevisiae* | 98 / 183 / 125 | 136 | 285 | 1,833 | 5,285 | 3,956 |
| *M. tuberculosis* H37Rv | 110 / 233 / 8 | 53 | 299 | 1,146 | 4,778 | 6,227 |

The rest of this page details the human build (2021 rerun of the 2012 pipeline; `Statistics`
and `IntegrationStatistics`). Gene counts below are pathway-gene rows unless stated otherwise.

## Sources after ID normalisation

| source | pathways (genes / gene pairs) | pathway-gene rows | gene pairs | genes per pathway | pairs per pathway |
|---|---|---|---|---|---|
| KEGG | 237 / 208 | 17,159 | 40,589 | 72.4 | 171.3 |
| BioCyc (HumanCyc) | 290 / 257 | 2,085 | 9,576 | 7.2 | 33.0 |
| WikiPathways | 135 / 126 | 6,247 | 22,435 | 46.3 | 166.2 |

Overlap of unique genes between sources: KEGG∩WikiPathways 2,485 (Jaccard 0.36),
KEGG∩BioCyc 824 (0.13), WikiPathways∩BioCyc 396 (0.10). Gene pairs overlap much less
(Jaccard 0.02–0.04). This is why full unification adds so many relations.

## Integration

| step | KEGG–Wiki | KEGG–BioCyc | Wiki–BioCyc | KEGG–KEGG | BioCyc–BioCyc | Wiki–Wiki |
|---|---|---|---|---|---|---|
| related name pairs | 29 | 3 | 12 | 5 | 34 | 4 |

- 87 related pairs → 57 integrated pathway groups covering 136 source pathways
- Result: **582 integrated pathways**, **7,134 unique genes**, 23,873 pathway-gene rows,
  **50,852 unique gene pairs**
- The curated mismatch list rejected look-alike names such as NOD-like vs Toll-like receptor
  signalling, Linoleic vs Lipoic acid metabolism, and Glycerolipid vs Glycerophospholipid
  metabolism

## PPI

STRING v9 human, combined score ≥ 750, mapped to IntPath symbols: 214,957 interactions. Old IntPath
used them for the "Analyze Distance" tool (Floyd-Warshall pathway distances).

## Code (`legacy/java/`)

| file | role |
|---|---|
| `KEGG.java`, `Wiki.java`, `BioCyc.java` | KGML / GPML / BioPAX extraction of pathway genes, gene pairs and groups |
| `Normalize.java` | ID mapping (NCBI, KEGG, UniProt, HGNC, BioMart) and relation unification |
| `Integration.java` | name alignment, union-find grouping, full unification |
| `GeneListsAnalyse.java`, `pathwayItem.java` | hypergeometric pathway identification and pathway distances |
| `API.java`, `stats.java`, `main.java`, `Debugging.java` | web-service API, statistics, pipeline driver |

`tools/` holds the shared helpers (`sequenceAlignment`, `calc`, `conts`, `createTables`,
`test`). `organisms/<org>/` holds the per-organism variants of the pipeline for mouse, yeast and
*M. tuberculosis*. These mainly differ in their curated mismatch lists and naming rules, which
`src/intpath/curation.py` now holds.
