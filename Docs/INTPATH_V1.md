# IntPath v1 (2012)

IntPath v1 is the database and methods published in 2012:

> Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L. IntPath: an integrated pathway gene
> relationship database for model organisms and important pathogens. *BMC Systems Biology*
> 2012, 6(Suppl 2):S2. https://doi.org/10.1186/1752-0509-6-S2-S2

It covered four organisms. The IntPathV2 Python port reproduces every one exactly in legacy
mode (see [METHODS.md](METHODS.md)):

| organism | source pathways | related pathway pairs | integrated pathways | genes | pathway-gene rows | gene pairs |
|---|---|---|---|---|---|---|
| *H. sapiens* | 661 | 87 | 582 | 7,134 | 23,873 | 50,852 |
| *M. musculus* | 678 | 129 | 555 | 8,013 | 24,878 | 64,029 |
| *S. cerevisiae* | 406 | 136 | 285 | 1,833 | 5,285 | 3,956 |
| *M. tuberculosis* H37Rv | 351 | 53 | 299 | 1,146 | 4,778 | 6,227 |

## Human build

- 87 related pathway pairs formed 57 integrated pathway groups covering 136 source pathways.
- Result: **582 integrated pathways**, **7,134 unique genes**, 23,873 pathway-gene rows and
  **50,852 unique gene pairs**.
- Genes overlapped much more between sources than gene pairs did, which is why full
  unification adds so many relations.
- The curated mismatch list rejected look-alike names such as NOD-like vs Toll-like receptor
  signalling, Linoleic vs Lipoic acid metabolism, and Glycerolipid vs Glycerophospholipid
  metabolism.

IntPath v1 also used a PPI network (human: 214,957 interactions) for its "Analyze Distance"
tool, which computed Floyd-Warshall distances between enriched pathways.

## Code (`legacy/java/`)

| file | role |
|---|---|
| `KEGG.java`, `Wiki.java`, `BioCyc.java` | extraction of pathway genes, gene pairs and groups from each source format |
| `Normalize.java` | gene ID mapping and relation unification |
| `Integration.java` | name alignment, union-find grouping, full unification |
| `GeneListsAnalyse.java`, `pathwayItem.java` | hypergeometric pathway identification and pathway distances |
| `API.java`, `stats.java`, `main.java`, `Debugging.java` | web-service API, statistics, pipeline driver |

`tools/` holds the shared helpers (`sequenceAlignment`, `calc`, `conts`, `createTables`,
`test`). `organisms/<org>/` holds the per-organism variants of the pipeline for mouse, yeast and
*M. tuberculosis*. They mainly differ in their curated mismatch lists and naming rules, which
`src/intpath/curation.py` now holds.
