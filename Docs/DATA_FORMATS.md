# IntPathV2 data files

`intpath build <organism>` writes one directory per organism (`Data/intpathv2/db/<organism>/`).
All tables are tab-delimited UTF-8 with a header row.

| file | content |
|---|---|
| `intpath.sqlite` | the database served by the web service: tables `meta`, `gene`, `alias`, `gset`, `set_member`, `set_gene`, `set_pair`, `set_link`, `ppi`, `string_ppi`, `msigdb_equivalent`, `match_log`, `diagram`, `diagram_of` (see `src/intpath/db.py`) |
| `intpath_genesets.tsv` | `set_id`, `name`, `collection` (pathway, GO:BP, GO:MF, GO:CC, msigdb:<collection>), `sources` (databases the set's member pathways come from), `n_genes`, `n_pairs`, `links` (related set IDs) |
| `intpath_set_genes.tsv` | `set_id`, `gene`, `sources` (supporting databases; for GO, the evidence codes) |
| `intpath_set_genepairs.tsv` | `set_id`, `gene_a`, `gene_b`, `relations` (ECrel/PPrel/GErel/GPrel), `sources` |
| `intpath_set_members.tsv` | `set_id`, `source`, `source_pathway` (original pathway ID, or GO term ID); names and descriptions of member pathways are in `intpath.sqlite` (`set_member`) |
| `intpath.gmt` | GMT export: `set_id`, `name`, genes… |
| `intpath_ppi.tsv` | merged physical PPI edges: `gene_a`, `gene_b`, `sources`, `n_sources`, `n_pmids` (publications), `methods`, `string_score`, `tier`, `pathway_sets` (IntPath pathways containing the pair) |
| `related_pathways.tsv` | the accepted related pathway pairs that were merged: `source_a`, `pathway_a`, `source_b`, `pathway_b`, `align_score`, `align_ratio`, `jaccard` (gene overlap, informational only) |
| `merge_review.tsv` | every candidate pathway pair with its decision (`accept`, `reject`, `pending`) and the reason |
| `stats.json` | build date and summary counts |

## Gene relationship types (IntPath v1, Table 1)

| type | meaning |
|---|---|
| ECrel | enzyme-enzyme relation: successive reaction steps |
| PPrel | protein-protein interaction: binding, modification, shared control |
| GErel | gene expression interaction: transcription factor → target |
| GPrel | same complex or group, not necessarily direct |

## IntPath v1 files (`Data/<organism>/integrated/IntPathData`)

- `*IntPathGenes`: `pathway`, `gene`, `sources`
- `*IntPathGenePairs`: `geneA`, `geneB`, `relation(s)`, `pathway`, `sources`
- `*GroupGenes`: `pathway`, `gene`, `group id`, `source`
