# IntPath V3 release files

All files are tab-delimited UTF-8 with a header row, one directory per organism
(`release/<organism>/`).

| file | columns |
|---|---|
| `intpath_genesets.tsv` | `set_id`, `name`, `collection` (pathway, GO:BP, GO:MF, GO:CC), `sources` (comma list), `n_genes`, `n_pairs`, `links` (related set IDs) |
| `intpath_set_genes.tsv` | `set_id`, `gene`, `sources`: the source databases that put the gene in the set (for GO, the evidence codes) |
| `intpath_set_genepairs.tsv` | `set_id`, `gene_a`, `gene_b`, `relations` (ECrel/PPrel/GErel/GPrel), `sources` |
| `intpath_set_members.tsv` | `set_id`, `source`, `source_pathway` (original pathway ID or name, or GO term ID) |
| `intpath.gmt` | GMT export: `set_id`, `name`, genes… |
| `intpath_ppi.tsv` | `gene_a`, `gene_b`, `sources`, `n_sources`, `n_pmids`, `methods`, `string_score`, `pathway_sets` |
| `related_pathways.tsv` | the accepted name matches: `source_a`, `pathway_a`, `source_b`, `pathway_b`, `align_score`, `align_ratio`, `overlap_coef` |
| `stats.json` | build date, per-source and integrated statistics, gene-mapping statistics |

## Unified gene relationships (from V2, Table 1)

| type | meaning | from |
|---|---|---|
| ECrel | enzyme-enzyme relation, successive reaction steps | KEGG ECrel, BioCyc SEQUENTIAL_CATALYSIS, GPML mim-conversion |
| PPrel | protein-protein interaction (binding, modification, shared control) | KEGG PPrel/PCrel, BioCyc CO_CONTROL/INTERACTS_WITH, GPML interactions |
| GErel | gene expression interaction (TF → target) | KEGG GErel, GPML mim-transcription-translation |
| GPrel | same complex or group, not necessarily direct | KEGG group entries, GPML groups, BioCyc IN_SAME_COMPONENT |

## V2 legacy files (`Data/<organism>/integrated/IntPathData`)

- `*IntPathGenes`: `pathway`, `gene`, `sources`
- `*IntPathGenePairs`: `geneA`, `geneB`, `relation(s)`, `pathway`, `sources`
- `*GroupGenes`: `pathway`, `gene`, `group id`, `source`
