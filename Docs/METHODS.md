# IntPathV2: methods

IntPathV2 keeps the integration method of IntPath v1 and applies it to much larger, current
pathway data. It adds a merged physical protein-protein interaction (PPI) network, Gene Ontology
(GO) and MSigDB gene sets, pathway maps, and an enrichment engine on the integrated data.

> Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L. **IntPath: an integrated pathway gene
> relationship database for model organisms and important pathogens.** *BMC Systems Biology*
> 2012, 6(Suppl 2):S2. doi:10.1186/1752-0509-6-S2-S2 (PMC3521174)

## 1. Organisms

IntPathV2 is built for the four IntPath v1 organisms: *H. sapiens*, *M. musculus*,
*S. cerevisiae* and *M. tuberculosis* H37Rv. One organism registry (`src/intpath/organisms.py`)
holds each organism's identifiers and annotation files, and every builder takes a registry entry
as a parameter. More organisms can be added with one registry entry, or with a JSON file passed
to `organisms.load_custom`.

## 2. Gene identifier normalisation

Every identifier type found in the source data (NCBI Gene IDs, UniProt, Ensembl, locus tags,
model-organism IDs, symbols and synonyms) is mapped to the organism's official gene symbol. The
base layer is NCBI `gene_info`, which exists for every organism. Human adds the HGNC layer, and
any organism can add a UniProt ID-mapping layer. Previous and alias symbols are used only when
they resolve to exactly one gene. Ambiguous identifiers are counted and reported, never guessed.
This follows the IntPath v1 rule: "no deletion, no introduced noise".

## 3. Pathway integration

Each source pathway gives a set of genes and a set of gene pairs. Gene pairs carry the four
unified relation types of IntPath v1:

| type | meaning |
|---|---|
| ECrel | enzyme-enzyme relation: successive reaction steps |
| PPrel | protein-protein interaction: binding, modification, shared control |
| GErel | gene expression interaction: transcription factor → target |
| GPrel | same complex or group, not necessarily direct |

Related pathways are found and merged in five steps.

1. **Name alignment (IntPath v1).** Each pair of pathway names is aligned case-insensitively
   by longest common subsequence (LCS): *score* = LCS length, *ratio* = 2·score / (len(a) +
   len(b)). Every name keeps only its best hit in each list. Comparisons run between every pair
   of sources and within each source. A best hit is a candidate when
   `(score > len(shorter) − 1 and ratio ≥ 0.5) or ratio > 0.91`, and it is not on the curated
   mismatch list ("T cell" vs "B cell", "NOD" vs "Toll", "Linoleic" vs "Lipoic", …).
2. **Review of candidates (replaces the manual review of IntPath v1).** Each candidate pair is
   judged on whether the two names describe the same biological process. Pairs with identical
   names after cleaning are accepted directly. Pairs whose names only look alike are rejected:
   different genes or subtypes ("FGFR1" vs "FGFR2", "IL-1" vs "IL-17"), a general process vs a
   part or variant of it ("Cell cycle" vs "Cell cycle checkpoints"), or a different meaning
   ("Prion disease" vs "Parkinson disease").
3. **Review of all pathway names by meaning.** Letter alignment misses synonyms with different
   spellings ("Citrate cycle (TCA cycle)" vs "Citric acid cycle"). For every pathway name, the
   most similar names are therefore also listed by shared meaningful words, and each pair is
   judged by meaning only. Gene overlap is not used. Accepted synonyms are merged even when the
   alignment never proposed them.
4. **Grouping.** Accepted pairs go into a disjoint set (union-find). Each connected component
   becomes one IntPath pathway, named after its shortest member name (a trailing variant
   numeral such as "II" is dropped).
5. **Full unification.** Every gene and every gene pair of every member pathway is kept. Each
   gene pair keeps the union of its relation types. Pathways that match nothing are carried
   over unchanged.

Every review decision is stored with its reason in `src/intpath/data/merge_curation.tsv`, and
every build writes a log of all candidate pairs with their decisions (`merge_review.tsv`). A
pair that has never been reviewed is not merged; it is logged as `pending` until it is reviewed.

Each IntPath pathway records:
- the original pathway names it was merged from, with a short description of what each means;
- links to related IntPath pathways (for example a pathway and its sub-pathways);
- a stable ID (`IP` + SHA-1 of its sorted member IDs) that stays the same across rebuilds as
  long as its members don't change.

**Validation.** IntPath v1 kept one hand-curated rule set per organism: mismatch list, name
tokens to drop, and name replacements (`src/intpath/curation.py`). With these rules,
`intpath legacy --organism <org>` reproduces every IntPath v1 release exactly: the related
pathway pairs, the integrated pathway names and every pathway-gene row (`tests/test_names.py`).

| IntPath v1 organism | related pairs | integrated pathways | genes | pathway-gene rows |
|---|---|---|---|---|
| *H. sapiens* | 87 | 582 | 7,134 | 23,873 |
| *M. musculus* | 129 | 555 | 8,013 | 24,878 |
| *S. cerevisiae* | 136 | 285 | 1,833 | 5,285 |
| *M. tuberculosis* H37Rv | 53 | 299 | 1,146 | 4,778 |

The v1 alignment (`tools/sequenceAlignment.java`) is a global alignment with match +1, mismatch
−1 and gap 0. A mismatch never beats two gaps, so its score is exactly the LCS length. When two
candidate names had equal length, v1 chose the integrated name by Java `HashMap` iteration
order. Legacy mode emulates the JDK ≤ 7 `HashMap` (`names.java_hashmap_order`) to reproduce
those choices.

## 4. Physical PPI network

Physical interactions from several experimental interaction databases are merged into one
network:
- **Identifiers:** every database's identifiers are mapped to official symbols, as in section 2.
- **Unification:** edges are undirected and fully unified. Each edge keeps the number of
  supporting databases, the number of distinct publications and the detection methods.
- **Confidence tiers:** **high** = supported by ≥ 2 databases or ≥ 2 publications;
  **medium** = one database. Tiers are stored, never used to drop edges; users choose a tier at
  analysis time.
- **Pathway overlay:** each edge records the IntPath pathways in which the pair is also a
  pathway gene pair.

A predicted functional-association network is kept separately and is never mixed into the
physical network.

## 5. GO

- GO is read from `go-basic.obo` (is_a and part_of edges) with the organism's annotation file.
  NOT-qualified and ND annotations are dropped. Annotations are propagated to all ancestors
  (true-path rule).
- **Lossless redundancy merge:** GO terms in the same namespace whose propagated gene sets are
  identical *and* that lie on one ancestor/descendant chain become one gene set. The set is
  named after the most specific term, and all term IDs are kept.
- **Links:** a GO set and an IntPath pathway are linked when their gene sets are very similar.
  GO terms stay separate gene sets; the links are used to group redundant enrichment hits.

## 6. MSigDB

MSigDB gene sets (human and mouse) are a separate library next to the pathways, not merged into
them. MSigDB copies of pathways and GO terms that IntPath already holds are not stored again;
each is linked to the IntPath set that holds the same pathway.

## 7. Storage

Each organism is one SQLite database (`Data/intpathv2/db/<organism>/intpath.sqlite`, schema in
`src/intpath/db.py` and [DATA_FORMATS.md](DATA_FORMATS.md)). It holds genes, aliases, gene
sets, memberships, gene pairs, links, PPI edges, pathway descriptions, the merge log and the
pathway drawings. The web service keeps a compact in-memory index of it and reads everything
else per request.

## 8. Enrichment analysis

| method | question | statistics |
|---|---|---|
| `ora` | Are my genes over-represented in a gene set? | One-sided hypergeometric test; universe = library genes ∩ optional background; BH-FDR; fold enrichment |
| `pairs` | Do the relations among my genes concentrate in a pathway's wiring? | Hypergeometric test over gene pairs: universe = all pathway gene pairs (optionally plus PPI edges inside each pathway); query pairs = pairs with both genes in the input |
| `gsea` | Is a gene set enriched at the top or bottom of a ranked list? | Preranked GSEA (weighted KS, p = 1); permutation nulls matched by set size and shared across sets; NES; nominal p; GSEA FDR and BH |
| `cluster` | Which hits are the same signal? | Significant hits are grouped into themes using IntPath merges, GO-pathway links and gene-set similarity |

## 9. Pathway maps

- **Map view:** each pathway with a curated diagram is drawn in the browser (Cytoscape.js)
  from the diagram's layout data: boxes, positions, arrows and compartments. The user's genes are
  highlighted, or coloured by score for GSEA. Genes of the IntPath pathway that are not on that
  diagram are listed beside it.
- **Network view:** the merged gene pairs of any IntPath pathway are drawn as a network, with
  optional PPI edges. It can be exported as SBML-qual or SBGN-ML.
