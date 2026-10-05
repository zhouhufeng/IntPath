# IntPathV2: methods

IntPathV2 keeps the integration method published for the old IntPath and applies the same
principles to two new data types, protein-protein interactions (PPIs) and Gene Ontology (GO).
It also adds an enrichment engine that uses the integrated, provenance-tagged data.

> Zhou H, Jin J, Zhang H, Yi B, Wozniak M, Wong L. **IntPath: an integrated pathway gene
> relationship database for model organisms and important pathogens.** *BMC Systems Biology*
> 2012, 6(Suppl 2):S2. doi:10.1186/1752-0509-6-S2-S2 (PMC3521174)

## 1. Organisms

IntPath has always covered more than one organism. Old IntPath shipped *H. sapiens*, *M. musculus*,
*S. cerevisiae* and *M. tuberculosis* H37Rv. IntPathV2 has a single organism registry
(`src/intpath/organisms.py`). Each entry holds the organism's NCBI taxid, KEGG code, Reactome
species, WikiPathways species, STRING taxid and GO annotation file. Every builder takes one of
these entries as a parameter. The registry ships with:

| key | organism | pathways | GO | PPI |
|---|---|---|---|---|
| sapiens (in old IntPath) | *Homo sapiens* | KEGG, Reactome, WikiPathways, (BioCyc*) | goa_human | STRING, BioGRID, IntAct |
| musculus (in old IntPath) | *Mus musculus* | KEGG, Reactome, WikiPathways | mgi | STRING, + |
| cerevisiae (in old IntPath) | *S. cerevisiae* | KEGG, Reactome, WikiPathways | sgd | STRING, + |
| tuberculosis (in old IntPath) | *M. tuberculosis* H37Rv | KEGG, WikiPathways | (UniProt GOA) | STRING |
| norvegicus | *Rattus norvegicus* | KEGG, Reactome, WikiPathways | rgd | STRING |
| rerio | *Danio rerio* | KEGG, Reactome, WikiPathways | zfin | STRING |
| melanogaster | *D. melanogaster* | KEGG, Reactome, WikiPathways | fb | STRING |
| elegans | *C. elegans* | KEGG, Reactome, WikiPathways | wb | STRING |
| thaliana | *A. thaliana* | KEGG, WikiPathways | tair | STRING |
| coli | *E. coli* K-12 MG1655 | KEGG | ecocyc | STRING |

\* BioCyc requires a subscription licence (since 2017). Pass a local `pathways.col` to use it.
To add another organism, add one registry entry or pass a JSON file to `organisms.load_custom`.

## 2. Gene identifier normalisation

Each source's identifiers (NCBI Gene IDs, UniProt, Ensembl, locus tags, MOD IDs, symbols and
synonyms) are mapped to the organism's official symbol. The base layer is NCBI `gene_info`,
which exists for every organism. Human adds the HGNC layer, and any organism can add a UniProt
`idmapping_selected` file. Previous and alias symbols are used only when they resolve to
exactly one gene. Ambiguous identifiers are counted and reported, never guessed. This follows
the old IntPath "no introduced noise" rule.

## 3. Pathway integration (old IntPath method, reproduced exactly)

Pathway sources: KEGG (REST + KGML), **Reactome** (current release, 97 at the first build:
all-levels membership, the pathway hierarchy, and gene pairs from Reactome's interaction
export), WikiPathways (GMT + GPML) and BioCyc (licensed). Reactome gene pairs are assigned to
pathways as follows:
- **Complex context:** the pair becomes GPrel, placed in that complex's pathways
  (`Complex_2_Pathway`).
- **Reaction context:** the pair becomes PPrel, placed in the lowest-level pathways that
  contain both genes.

Either way, the pair is then propagated to every ancestor pathway.

1. **Alignment.** Each pair of pathway names is aligned case-insensitively by longest common
   subsequence (LCS). *score* = LCS length;
   *ratio* = 2·score / (len(a) + len(b)). The LCS uses a bit-parallel implementation
   (Hyyrö 2004).
2. **Best hit.** Every name in list *x* is compared with all names in list *y*, and only the best
   hit (the highest ratio, with the first one winning ties) is kept. Comparisons run between
   every pair of sources (old IntPath order: KEGG→BioCyc, KEGG→WikiPathways, BioCyc→WikiPathways) and
   within each source.
3. **Acceptance.** A best hit is accepted when
   `(score > len(shorter) − 1 and ratio ≥ 0.5) or ratio > 0.91`.
   The pair is then rejected if it is on the curated mismatch list ("T cell" vs "B cell",
   "NOD" vs "Toll", "Linoleic" vs "Lipoic", …).
4. **Grouping.** Accepted pairs go into a disjoint set (union-find). Each connected component
   becomes one integrated pathway, named after its shortest member name with variant tokens
   (I, II, Type, …) removed.
5. **Full unification.** Every gene and every gene pair of every member pathway is kept.
   Each gene membership records the sources that support it. Each gene pair keeps the union of
   its unified relation types (ECrel, PPrel, GErel, GPrel) and the sources that support it.
   Pathways that match nothing are carried over unchanged.

**Validation.** Old IntPath kept one hand-curated rule set per organism: mismatch list, name tokens to
drop, and name replacements, in `src/intpath/curation.py`. With these rules,
`intpath legacy --organism <org>` reproduces every old IntPath release exactly: the related-pathway
pairs, the integrated pathway names and every pathway-gene row (`tests/test_names.py`).

| Old IntPath organism | related pairs | integrated pathways | genes | pathway-gene rows |
|---|---|---|---|---|
| *H. sapiens* | 87 | 582 | 7,134 | 23,873 |
| *M. musculus* | 129 | 555 | 8,013 | 24,878 |
| *S. cerevisiae* | 136 | 285 | 1,833 | 5,285 |
| *M. tuberculosis* H37Rv | 53 | 299 | 1,146 | 4,778 |

The original `tools/sequenceAlignment.java` is a global alignment (match +1, case-insensitive;
mismatch −1; gap 0). A mismatch never beats two gaps, so its score is exactly the LCS length.
When two candidate names had equal length, old IntPath chose the integrated name by Java `HashMap`
iteration order. Legacy mode emulates the JDK ≤ 7 `HashMap` (`names.java_hashmap_order`),
which reproduces those choices. IntPathV2 mode breaks such ties deterministically instead.

### IntPathV2 extensions (default for new builds; `legacy=True` turns them off)

- **Name cleaning** strips HTML (`<i>myo</i>-inositol`) and species suffixes
  (" - Homo sapiens (human)").
- **Mismatch list = the union of all old IntPath organisms' curated lists** (49 pairs), plus IntPathV2 additions
  (VEGF/EGF, biosynthesis/degradation, positive/negative regulation, …).
- **Entity guard.** Names that refer to different molecular entities or subtypes are not
  merged: "Signaling by FGFR1" vs "FGFR2", "IL-1" vs "IL-17", "RHOC" vs "RHOG",
  "hyperlipidemia type 1" vs "type 2", "Type I" vs "Type II diabetes". Free-standing variant
  suffixes ("… I" vs "… II") are not affected, so the old IntPath BioCyc-variant merges still happen.
- **Gene-overlap guard.** A name match is kept only if the gene sets have Jaccard ≥ 0.1 or the
  cleaned names are identical. This is needed now that Reactome adds about 2,800 human pathway
  names. Jaccard, not the overlap coefficient, is used so that a superpathway (Reactome
  "Disease") can't absorb its children. In old IntPath, for example, "Steroid biosynthesis" was merged
  with "Steroid hormone biosynthesis".
- **No within-source merging for hierarchical sources** (Reactome, GO). Similar names inside
  Reactome are curated siblings ("RHOC/RHOG GTPase cycle", "RNA Polymerase I/II/III
  Transcription"), not duplicates.
- **Hierarchy-aware merging.** If a Reactome pathway and one of its ancestors end up in the
  same group through name chaining, the descendant is detached into its own set. Every
  Reactome parent/child pair is recorded as a link between the integrated sets.
- **Stable IDs.** Each integrated pathway's ID (`IP` + SHA-1 of its sorted member IDs) stays the
  same across rebuilds as long as its members don't change.

## 4. PPI integration

| source | file | used |
|---|---|---|
| STRING v12.0 | `protein.physical.links` (+ `protein.info`) | physical subnetwork, combined score ≥ 700 |
| BioGRID (5.0.262 in the first builds) | `BIOGRID-ORGANISM-LATEST.tab3.zip`, organism member | physical experimental systems, same-species pairs |
| IntAct and MINT | `intact.zip` (PSI-MITAB 2.7, all IMEx evidence) | same-species pairs. Rows whose source database is MINT are kept as source **MINT**; MINT curates into IntAct under IMEx, so this is the current form of MINT. |
| HuRI | `interactome-atlas.org/data/HuRI.tsv` | human binary reference interactome (Luck et al. 2020) |

- **Identifiers:** STRING ENSP ids, BioGRID Entrez ids, IntAct UniProt accessions and HuRI
  Ensembl genes are mapped to official symbols. UniProt accessions resolve through UniProt's
  per-organism ID mapping (or its REST stream when no file exists).
- **Unification:** edges are stored undirected and *fully unified*. Each edge keeps:
  - its supporting sources;
  - the number of distinct publications;
  - the detection methods;
  - the best STRING score;
  - a **confidence tier**:
    - **high:** ≥ 2 experimental sources (BioGRID / IntAct / MINT / HuRI), ≥ 2 publications,
      or STRING ≥ 900;
    - **medium:** one experimental source or STRING ≥ 700.
- **Tier use:** tiers are stored, never used to drop edges; users choose one at analysis time.
- **Pathway overlay:** each edge also records the integrated pathways in which the pair is a
  curated relation ("pathway-supported" edges).

## 5. GO integration

- GO is parsed from `go-basic.obo` (is_a and part_of edges) together with the organism's GAF.
  NOT-qualified and ND annotations are dropped. Annotations are propagated to all ancestors
  (true-path rule), and evidence codes are kept as provenance.
- **Lossless redundancy merge.** GO terms in the same namespace whose propagated gene sets are
  identical *and* that lie on one ancestor/descendant chain are merged into one gene set. The
  set is named after the most specific term, and all term IDs are kept as members. Identical
  terms that are not related by lineage stay separate.
- **Cross-links.** A GO set and an integrated pathway are linked when Jaccard ≥ 0.5, or when
  their names pass the IntPath alignment rule and Jaccard ≥ 0.2. GO terms and pathways describe
  different kinds of concept, so they stay separate gene sets. The links are used to collapse
  redundant enrichment hits.

## 6. MSigDB (companion library)

MSigDB (v2026.1 Hs and Mm at the first build) is used next to the integrated data, not merged
into it:
- **Gene sets:** collections IntPathV2 has no other source for become gene sets with
  collection `msigdb:<collection>`. Human: H, C1, C2:CGP, C2:CP:PID, C3, C4, C5:HPO, C6, C7,
  C8, C9. Mouse: MH, M1, M2:CGP, M3, M5:MPT, M7, M8.
- **Equivalence links only:** MSigDB's copies of Reactome, WikiPathways, KEGG and GO are not
  stored again. Each one is linked through its `exactSource` to the IntPathV2 set holding that
  source pathway, with the Jaccard between the two versions. For human Reactome, 1,830 of
  1,839 sets map, with median Jaccard 1.0.
- **Restricted collections:** BioCarta, KEGG_LEGACY and KEGG_MEDICUS have extra licence terms
  and are built only into the full tier.

## 7. Storage

- **Release database:** each organism release is one SQLite file (`intpath.sqlite`, schema in
  `src/intpath/db.py` and DATA_FORMATS.md). It holds genes, aliases, sets, memberships with
  provenance, gene pairs, links, PPI edges with evidence, MSigDB equivalents and the merge log.
- **Web service:** loads only a compact in-memory library from it (memberships, pathway source
  support, gene pairs, links, PPI adjacency per tier). Everything else is read per request.
- **Tiers:** there are two release trees. `release-open/` has no KEGG, BioCyc or restricted
  MSigDB collections and is redistributable; it's what intpath.genohub.org serves.
  `release/` is the full tier, for internal use.

## 8. Enrichment analysis

| method | question | statistics |
|---|---|---|
| `ora` | Are my genes over-represented in a set? | One-sided hypergeometric test; universe = library genes ∩ optional background; BH-FDR; fold enrichment. Each hit also reports **consensus**: the share of overlapping genes that ≥ 2 source databases place in the pathway. |
| `pairs` | Do the curated relations among my genes concentrate in a pathway's wiring? | Hypergeometric test over gene pairs: universe = all pathway gene pairs (optionally plus merged PPI edges inside each pathway); "query pairs" = pairs with both genes in the input. This is the gene-pair level that only IntPath's integrated relations allow. |
| `gsea` | Is a set enriched at the top or bottom of a ranked list? | Preranked GSEA (weighted KS, p = 1); permutation nulls matched by gene-set size and shared across sets; NES; nominal p; GSEA FDR (Subramanian 2005) and BH. |
| `cluster` | Which hits are the same signal? | Significant hits are grouped greedily into themes, using IntPath merge groups, GO↔pathway links and gene-set Jaccard ≥ 0.5. |

The 2012 "Analyze Distance" tool (Floyd-Warshall distances between enriched pathways over the
STRING network) is planned for IntPathV2 on the merged PPI network (see ROADMAP).
