# IntPathV2 roadmap (intpath.genohub.org)

## Done in this repository
- [x] Python port of the old IntPath merge method (`intpath.names`, `intpath.unify`), validated to
      reproduce the old IntPath human release exactly (`tests/test_names.py`, `intpath legacy`)
- [x] Multi-organism registry with 10 organisms, including the 4 old IntPath organisms
      (`intpath.organisms`)
- [x] Gene ID mapping for any organism (NCBI gene_info, plus HGNC and UniProt layers)
- [x] Source extractors: KEGG REST+KGML (gene pairs), Reactome, WikiPathways GMT+GPML
      (gene pairs), BioCyc `pathways.col` (licensed)
- [x] PPI integration: STRING, BioGRID, PSI-MITAB, TSV; evidence-preserving union;
      pathway overlay
- [x] GO integration: propagation, lossless redundancy merge, GO↔pathway cross-links
- [x] Enrichment: ORA with source consensus, gene-pair/network enrichment, preranked GSEA,
      redundancy-aware themes
- [x] Release format (TSV + GMT + stats.json), CLI, FastAPI service with a single-page UI
- [x] Public-repo sync script (`scripts/sync_public.sh`)
- [x] First IntPathV2 human build (2026-10-05): KEGG 372 + Reactome 2,835 + WikiPathways 991 pathways
      → 4,003 integrated sets (130 merged from 315 source pathways), 14,849 genes,
      173,645 gene pairs; 10,721 GO sets (1,020 GO↔pathway links); 85,998 STRING physical
      edges (score ≥ 700), 23,175 of them also curated pathway relations
      (`Data/intpathv2/release/sapiens/stats.json`)

- [x] Reactome topology (interaction export) and hierarchy-aware merging
- [x] PPIs: STRING, BioGRID, IntAct, MINT (via IntAct/IMEx), HuRI, with confidence tiers
- [x] MSigDB v2026.1 companion library and equivalence links
- [x] SQLite release database and compact serving library
- [x] Deployed at https://intpath.genohub.org (open tier: human, mouse, yeast, *M. tuberculosis*),
      2026-10-05

## Known gaps
- *M. tuberculosis* has no open-licence pathway source (KEGG is full tier only). Its public
  release has GO (from UniProt) and PPIs only.
- Mouse and yeast Reactome are inferred by orthology and have no interaction export, so they
  carry membership only. Yeast WikiPathways is small.
- The full tier (with KEGG) is built for human only and isn't served publicly. Serving it
  needs a KEGG licence (plan D4).
- Not done yet: CORUM/Complex Portal complexes, SIGNOR, GO evidence variants (no-IEA), and
  multilevel GSEA p-values.

## Next
1. **Builds for all registry organisms.** Run `intpath build <org>` for each organism, then
   review `related_pathways.tsv` by hand. Old IntPath also needed moderate manual curation; record new
   cases in the mismatch lists.
2. **Curation file.** Move the mismatch lists and manual merge/split decisions into a versioned
   `curation/<organism>.tsv` that the builder reads.
3. **Licensing.** KEGG- and BioCyc-derived files can't be redistributed without a licence.
   Either get a KEGG FTP/redistribution licence, or serve KEGG-derived content only through
   the web analysis and publish the downloadable release from the `--public` build
   (Reactome + WikiPathways + GO + PPIs).
4. **More sources.** Pathway Commons / PathBank / PANTHER pathways, MSigDB Hallmarks as a
   reference collection, CORUM complexes (GPrel), SIGNOR causal relations (directed PPrel/GErel).
5. **Analyze Distance (old IntPath tool).** Shortest-path distances between pathways on the merged PPI
   network, plus network propagation (random walk with restart) from the query genes.
6. **Topology-aware enrichment.** Weight pathway genes by betweenness in the unified gene-pair
   graph (SPIA/ROntoTools-like) and compare the results with the plain pair test.
7. **Ortholog projection** for organisms without curated pathways (e.g. Reactome-inferred,
   or orthologs via NCBI gene_orthologs), with every projected membership labelled with its
   source.
8. **Web.** Interactive pathway graph (Cytoscape.js), downloadable result tables, job queue
   for large GSEA runs, result permalinks, and an API key or rate limit for the public server.
9. **Benchmark.** Compare against single-source libraries, MSigDB C2/C5 and Enrichr libraries
   on GEO/KEGG disease benchmark sets (e.g. GSEABenchmarkeR), measuring sensitivity,
   prioritisation and redundancy. Results go into the IntPathV2 manuscript.
10. **Release cadence.** Quarterly automated builds, with the date and source versions recorded
    in `stats.json`, and a changelog of added, removed and re-grouped pathways keyed by stable
    ID.
