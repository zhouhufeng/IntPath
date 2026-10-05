"""Gene identifier normalisation to official gene symbols (IntPathV2 rewrite of Normalize.java).

Works for any organism: the base layer is NCBI ``gene_info`` for the species
(symbol, GeneID, locus tag, MOD/Ensembl cross-references, synonyms). Human adds
the HGNC table; any species can add a UniProt ``idmapping_selected`` file.

Resolution priority (first unambiguous hit wins):
  1. official symbol / stable ids (GeneID, HGNC/MGI/SGD/... ids, Ensembl,
     locus tag, UniProt)
  2. previous symbol (only if it maps to one gene)
  3. alias / synonym (only if it maps to one gene)
Ambiguous previous/alias symbols are reported, never guessed - the same
"no introduced noise" principle used for pathways.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from .io import open_text


class GeneMapper:
    def __init__(self) -> None:
        self.exact: dict[str, str] = {}  # approved symbol and stable ids -> symbol
        self.prev: dict[str, set[str]] = {}
        self.alias: dict[str, set[str]] = {}
        self.entrez: dict[str, str] = {}  # symbol -> NCBI Gene id
        self.stats: Counter = Counter()

    def _add_exact(self, key: str, sym: str) -> None:
        self.exact.setdefault(key.strip().upper(), sym)

    @classmethod
    def from_ncbi_gene_info(cls, path: str | Path, taxid: str | None = None) -> "GeneMapper":
        """Base layer for any organism: NCBI ``<Species>.gene_info(.gz)``."""
        m = cls()
        with open_text(path) as fh:
            header = fh.readline().lstrip("#").rstrip("\n").split("\t")
            c = {h: i for i, h in enumerate(header)}
            rows = []
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if taxid and f[c["tax_id"]] != taxid:
                    continue
                rows.append(f)
        for f in rows:  # official symbols first so they win over ids/synonyms
            m.exact[f[c["Symbol"]].upper()] = f[c["Symbol"]]
            m.entrez[f[c["Symbol"]]] = f[c["GeneID"]]
        for f in rows:
            sym = f[c["Symbol"]]
            m._add_exact(f[c["GeneID"]], sym)
            if f[c["LocusTag"]] not in ("-", ""):
                m._add_exact(f[c["LocusTag"]], sym)
            for x in f[c["dbXrefs"]].split("|"):
                if ":" in x:
                    db, val = x.split(":", 1)
                    m._add_exact(val if db in ("Ensembl", "FLYBASE", "WormBase", "TAIR", "ZFIN") else x, sym)
                    if db in ("MGI", "HGNC", "RGD", "SGD"):  # "MGI:MGI:12345" -> also "MGI:12345"
                        m._add_exact(val, sym)
            nom = f[c.get("Symbol_from_nomenclature_authority", -1)] if "Symbol_from_nomenclature_authority" in c else "-"
            if nom not in ("-", "") and nom != sym:
                m.alias.setdefault(nom.upper(), set()).add(sym)
            for syn in f[c["Synonyms"]].split("|"):
                if syn not in ("-", ""):
                    m.alias.setdefault(syn.upper(), set()).add(sym)
        return m

    @classmethod
    def from_kegg_list(cls, path: str | Path) -> "GeneMapper":
        """KEGG ``list/<org>``: "eco:b0001<TAB>CDS<TAB>pos<TAB>thrL; thr operon leader peptide".

        Used where NCBI Gene has no records (most prokaryotes). The symbol is the
        first name before ';'; genes without one keep their KEGG id (locus tag).
        """
        m = cls()
        with open_text(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if not f or ":" not in f[0]:
                    continue
                locus = f[0].split(":", 1)[1]
                desc = f[-1] if len(f) > 1 else ""
                names = [n.strip() for n in desc.split(";", 1)[0].split(",")] if ";" in desc else []
                sym = names[0] if names and names[0] and " " not in names[0] else locus
                m.exact.setdefault(sym.upper(), sym)
                m.exact[locus.upper()] = sym
                for n in names[1:]:
                    if n and " " not in n:
                        m.alias.setdefault(n.upper(), set()).add(sym)
        return m

    def add_hgnc(self, path: str | Path) -> "GeneMapper":
        """Human layer: HGNC ``hgnc_complete_set.txt`` (approved entries only)."""
        with open_text(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(header)}

            def get(f: list[str], name: str) -> list[str]:
                i = col.get(name)
                if i is None or i >= len(f) or not f[i]:
                    return []
                return [x.strip().strip('"') for x in f[i].split("|") if x.strip()]

            for line in fh:
                f = line.rstrip("\n").split("\t")
                status = get(f, "status")
                if status and status[0] != "Approved":
                    continue
                sym = get(f, "symbol")[0]
                self.exact[sym.upper()] = sym  # HGNC approved symbol is authoritative
                for key in ("hgnc_id", "entrez_id", "ensembl_gene_id", "uniprot_ids"):
                    for v in get(f, key):
                        self.exact[v.upper()] = sym
                for v in get(f, "entrez_id"):
                    self.entrez[sym] = v
                for v in get(f, "prev_symbol"):
                    self.prev.setdefault(v.upper(), set()).add(sym)
                for v in get(f, "alias_symbol"):
                    self.alias.setdefault(v.upper(), set()).add(sym)
        return self

    @classmethod
    def from_hgnc(cls, path: str | Path) -> "GeneMapper":
        return cls().add_hgnc(path)

    def add_uniprot_idmapping(self, path: str | Path) -> "GeneMapper":
        """UniProt ``*_idmapping_selected.tab.gz``: accession (col 1) -> GeneID (col 3)."""
        by_geneid = {v: k for k, v in self.entrez.items()}
        with open_text(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 2 and f[2]:
                    sym = by_geneid.get(f[2].split(";")[0].strip())
                    if sym:
                        self._add_exact(f[0], sym)
                        self._add_exact(f[1], sym)  # UniProtKB entry name, e.g. P53_HUMAN
        return self

    @classmethod
    def for_organism(cls, org, raw_dir: str | Path) -> "GeneMapper":
        """Build the mapper from files fetched by :mod:`intpath.sources`."""
        from .sources import gene_namespace

        raw = Path(raw_dir)
        gi = gene_namespace(raw, org)
        m = cls.from_kegg_list(gi) if gi.name.endswith("kegg_genes.tsv") else cls.from_ncbi_gene_info(gi)
        if org.hgnc and (raw / "hgnc" / "hgnc_complete_set.txt").exists():
            m.add_hgnc(raw / "hgnc" / "hgnc_complete_set.txt")
        up = raw / "uniprot" / f"{org.key}_idmapping_selected.tab.gz"
        if up.exists():
            m.add_uniprot_idmapping(up)
        rest = raw / "uniprot" / f"{org.key}_uniprot_rest.tsv"
        if rest.exists():
            m.add_uniprot_rest(rest)
        return m

    def add_uniprot_rest(self, path: str | Path) -> "GeneMapper":
        """UniProt REST TSV (Entry, Entry Name, GeneID, ordered locus): accession -> symbol."""
        by_geneid = {v: k for k, v in self.entrez.items()}
        with open_text(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            col = {h: i for i, h in enumerate(header)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                sym = None
                for gid in f[col["GeneID"]].split(";") if "GeneID" in col else ():
                    sym = by_geneid.get(gid.strip())
                    if sym:
                        break
                if sym is None and "Gene Names (ordered locus)" in col:  # e.g. Rv0153c
                    for tag in f[col["Gene Names (ordered locus)"]].replace("/", " ").split():
                        sym = self.exact.get(tag.upper())
                        if sym:
                            break
                if sym:
                    self._add_exact(f[col["Entry"]], sym)
                    self._add_exact(f[col["Entry Name"]], sym)
        return self

    def map(self, ident: str) -> str | None:
        key = ident.strip().upper()
        if key.startswith("ENS") and "." in key:  # strip Ensembl version, any species
            key = key.split(".")[0]
        if key.startswith("UNIPROTKB:"):
            key = key[10:]
        if key.startswith("ENTREZ GENE/LOCUSLINK:") or key.startswith("ENTREZGENE/LOCUSLINK:"):
            key = key.split(":", 1)[1]
        if "-" in key and key.split("-")[0] in self.exact and key not in self.exact:  # UniProt isoform P12345-2
            key = key.split("-")[0]
        if key in self.exact:
            self.stats["exact"] += 1
            return self.exact[key]
        for table, label in ((self.prev, "previous"), (self.alias, "alias")):
            hits = table.get(key)
            if hits:
                if len(hits) == 1:
                    self.stats[label] += 1
                    return next(iter(hits))
                self.stats[f"ambiguous_{label}"] += 1
                return None
        self.stats["unmapped"] += 1
        return None

    def write_aliases(self, path: str | Path, genes: set[str]) -> int:
        """Unambiguous identifier/alias -> symbol table for genes in a release (used by the web UI)."""
        rows = self.alias_table(genes)
        with open(path, "w") as fh:
            fh.write("alias\tsymbol\n")
            for k in sorted(rows):
                fh.write(f"{k}\t{rows[k]}\n")
        return len(rows)

    def alias_table(self, genes: set[str]) -> dict[str, str]:
        rows: dict[str, str] = {}
        for table in (self.alias, self.prev):  # weakest first, so stronger keys overwrite
            for k, syms in table.items():
                if len(syms) == 1 and next(iter(syms)) in genes:
                    rows[k] = next(iter(syms))
        for k, sym in self.exact.items():
            if sym in genes:
                rows[k] = sym
        return rows

    def map_many(self, idents) -> dict[str, str]:
        out = {}
        for i in idents:
            s = self.map(i)
            if s:
                out[i] = s
        return out
