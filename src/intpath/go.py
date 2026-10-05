"""Gene Ontology integration (V3): GO terms as IntPath gene sets, merged the IntPath way.

1. Parse the ontology (go-basic.obo) and the organism's GAF.
2. Propagate annotations up is_a / part_of (true-path rule); NOT-qualified and
   ND annotations are dropped; evidence codes are kept as provenance.
3. *Lossless redundancy merge*: GO terms of one namespace whose propagated gene
   sets are identical and that sit on one ancestor/descendant chain are merged
   into one gene set (named after the most specific term; all term ids kept as
   members). No gene or term is dropped - the full-unification principle.
4. *Cross-linking with pathways*: a GO set and an integrated pathway are linked
   when their gene sets are near-identical (Jaccard >= 0.5) or their names pass
   the IntPath alignment rule and genes overlap (Jaccard >= 0.2). GO and
   pathways stay separate sets (different concept types) but carry ``links``,
   which the enrichment report uses to collapse redundant hits.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .io import open_text
from .mapping import GeneMapper
from .model import GeneSet
from .names import align, clean_name, is_related
from .unify import jaccard

NAMESPACE = {"biological_process": "GO:BP", "molecular_function": "GO:MF", "cellular_component": "GO:CC"}
ASPECT = {"P": "GO:BP", "F": "GO:MF", "C": "GO:CC"}


@dataclass
class Term:
    id: str
    name: str = ""
    namespace: str = ""
    parents: set[str] = field(default_factory=set)
    alt_ids: list[str] = field(default_factory=list)
    obsolete: bool = False
    replaced_by: str = ""


def parse_obo(path: str | Path) -> dict[str, Term]:
    terms: dict[str, Term] = {}
    cur: Term | None = None
    in_term = False
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                in_term = line == "[Term]"
                cur = None
                continue
            if not in_term or ":" not in line:
                continue
            tag, val = line.split(":", 1)
            val = val.strip()
            if tag == "id":
                cur = terms.setdefault(val, Term(val))
            elif cur is None:
                continue
            elif tag == "name":
                cur.name = val
            elif tag == "namespace":
                cur.namespace = val
            elif tag == "is_a":
                cur.parents.add(val.split()[0])
            elif tag == "relationship" and val.startswith("part_of "):
                cur.parents.add(val.split()[1])
            elif tag == "alt_id":
                cur.alt_ids.append(val)
            elif tag == "is_obsolete" and val == "true":
                cur.obsolete = True
            elif tag == "replaced_by":
                cur.replaced_by = val
    return terms


def ancestors(terms: dict[str, Term]) -> dict[str, set[str]]:
    memo: dict[str, set[str]] = {}

    def anc(t: str) -> set[str]:
        if t in memo:
            return memo[t]
        memo[t] = set()  # cycle guard
        out: set[str] = set()
        for p in terms[t].parents if t in terms else ():
            out.add(p)
            out |= anc(p)
        memo[t] = out
        return out

    for t in terms:
        anc(t)
    return memo


def parse_gaf(
    path: str | Path,
    terms: dict[str, Term],
    mapper: GeneMapper | None = None,
    exclude_evidence: frozenset[str] = frozenset({"ND"}),
    taxid: str | None = None,
) -> dict[str, dict[str, set[str]]]:
    """Direct annotations: term -> gene symbol -> evidence codes."""
    alt = {a: t.id for t in terms.values() for a in t.alt_ids}
    ann: dict[str, dict[str, set[str]]] = defaultdict(dict)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 15 or "NOT" in f[3] or f[6] in exclude_evidence:
                continue
            if taxid and f"taxon:{taxid}" not in f[12].split("|"):
                continue
            go = alt.get(f[4], f[4])
            t = terms.get(go)
            if t is None or t.obsolete:
                if t is not None and t.replaced_by:
                    go = t.replaced_by
                else:
                    continue
            gene = f[2]
            if mapper is not None:
                gene = mapper.map(f[2]) or mapper.map(f[1]) or ""
            if gene:
                ann[go].setdefault(gene, set()).add(f[6])
    return ann


def propagate(ann: dict[str, dict[str, set[str]]], terms: dict[str, Term]) -> dict[str, dict[str, set[str]]]:
    anc = ancestors(terms)
    full: dict[str, dict[str, set[str]]] = defaultdict(dict)
    for t, genes in ann.items():
        for a in {t} | anc.get(t, set()):
            dst = full[a]
            for g, ev in genes.items():
                dst.setdefault(g, set()).update(ev)
    return full


def build_go_sets(
    terms: dict[str, Term],
    full: dict[str, dict[str, set[str]]],
    min_size: int = 5,
    max_size: int = 2000,
    merge_identical: bool = True,
) -> list[GeneSet]:
    anc = ancestors(terms)
    by_content: dict[tuple, list[str]] = defaultdict(list)
    for t, genes in full.items():
        term = terms.get(t)
        if term is None or term.namespace not in NAMESPACE or not (min_size <= len(genes) <= max_size):
            continue
        key = (term.namespace, frozenset(genes)) if merge_identical else (t,)
        by_content[key].append(t)

    sets: list[GeneSet] = []
    for ids in by_content.values():
        # identical gene sets on one lineage merge; unrelated identical terms stay apart
        ds_groups: list[list[str]] = []
        for t in sorted(ids):
            for grp in ds_groups:
                if any(t in anc.get(o, ()) or o in anc.get(t, ()) for o in grp):
                    grp.append(t)
                    break
            else:
                ds_groups.append([t])
        for grp in ds_groups:
            specific = [t for t in grp if not any(t in anc.get(o, ()) for o in grp if o != t)]
            rep = sorted(specific)[0]
            genes = full[rep]
            sets.append(
                GeneSet(
                    id=rep,
                    name=terms[rep].name,
                    collection=NAMESPACE[terms[rep].namespace],
                    genes={g: {"GO:" + ",".join(sorted(ev))} for g, ev in genes.items()},
                    members=[("GO", t) for t in sorted(grp)],
                )
            )
    return sets


def link_sets(
    a_sets: list[GeneSet], b_sets: list[GeneSet], min_jaccard_genes: float = 0.5, min_jaccard_named: float = 0.2
) -> int:
    """Cross-link two collections (e.g. GO vs pathways). Returns number of links."""
    index: dict[str, list[int]] = defaultdict(list)
    for j, b in enumerate(b_sets):
        for g in b.genes:
            index[g].append(j)
    n = 0
    for a in a_sets:
        shared: dict[int, int] = defaultdict(int)
        for g in a.genes:
            for j in index.get(g, ()):
                shared[j] += 1
        for j, k in shared.items():
            b = b_sets[j]
            jac = k / (len(a.genes) + len(b.genes) - k)
            if jac < min_jaccard_named:
                continue
            ok = jac >= min_jaccard_genes
            if not ok:
                na, nb = clean_name(a.name), clean_name(b.name)
                ok = is_related(na, nb, align(na, nb))
            if ok:
                a.links.append(b.id)
                b.links.append(a.id)
                n += 1
    return n


__all__ = ["parse_obo", "parse_gaf", "propagate", "build_go_sets", "link_sets", "jaccard"]
