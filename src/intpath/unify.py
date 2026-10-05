"""Full unification of related pathways (BuildIntPathGEN / BuildIntPathGPR).

Every gene and every gene pair of every member pathway is kept ("no deletion,
no introduced noise"); each membership records which source(s) support it, and
each gene pair records the union of unified relation types.
"""

from __future__ import annotations

from typing import Callable, Iterable

from .model import GeneSet, SourcePathway, stable_id
from .curation import LEGACY_V2, rules_for
from .names import (
    V3_MISMATCHES,
    Match,
    PathwayGroup,
    clean_name,
    find_related_pairs,
    group_related,
    integrated_name,
)


def overlap_coefficient(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


# Sources whose pathways are curated as distinct entities with stable ids and an
# explicit hierarchy: similar names inside them are siblings, not duplicates.
NO_WITHIN_SOURCE_MERGE = frozenset({"Reactome", "GO"})


def merge_pathways(
    pathways: Iterable[SourcePathway],
    *,
    legacy: bool = False,
    min_jaccard: float = 0.1,
    extra_mismatches: Iterable[tuple[str, str]] = (),
    no_within: Iterable[str] = NO_WITHIN_SOURCE_MERGE,
    organism: str | None = None,
) -> tuple[list[GeneSet], list[Match]]:
    """Group related pathways by name and fully unify their genes and gene pairs.

    legacy=True  -> exact V2 rule set of ``organism`` (intpath.curation) on raw
                    names, no gene-overlap guard.
    legacy=False -> names cleaned; extended mismatch list; entity guard; no
                    within-source matching for hierarchical sources (Reactome);
                    and a name match is kept only if the gene sets have Jaccard
                    >= ``min_jaccard`` or the cleaned names are identical.
                    Jaccard (not the overlap coefficient) stops a superpathway
                    such as Reactome "Disease" absorbing its children.
    """
    by_key: dict[tuple[str, str], SourcePathway] = {}
    names: dict[str, list[str]] = {}
    for p in pathways:
        if p.key in by_key:  # same name twice in one source: union them
            q = by_key[p.key]
            q.genes |= p.genes
            for pr, rel in p.pairs.items():
                q.pairs.setdefault(pr, set()).update(rel)
            continue
        by_key[p.key] = p
        names.setdefault(p.source, []).append(p.name)

    rules = rules_for(organism, legacy)
    mismatches = tuple(rules.mismatches) + tuple(extra_mismatches)
    display: Callable[[str], str] | None = None
    overlap = None
    within = None
    if not legacy:
        mismatches += V3_MISMATCHES
        display = clean_name
        overlap = lambda a, b: jaccard(by_key[a].genes, by_key[b].genes)  # noqa: E731
        within = [s for s in names if s not in set(no_within)]

    matches = find_related_pairs(
        names,
        mismatches=mismatches,
        entity_guard=not legacy,
        within_sources=within,
        display=display,
        overlap=overlap,
        min_overlap=min_jaccard,
    )
    groups = group_related(
        matches, legacy=legacy, display=display, rules=rules, merge_same_name=legacy and organism == "musculus"
    )

    grouped = {k for g in groups for k in g.members}
    for key, p in by_key.items():  # every unmatched pathway is carried over unchanged
        if key not in grouped:
            nm = p.name if legacy else clean_name(p.name)
            groups.append(PathwayGroup(nm, [key]))

    return [unify_group(g, by_key) for g in groups], matches


def unify_group(group: PathwayGroup, by_key: dict[tuple[str, str], SourcePathway]) -> GeneSet:
    members = [(s, by_key[(s, n)].source_id or n) for s, n in group.members]
    gs = GeneSet(id=stable_id("IP", members), name=group.name.strip(), collection="pathway", members=members)
    for key in group.members:
        p = by_key[key]
        for g in p.genes:
            gs.genes.setdefault(g, set()).add(p.source)
        for pr, rels in p.pairs.items():
            entry = gs.pairs.setdefault(pr, {"rel": set(), "src": set()})
            entry["rel"].update(rels)
            entry["src"].add(p.source)
    return gs


__all__ = ["merge_pathways", "unify_group", "overlap_coefficient", "jaccard", "integrated_name", "LEGACY_V2"]
