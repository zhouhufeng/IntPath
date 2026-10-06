"""Review of candidate pathway merges (replaces the manual review of the 2012 IntPath).

The IntPath algorithm proposes related pathway pairs by name alignment (best hit,
acceptance rule). In 2012 the false pairs were removed by hand (the curated
mismatch list). In IntPathV2 every proposed pair is reviewed by an LLM for whether
the two names describe the same pathway; decisions and reasons are kept in
``data/merge_curation.tsv`` (versioned, reusable across organisms and releases):

* identical names (after clean_name)  -> accept (same pathway in two sources)
* reviewed "merge"                    -> accept
* reviewed "separate"                 -> reject (logged as a prevented wrong merge)
* not yet reviewed                    -> pending: not merged; listed for review

Name review: the reviewed "merge" pairs also apply when the name alignment never
proposed them (synonyms such as "Citrate cycle (TCA cycle)" / "Citric acid cycle"):
``name_review_matches`` adds them as accepted matches, so pathways merge purely on
whether their names mean the same biological process.

The review procedure and criteria are documented in Docs/Curation/LLM_REVIEW.md.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from .names import Match, clean_name

CURATION = Path(__file__).parent / "data" / "merge_curation.tsv"


@lru_cache(maxsize=4)
def load_curation(path: str | Path = CURATION) -> dict[frozenset, tuple[str, str]]:
    out: dict[frozenset, tuple[str, str]] = {}
    with open(path) as fh:
        rows = csv.reader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for r in rows:
            if len(r) >= 3 and r[0] != "name_a":
                out[frozenset((r[0].strip().lower(), r[1].strip().lower()))] = (r[2].strip(), r[3].strip() if len(r) > 3 else "")
    return out


def review(matches: list[Match], curation: dict[frozenset, tuple[str, str]] | None = None) -> list[Match]:
    """Set decision/reason on every candidate match (in place) and return the list."""
    cur = load_curation() if curation is None else curation
    for m in matches:
        a, b = clean_name(m.a[1]), clean_name(m.b[1])
        if a.lower() == b.lower():
            m.decision, m.reason = "accept", "identical names"
            continue
        hit = cur.get(frozenset((a.lower(), b.lower())))
        if hit is None:
            m.decision, m.reason = "pending", "not yet reviewed"
        elif hit[0] == "merge":
            m.decision, m.reason = "accept", "reviewed: " + hit[1]
        else:
            m.decision, m.reason = "reject", "reviewed: " + hit[1]
    return matches


def name_review_matches(names: dict[str, list[str]], matches: list[Match],
                        curation: dict[frozenset, tuple[str, str]] | None = None) -> list[Match]:
    """Accepted matches for reviewed "merge" pairs not already proposed by the alignment."""
    cur = load_curation() if curation is None else curation
    keys: dict[str, list[tuple[str, str]]] = {}
    for src, nms in names.items():
        for n in nms:
            keys.setdefault(clean_name(n).lower(), []).append((src, n))
    seen = {frozenset((m.a, m.b)) for m in matches}
    out: list[Match] = []
    for pair, (decision, reason) in cur.items():
        if decision != "merge" or len(pair) != 2:
            continue
        a, b = tuple(pair)
        for ka in keys.get(a, ()):
            for kb in keys.get(b, ()):
                if frozenset((ka, kb)) not in seen:
                    seen.add(frozenset((ka, kb)))
                    out.append(Match(ka, kb, 0, 0.0, decision="accept", reason="name review: " + reason))
    return out


def summary(matches: list[Match]) -> dict:
    from collections import Counter

    c = Counter(m.decision for m in matches)
    return {"candidates": len(matches), "accepted": c.get("accept", 0), "rejected": c.get("reject", 0),
            "pending_review": c.get("pending", 0),
            "accepted_identical_names": sum(1 for m in matches if m.reason == "identical names"),
            "accepted_by_name_review_only": sum(1 for m in matches if m.reason.startswith("name review"))}


def write_log(matches: list[Match], path: str | Path) -> None:
    with open(path, "w") as fh:
        fh.write("decision\tsource_a\tpathway_a\tsource_b\tpathway_b\talign_ratio\treason\n")
        for m in sorted(matches, key=lambda m: (m.decision, clean_name(m.a[1]).lower())):
            fh.write(f"{m.decision}\t{m.a[0]}\t{clean_name(m.a[1])}\t{m.b[0]}\t{clean_name(m.b[1])}\t{m.ratio:.3f}\t{m.reason}\n")
