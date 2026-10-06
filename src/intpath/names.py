"""Related-pathway-name identification (port of Integration.java, IntPath 2012).

Method (Zhou et al., BMC Syst Biol 2012, 6(Suppl 2):S2):

* alignment score  = length of the longest common subsequence (LCS) of the two
  names, compared case-insensitively;
* alignment ratio  = 2 * score / (len(a) + len(b));
* for every name in list x, only the *best hit* in list y (highest ratio, first
  one wins on ties) is considered;
* a best hit is accepted when
      (score > len(shorter) - 1 and ratio >= 0.5)  or  ratio > 0.91
  and the pair is not on the curated mismatch list (e.g. "T cell" vs "B cell");
* names are compared between every pair of sources *and* within each source;
* accepted pairs are grouped with a disjoint set (union-find); every connected
  component becomes one integrated pathway, named after its shortest member
  with roman-numeral style suffixes removed.

``legacy=True`` reproduces the 2012/2021 behaviour exactly on raw names.
``legacy=False`` (the IntPathV2 default) cleans names first (HTML tags, species
suffixes such as " - Homo sapiens (human)") and can require gene-set overlap as
a second line of evidence (see :func:`find_related_pairs`).
"""

from __future__ import annotations

import html
import re
from dataclasses import dataclass, field
from typing import Callable, Iterable, Sequence

# Curated old IntPath rules live in intpath.curation (one set per old IntPath organism).
from .curation import LEGACY_RULES, Rules  # noqa: E402

LEGACY_MISMATCHES = LEGACY_RULES["sapiens"].mismatches

# Additional IntPathV2 guards, found while extending to Reactome/GO vocabularies.
INTPATHV2_MISMATCHES: tuple[tuple[str, str], ...] = (
    ("VEGF", "EGF"),
    ("EGFR", "VEGFR"),
    ("Type I ", "Type II "),
    ("type I ", "type II "),
    ("MAPK", "MAP2K"),
    ("IL-1", "IL-17"),
    ("IL-2", "IL-12"),
    ("IL-3", "IL-13"),
    ("IL-4", "IL-14"),
    ("IL-6", "IL-16"),
    ("mTORC1", "mTORC2"),
    ("Class I ", "Class II "),
    ("biosynthesis", "degradation"),
    ("anabolism", "catabolism"),
    ("synthesis", "degradation"),
    ("activation", "inhibition"),
    ("positive regulation", "negative regulation"),
    # IntPathV2 curation (human KEGG + Reactome + WikiPathways review, 2026-10-06): pathways
    # whose names align but whose meanings differ are never merged
    ("Non-alcoholic", "Alcoholic liver"),
    ("Nonalcoholic", "Alcoholic liver"),
    ("Hepatitis B", "Hepatitis C"),
    ("Prion", "Parkinson"),
    ("EPO receptor", "Hepatocyte growth factor receptor"),
    ("Rap1", "Ras signaling"),
    ("Renin secretion", "Insulin secretion"),
    ("Osteoblast", "Osteoclast"),
    ("Omega 3", "Omega 6"),
    ("Carbohydrate digestion", "Fat digestion"),
    ("Non-small cell", "Small cell lung"),
    ("Non small cell", "Small cell lung"),
    ("Integrated breast cancer", "Integrated cancer pathway"),
    ("Focal adhesion PI3K", "PI3K Akt signaling"),
    ("Focal adhesion PI3K", "PI3K-Akt signaling"),
    ("Vitamins", "Vitamin D"),
    ("Neolacto series", "Lacto series"),
    ("biosynthesis - ganglio series", "Glycosphingolipid biosynthesis"),
    ("biosynthesis - globo and isoglobo series", "Glycosphingolipid biosynthesis"),
    ("biosynthesis - lacto and neolacto series", "Glycosphingolipid biosynthesis"),
    ("Alternative complement", "Classical antibody-mediated complement"),
    ("Dengue 2 interactions", "Complement and coagulation cascades"),
    ("Hippo Merlin signaling dysregulation", "Hippo signaling regulation"),
)

# --------------------------------------------------------------------------- #
# Similarity
# --------------------------------------------------------------------------- #
def lcs_length(a: str, b: str) -> int:
    """Length of the longest common subsequence (bit-parallel, Hyyro 2004)."""
    if not a or not b:
        return 0
    if len(a) < len(b):
        a, b = b, a
    masks: dict[str, int] = {}
    for i, ch in enumerate(a):
        masks[ch] = masks.get(ch, 0) | (1 << i)
    full = (1 << len(a)) - 1
    v = full
    for ch in b:
        u = v & masks.get(ch, 0)
        v = ((v + u) | (v - u)) & full
    return len(a) - bin(v).count("1")


@dataclass(frozen=True)
class Alignment:
    score: int
    ratio: float


def align(a: str, b: str) -> Alignment:
    """IntPath alignment of two names (case-insensitive LCS)."""
    a, b = a.strip(), b.strip()
    score = lcs_length(a.lower(), b.lower())
    total = len(a) + len(b)
    return Alignment(score, 2.0 * score / total if total else 0.0)


def is_related(a: str, b: str, aln: Alignment) -> bool:
    """The empirically determined acceptance rule of the 2012 paper."""
    shorter = min(len(a.strip()), len(b.strip()))
    return (aln.score > shorter - 1 or aln.ratio > 0.91) and aln.ratio >= 0.5


def is_mismatch(a: str, b: str, mismatches: Iterable[tuple[str, str]]) -> bool:
    for x, y in mismatches:
        if (x in a and y in b) or (x in b and y in a):
            return True
    return False


_TOKEN = re.compile(r"[A-Za-z0-9]+")
_ROMAN = re.compile(r"^[IVX]+$")
_SUBTYPE_WORDS = frozenset({"type", "class", "group", "complex", "phase", "stage", "grade"})


def _roman_to_int(t: str) -> int:
    if t.isdigit():
        return int(t)
    vals = {"I": 1, "V": 5, "X": 10}
    total = 0
    for i, c in enumerate(t):
        v = vals[c]
        total += -v if i + 1 < len(t) and vals[t[i + 1]] > v else v
    return total


def _entities(name: str) -> set[str]:
    """Gene/protein-like tokens: >= 2 capitals (not a roman numeral) or letters+digits.

    A bare number directly after such a token is joined to it ("IL-1" -> IL1,
    "TGF-beta" -> TGF), so hyphenation differences do not matter. A number or
    roman numeral after a subtype word is an entity too ("type 2" == "Type II").
    """
    toks = _TOKEN.findall(name)
    out: set[str] = set()
    prev = None
    for i, t in enumerate(toks):
        if i and toks[i - 1].lower() in _SUBTYPE_WORDS and (t.isdigit() or _ROMAN.match(t)):
            out.add(f"{toks[i - 1].upper()}{_roman_to_int(t)}")  # "type 2" == "Type II"
            prev = None
            continue
        if t.isdigit() and prev is not None:
            out.discard(prev)
            out.add(prev + t)
            prev = None
            continue
        caps = sum(c.isupper() for c in t)
        is_ent = (caps >= 2 and not _ROMAN.match(t)) or (any(c.isalpha() for c in t) and any(c.isdigit() for c in t))
        prev = t.upper() if is_ent else None
        if is_ent:
            out.add(prev)
    return out


def numbered_entity_mismatch(a: str, b: str) -> bool:
    """IntPathV2 guard: the two names are about different molecular entities.

    "Signaling by FGFR1" vs "Signaling by FGFR2", "IL-1 signaling" vs "IL-17
    signaling", "RHOC GTPase cycle" vs "RHOG GTPase cycle": each name carries a
    gene-like token the other lacks. Free-standing roman/arabic variant suffixes
    ("... I" vs "... II", the BioCyc variant convention merged on purpose in
    2012) are not entities, so those merges are unaffected.
    """
    ea, eb = _entities(a), _entities(b)
    return bool(ea - eb) and bool(eb - ea)


# --------------------------------------------------------------------------- #
# Name cleaning (IntPathV2)
# --------------------------------------------------------------------------- #
_TAG = re.compile(r"<[^>]+>")
# KEGG REST: " - Homo sapiens (human)", " - Mycobacterium tuberculosis H37Rv" (strain tokens carry
# a capital or a digit, so pathway subtitles such as " - multiple species" are kept);
# disease subtitles ("MPS I - Hurler syndrome (...)") are not species names
_SPECIES_SUFFIX = re.compile(r"\s+-\s+[A-Z][a-z]+ (?!syndrome\b|disease\b)[a-z]+(?:\s+(?=\S*[A-Z0-9])\S+)*(\s+\([^)]*\))?\s*$")
_WS = re.compile(r"\s+")


def clean_name(name: str) -> str:
    """Strip markup and source-specific decorations; keeps case and wording."""
    name = html.unescape(_TAG.sub("", name))
    name = _SPECIES_SUFFIX.sub("", name)
    return _WS.sub(" ", name).strip()


def integrated_name(shortest: str, legacy: bool = True, rules: Rules | None = None) -> str:
    """Derive the integrated pathway name from the shortest member name (processIntPathNames)."""
    rules = rules or LEGACY_RULES["sapiens"]
    drop = {t.lower() for t in rules.drop_tokens} if rules.drop_case_insensitive else rules.drop_tokens
    tokens = []
    for tok in shortest.split(" "):
        if (tok.lower() if rules.drop_case_insensitive else tok) in drop:
            continue
        tokens.append(rules.replace.get(tok.lower(), tok))
    name = " ".join(tokens)
    if legacy:
        return name + " "  # the Java implementation leaves a trailing blank
    # IntPathV2: keep every word of the shortest name (the old per-organism drop lists removed
    # words such as "Small" and broke names like "Small cell lung cancer"); drop only a trailing
    # variant numeral ("... biosynthesis II" -> "... biosynthesis")
    words = shortest.strip().split()
    while len(words) > 1 and _ROMAN.match(words[-1]):
        words.pop()
    return " ".join(words) or shortest.strip()


# --------------------------------------------------------------------------- #
# Best-hit search and grouping
# --------------------------------------------------------------------------- #
PathwayKey = tuple[str, str]  # (source, pathway name)


@dataclass
class Match:
    a: PathwayKey
    b: PathwayKey
    score: int
    ratio: float
    overlap: float | None = None  # gene-set overlap coefficient, when checked
    decision: str = "accept"  # accept | reject | pending (IntPathV2 review, see intpath.review)
    reason: str = ""


def _best_hit(name: str, candidates: Sequence[str]) -> tuple[str | None, Alignment]:
    best, best_aln = None, Alignment(0, 0.0)
    n = len(name.strip())
    for cand in candidates:
        m = len(cand.strip())
        # ratio can never exceed 2*min/(n+m); skip hopeless candidates early
        if n + m == 0 or 2.0 * min(n, m) / (n + m) <= best_aln.ratio:
            continue
        aln = align(name, cand)
        if aln.ratio > best_aln.ratio:
            best, best_aln = cand, aln
    return best, best_aln


def find_related_pairs(
    names: dict[str, Sequence[str]],
    *,
    mismatches: Iterable[tuple[str, str]] = LEGACY_MISMATCHES,
    entity_guard: bool = False,
    within_sources: Iterable[str] | None = None,
    display: Callable[[str], str] | None = None,
    overlap: Callable[[PathwayKey, PathwayKey], float] | None = None,
    min_overlap: float = 0.0,
) -> list[Match]:
    """Pairwise best-hit comparison between and within sources.

    ``names`` maps a source label to its ordered list of pathway names.
    ``within_sources`` limits the within-source comparison (default: all).
    ``entity_guard`` enables :func:`numbered_entity_mismatch` (IntPathV2).
    ``display`` optionally maps a raw name to the string that is aligned
    (e.g. :func:`clean_name`). If ``overlap`` is given, an accepted name match
    is kept only when ``overlap(a, b) >= min_overlap`` *or* the names are
    identical after cleaning - a guard against look-alike names.
    """
    mismatches = tuple(mismatches)
    show = display or (lambda s: s)
    sources = list(names)
    shown = {s: [show(n) for n in names[s]] for s in sources}
    back = {s: dict(zip(shown[s], names[s])) for s in sources}
    out: list[Match] = []

    def consider(sa: str, na: str, sb: str, cands: Sequence[str]) -> None:
        hit, aln = _best_hit(na, cands)
        if hit is None or not is_related(na, hit, aln) or is_mismatch(na, hit, mismatches):
            return
        if entity_guard and numbered_entity_mismatch(na, hit):
            return
        ka, kb = (sa, back[sa][na]), (sb, back[sb][hit])
        ov = None
        if overlap is not None:
            ov = overlap(ka, kb)
            if ov < min_overlap and na.strip().lower() != hit.strip().lower():
                return
        out.append(Match(ka, kb, aln.score, aln.ratio, ov))

    # between sources, in the legacy order K-W, K-C, C-W generalised to i<j
    for i, sa in enumerate(sources):
        for sb in sources[i + 1 :]:
            for na in shown[sa]:
                consider(sa, na, sb, shown[sb])
    # within each source: compare each name to the names after it
    within = set(sources if within_sources is None else within_sources)
    for s in sources:
        if s not in within:
            continue
        lst = shown[s]
        for k, na in enumerate(lst):
            consider(s, na, s, lst[k + 1 :])
    return out


def java_string_hash(s: str) -> int:
    """java.lang.String.hashCode (UTF-16 code units, 32-bit overflow)."""
    h = 0
    for unit in s.encode("utf-16-be").hex(" ", 2).split():
        h = (31 * h + int(unit, 16)) & 0xFFFFFFFF
    return h


def java_hashmap_order(keys_in_insertion_order: Sequence[str]) -> list[str]:
    """Iteration order of a JDK <= 7 java.util.HashMap<String, ?> holding these keys.

    Old IntPath picked integrated pathway names from HashMap iteration order when two
    candidate names had equal length; emulating it makes the legacy rebuild
    reproduce old IntPath labels. The JDK 7 supplemental hash matches the member order
    printed in the old IntPath RelPthNamsGEN files (human 57/57 groups, M. tuberculosis
    35/35, yeast 76/76, mouse 84/85); the JDK 8 hash does not.
    Order = bucket index, then insertion order.
    """
    keys = list(dict.fromkeys(keys_in_insertion_order))
    cap = 16
    while len(keys) > cap * 0.75:
        cap *= 2
    def bucket(k: str) -> int:
        h = java_string_hash(k)
        h ^= (h >> 20) ^ (h >> 12)
        h = (h ^ (h >> 7) ^ (h >> 4)) & 0xFFFFFFFF
        return h & (cap - 1)
    return sorted(keys, key=lambda k: bucket(k))  # stable sort keeps insertion order within a bucket


class DisjointSet:
    def __init__(self) -> None:
        self.parent: dict = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:  # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb

    def groups(self) -> list[list]:
        comp: dict = {}
        for x in self.parent:
            comp.setdefault(self.find(x), []).append(x)
        return list(comp.values())


@dataclass
class PathwayGroup:
    name: str
    members: list[PathwayKey] = field(default_factory=list)


def group_related(
    matches: Iterable[Match],
    legacy: bool = True,
    display: Callable[[str], str] | None = None,
    rules: Rules | None = None,
    merge_same_name: bool = False,
) -> list[PathwayGroup]:
    """Union-find over accepted matches; name each component by its shortest member.

    ``merge_same_name`` joins components that end up with the same integrated
    name (old IntPath mouse behaviour); otherwise they stay separate sets.
    """
    show = display or (lambda s: s)
    matches = list(matches)
    ds = DisjointSet()
    for m in matches:
        ds.union(m.a, m.b)
    java_rank: dict = {}
    if legacy:  # Old IntPath node strings were "<name>+<K|C|W>" keys of a HashMap
        node = lambda k: f"{k[1]}+{k[0]}"  # noqa: E731
        order = java_hashmap_order([node(k) for m in matches for k in (m.a, m.b)])
        java_rank = {n: i for i, n in enumerate(order)}
    by_name: dict[str, PathwayGroup] = {}
    groups: list[PathwayGroup] = []
    for members in ds.groups():
        members = sorted(members, key=lambda k: (k[0], k[1]))
        if legacy:  # first strictly-shortest in Java iteration order
            shortest = min(members, key=lambda k: (len(k[1]), java_rank[node(k)]))
        else:
            shortest = min(members, key=lambda k: (len(show(k[1])), k))
        name = integrated_name(show(shortest[1]), legacy=legacy, rules=rules)
        if merge_same_name and name in by_name:
            by_name[name].members = sorted(by_name[name].members + members)
            continue
        g = PathwayGroup(name, members)
        by_name[name] = g
        groups.append(g)
    return sorted(groups, key=lambda g: g.name.lower())
