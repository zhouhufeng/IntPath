"""Curated merge rules (from the old IntPath per-organism Integration.java files).

Old IntPath kept a hand-curated rule set per organism, written while reviewing each
organism's merge results:

* ``mismatches``  - pairs of name fragments that look alike but name different
                    pathways (Integration.ckmismatch)
* ``drop_tokens`` - tokens removed when naming an integrated pathway (ckNameParts),
                    compared case-insensitively when ``drop_case_insensitive``
* ``replace``     - whole-token replacements in the integrated name (SpecialNameReplace)

``rules_for(organism)`` returns the old IntPath rules for the four old IntPath organisms. New
organisms, and every IntPathV2 build, use ``INTPATHV2_RULES``: the union of all old IntPath mismatch
lists plus the IntPathV2 additions in :mod:`intpath.names`. Add new curation here.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Rules:
    mismatches: tuple[tuple[str, str], ...]
    drop_tokens: frozenset[str] = field(default_factory=frozenset)
    drop_case_insensitive: bool = False
    replace: dict[str, str] = field(default_factory=dict)  # lower-case token -> replacement


LEGACY_RULES: dict[str, Rules] = {
    "sapiens": Rules(
        mismatches=(
            ('NOD', 'Toll'),
            ('Linoleic', 'Lipoic'),
            ('T cell', 'B cell'),
            ('EPO', 'TOR'),
            ('L-cysteine', 'lysine'),
            ('spermine', 'serine'),
            ('serotonin', 'serine'),
            ('Steroid', 'thyroid'),
            ('isoleucine', 'leucine'),
            ('dermatan', 'heparan'),
            ('ribonucleotides', 'deoxyribonucleotides'),
            ('guanosine', 'adenosine'),
            ('sulfation', 'oxidation'),
            ('Glycerolipid', 'Glycerophospholipid'),
        ),
        drop_tokens=frozenset(['', '-', 'I', 'II', 'III', 'IV', 'Small', 'Type', 'V', 'VI', 'X', 'adenosine', 'globo', 'heparan', 'keratan']),
        drop_case_insensitive=False,
        replace={'glycolysis': 'Glycolysis and Gluconeogenesis'},
    ),
    "musculus": Rules(
        mismatches=(
            ('NOD', 'Toll'),
            ('Linoleic', 'Lipoic'),
            ('T cell', 'B cell'),
            ('L-cysteine', 'lysine'),
            ('spermine', 'serine'),
            ('serotonin', 'serine'),
            ('Steroid', 'thyroid'),
            ('dermatan', 'heparan'),
            ('ribonucleotides', 'deoxyribonucleotides'),
            ('guanosine', 'adenosine'),
            ('phenylalanine', 'alanine'),
            ('serine', 'Steroid'),
            ('T Cell', 'B Cell'),
            ('lactose', 'galactose'),
            ('galactosamine', 'glucosamine'),
            ('purine', 'pyrimidine'),
            ('deoxy', 'ribose'),
            ('GMP', 'AMP'),
            ('VEGF', 'EGFR1'),
            ('methylglyoxal', 'methylation'),
            ('Glycerolipid', 'Glycerophospholipid'),
        ),
        drop_tokens=frozenset(['', '(from', '(unsaturated,', 'ADP-D-Glucose)', 'I', 'II', 'III', 'IV', 'Small', 'Type', 'V', 'VI', 'X', 'adenosine', 'branch', 'even', 'globo', 'heparan', 'keratan', 'number)', 'pathway']),
        drop_case_insensitive=False,
        replace={'o-glycan': 'Glycan', 'il-2': 'Interleukin', 'ubiquinone-8': 'ubiquinone', 'glycolysis': 'Glycolysis and Gluconeogenesis', '3-hydroxypropionate': '3-/4-hydroxypropionate', 'leucine': 'Valine, leucine and isoleucine'},
    ),
    "cerevisiae": Rules(
        mismatches=(
            ('NOD', 'Toll'),
            ('Linoleic', 'Lipoic'),
            ('T cell', 'B cell'),
            ('EPO', 'TOR'),
            ('L-cysteine', 'lysine'),
            ('spermine', 'serine'),
            ('serotonin', 'serine'),
            ('Steroid', 'thyroid'),
            ('dermatan', 'heparan'),
            ('ribonucleotides', 'deoxyribonucleotides'),
            ('guanosine', 'adenosine'),
            ('sulfation', 'oxidation'),
            ('Lysine', 'Glycine'),
            ('Serine', 'homoserine'),
            ('galactose', 'Lactose'),
            ('alanine biosynthesis', 'Valine Biosynthesis'),
            ('valine biosynthesis', 'alanine biosynthesis'),
            ('oleate', 'folate'),
            ('glutaredoxin', 'thioredoxin'),
            ('homoserine', 'homocysteine'),
            ('valine degradation', 'alanine degradation'),
            ('phenylalanine biosynthesis', 'alanine biosynthesis'),
            ('glutathione', 'glutamine'),
            ('S-adenosylmethionine', 'Methionine'),
            ('Glycerolipid', 'Glycerophospholipid'),
            ('guanine, xanthine', 'adenine, hypoxanthine'),
        ),
        drop_tokens=frozenset(['', ',', '-', 'I', 'II', 'III', 'IV', 'Small', 'Type', 'V', 'VI', 'X', 'adenosine', 'globo', 'heparan', 'keratan', 'saturated']),
        drop_case_insensitive=True,
        replace={'o-glycan': 'Glycan', 'il-2': 'Interleukin', 'ubiquinone-8': 'ubiquinone', 'palmitate': 'palmitoleate and palmitate', 'glycolysis': 'Glycolysis and Gluconeogenesis', 'deoxyribose': 'Ribose and Deoxyribose', 'leucine': 'Valine, leucine and isoleucine', 'spermine': 'spermine, spermidine'},
    ),
    "tuberculosis": Rules(
        mismatches=(
            ('NOD', 'Toll'),
            ('Linoleic', 'Lipoic'),
            ('T cell', 'B cell'),
            ('EPO', 'TOR'),
            ('L-cysteine', 'lysine'),
            ('spermine', 'serine'),
            ('serotonin', 'serine'),
            ('Steroid', 'thyroid'),
            ('spermidine', 'serine'),
            ('dermatan', 'heparan'),
            ('ribonucleotides', 'deoxyribonucleotides'),
            ('guanosine', 'adenosine'),
            ('sulfation', 'oxidation'),
            ('Lysine', 'Glycine'),
            ('glycine biosynthesis I', 'lysine biosynthesis I'),
            ('Serine', 'homoserine'),
            ('galactose', 'Lactose'),
            ('Nitrotoluene', 'Toluene'),
            ('alanine biosynthesis', 'Valine Biosynthesis'),
            ('Fluorobenzoate degradation', 'Benzoate'),
            ('valine biosynthesis', 'alanine biosynthesis'),
            ('oleate', 'folate'),
            ('glutaredoxin', 'thioredoxin'),
            ('homoserine', 'homocysteine'),
            ('valine degradation', 'alanine degradation'),
            ('Fluorobenzoate', 'Benzoate'),
            ('enterobactin', 'antigen'),
            ('glutathione', 'glutamine'),
            ('S-adenosylmethionine', 'Methionine'),
            ('Aminobenzoate degradation', 'Benzoate degradation'),
            ('deoxy', 'ribose degradation'),
            ('Glycerolipid', 'Glycerophospholipid'),
            ('demethylmenaquinone-8 biosynthesis I', 'menaquinone-8 biosynthesis'),
        ),
        drop_tokens=frozenset(['', ',', '-', 'I', 'II', 'III', 'IV', 'Small', 'Type', 'V', 'VI', 'X', 'adenosine', 'globo', 'heparan', 'keratan', 'saturated']),
        drop_case_insensitive=True,
        replace={'o-glycan': 'Glycan', 'il-2': 'Interleukin', 'glycolysis': 'Glycolysis and Gluconeogenesis', '4-aminobutyrate': 'arginine, putrescine, and 4-aminobutyrate', 'alanine': 'alanine and phenylalanine', 'leucine': 'leucine and isoleucine', 'ornithine': 'arginine and ornithine', 'spermine': 'spermine, spermidine', 'folate': 'Folate and tetrahydrofolate'},
    ),
}



def _union_mismatches() -> tuple[tuple[str, str], ...]:
    seen: dict[tuple[str, str], None] = {}
    for rules in LEGACY_RULES.values():
        for p in rules.mismatches:
            seen.setdefault(p, None)
    return tuple(seen)


ALL_LEGACY_MISMATCHES = _union_mismatches()
INTPATHV2_RULES = Rules(
    mismatches=ALL_LEGACY_MISMATCHES,
    drop_tokens=LEGACY_RULES["sapiens"].drop_tokens,
    replace={"glycolysis": "Glycolysis and Gluconeogenesis"},
)


def rules_for(organism: str | None, legacy: bool) -> Rules:
    if legacy and organism in LEGACY_RULES:
        return LEGACY_RULES[organism]
    if legacy:
        return LEGACY_RULES["sapiens"]
    return INTPATHV2_RULES
