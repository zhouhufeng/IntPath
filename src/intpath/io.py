"""Readers and writers: legacy IntPath text files, GMT, and the IntPathV2 release tables."""

from __future__ import annotations

import gzip
import io
import json
import zipfile
from pathlib import Path
from typing import Iterable, Iterator

from .model import GeneSet, SourcePathway


def open_text(path: str | Path) -> io.TextIOBase:
    """Open plain, .gz or single-member .zip files as text."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if path.suffix == ".zip":
        zf = zipfile.ZipFile(path)
        name = next(n for n in zf.namelist() if not n.endswith("/"))
        return io.TextIOWrapper(zf.open(name), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------- #
# Legacy (old IntPath, 2012-2021) normalized files
#   *NormPthGEN : pathway \t gene \t source
#   *NormPthGPR : geneA \t geneB \t relation \t pathway \t source
# --------------------------------------------------------------------------- #
def read_legacy_source(gen_file: str | Path, gpr_file: str | Path | None, source: str) -> list[SourcePathway]:
    """Read one source's normalized pathway-gene (+ gene pair) files, keeping file order."""
    pws: dict[str, SourcePathway] = {}
    with open_text(gen_file) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2 or not f[0]:
                continue
            pws.setdefault(f[0], SourcePathway(source, f[0])).genes.add(f[1])
    if gpr_file:
        with open_text(gpr_file) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                p = pws.setdefault(f[3], SourcePathway(source, f[3]))
                p.pairs.setdefault((f[0], f[1]), set()).update(f[2].split())
    return list(pws.values())


# --------------------------------------------------------------------------- #
# GMT
# --------------------------------------------------------------------------- #
def read_gmt(path: str | Path) -> Iterator[tuple[str, str, list[str]]]:
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                yield f[0], f[1], [g for g in f[2:] if g]


def write_gmt(sets: Iterable[GeneSet], path: str | Path) -> None:
    with open(path, "w") as fh:
        for gs in sets:
            fh.write("\t".join([gs.id, gs.name, *sorted(gs.genes)]) + "\n")


# --------------------------------------------------------------------------- #
# IntPathV2 release tables (tab-delimited, documented in Docs/DATA_FORMATS.md)
# --------------------------------------------------------------------------- #
def write_release(sets: list[GeneSet], outdir: str | Path, prefix: str = "intpath") -> dict[str, str]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = {
        "sets": out / f"{prefix}_genesets.tsv",
        "genes": out / f"{prefix}_set_genes.tsv",
        "pairs": out / f"{prefix}_set_genepairs.tsv",
        "members": out / f"{prefix}_set_members.tsv",
        "gmt": out / f"{prefix}.gmt",
    }
    with open(files["sets"], "w") as s, open(files["genes"], "w") as g, open(files["pairs"], "w") as p, open(
        files["members"], "w"
    ) as m:
        s.write("set_id\tname\tcollection\tsources\tn_genes\tn_pairs\tlinks\n")
        g.write("set_id\tgene\tsources\n")
        p.write("set_id\tgene_a\tgene_b\trelations\tsources\n")
        m.write("set_id\tsource\tsource_pathway\n")
        for gs in sets:
            s.write(
                f"{gs.id}\t{gs.name}\t{gs.collection}\t{','.join(gs.sources)}\t{gs.size}\t{len(gs.pairs)}\t{','.join(gs.links)}\n"
            )
            for gene in sorted(gs.genes):
                g.write(f"{gs.id}\t{gene}\t{','.join(sorted(gs.genes[gene]))}\n")
            for (a, b), e in sorted(gs.pairs.items()):
                p.write(f"{gs.id}\t{a}\t{b}\t{','.join(sorted(e['rel']))}\t{','.join(sorted(e['src']))}\n")
            for src, name in gs.members:
                m.write(f"{gs.id}\t{src}\t{name}\n")
    write_gmt(sets, files["gmt"])
    return {k: str(v) for k, v in files.items()}


def read_release(outdir: str | Path, prefix: str = "intpath") -> list[GeneSet]:
    d = Path(outdir)
    sets: dict[str, GeneSet] = {}
    with open(d / f"{prefix}_genesets.tsv") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sets[f[0]] = GeneSet(f[0], f[1], f[2], links=[x for x in f[6].split(",") if x] if len(f) > 6 else [])
    with open(d / f"{prefix}_set_genes.tsv") as fh:
        next(fh)
        for line in fh:
            sid, gene, src = line.rstrip("\n").split("\t")
            sets[sid].genes[gene] = set(src.split(","))
    pf = d / f"{prefix}_set_genepairs.tsv"
    if pf.exists():
        with open(pf) as fh:
            next(fh)
            for line in fh:
                sid, a, b, rel, src = line.rstrip("\n").split("\t")
                sets[sid].pairs[(a, b)] = {"rel": set(rel.split(",")), "src": set(src.split(","))}
    with open(d / f"{prefix}_set_members.tsv") as fh:
        next(fh)
        for line in fh:
            sid, src, name = line.rstrip("\n").split("\t")
            sets[sid].members.append((src, name))
    return list(sets.values())


def write_json(obj, path: str | Path) -> None:
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=lambda o: sorted(o) if isinstance(o, set) else str(o))
