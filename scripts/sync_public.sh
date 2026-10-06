#!/usr/bin/env bash
# Copy the shareable part of this private repo into the public IntPath repo.
#   scripts/sync_public.sh [PUBLIC_DIR]      (default: ../IntPath)
# Only whitelisted paths are copied; Data/ (raw, KEGG/BioCyc-derived) never is.
# Review with `git -C ../IntPath status` / `diff`, then commit and push there.
set -euo pipefail
SRC="$(cd "$(dirname "$0")/.." && pwd)"
DST="$(cd "${1:-$SRC/../IntPath}" && pwd)"
[ -d "$DST/.git" ] || { echo "not a git repo: $DST" >&2; exit 1; }

R="rsync -a --delete --exclude __pycache__ --exclude .pytest_cache --exclude *.egg-info"
$R "$SRC/src/" "$DST/src/"
$R "$SRC/web/" "$DST/web/"
$R "$SRC/tests/" "$DST/tests/"
mkdir -p "$DST/Docs" "$DST/legacy/java" "$DST/stats" "$DST/scripts"
# public docs are the source-free versions in Docs/Public/ (the private Docs/ keep full detail)
rm -f "$DST"/Docs/*.md
cp "$SRC"/Docs/Public/*.md "$DST/Docs/"
$R "$SRC/Scripts/" "$DST/legacy/java/"
# never publish database credentials from the legacy code
find "$DST/legacy/java" -name '*.java' -exec sed -i -E 's/(user=)[^&"]*/\1INTPATH_DB_USER/; s/(password=)[^&"]*/\1INTPATH_DB_PASSWORD/' {} +
if grep -rqE 'password=[^I&"]' "$DST/legacy/java"; then echo "credential redaction failed" >&2; exit 1; fi
cp "$SRC/Docs/PUBLIC_README.md" "$DST/README.md"
cp "$SRC/pyproject.toml" "$SRC/Dockerfile" "$SRC/.gitignore" "$DST/"
cp "$SRC/scripts/sync_public.sh" "$DST/scripts/"

# build summaries only: totals, no per-source breakdown, no data rows
rm -f "$DST"/stats/*.json
for s in "$SRC"/Data/intpathv2/db/*/stats.json; do
  [ -f "$s" ] || continue
  org=$(basename "$(dirname "$s")")
  python3 - "$s" "$DST/stats/${org}_stats.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
p, r = d["pathways"], d.get("merge_review", {})
out = {
    "organism": d["organism"], "taxid": d["taxid"], "built": d["built"],
    "source_pathways": sum(v.get("pathways") or 0 for v in d.get("sources", {}).values()),
    "related_pathway_pairs": r.get("accepted"),
    "integrated_pathways": p["sets"], "merged_pathways": p["merged_sets"],
    "source_pathways_in_merged": p["member_pathways_in_merged_sets"],
    "genes": p["genes"], "gene_pairs": p["gene_pairs"],
    "physical_ppi_edges": d.get("ppi", {}).get("edges"),
    "go_sets": d.get("go", {}).get("go_sets"), "msigdb_sets": d.get("msigdb", {}).get("sets"),
}
json.dump({k: v for k, v in out.items() if v is not None}, open(sys.argv[2], "w"), indent=2)
PY
done
echo "synced -> $DST"
git -C "$DST" status --short | head -50
