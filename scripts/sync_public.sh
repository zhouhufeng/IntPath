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
mkdir -p "$DST/docs" "$DST/legacy/java" "$DST/stats" "$DST/scripts"
for f in METHODS.md DATA_FORMATS.md ROADMAP.md DEPLOY.md OLD_INTPATH.md; do cp "$SRC/docs/$f" "$DST/docs/$f"; done
$R "$SRC/Scripts/" "$DST/legacy/java/"
# never publish database credentials from the legacy code
find "$DST/legacy/java" -name '*.java' -exec sed -i -E 's/(user=)[^&"]*/\1INTPATH_DB_USER/; s/(password=)[^&"]*/\1INTPATH_DB_PASSWORD/' {} +
if grep -rqE 'password=[^I&"]' "$DST/legacy/java"; then echo "credential redaction failed" >&2; exit 1; fi
cp "$SRC/docs/PUBLIC_README.md" "$DST/README.md"
cp "$SRC/pyproject.toml" "$SRC/Dockerfile" "$SRC/.gitignore" "$DST/"
cp "$SRC/scripts/sync_public.sh" "$DST/scripts/"

# build summaries only (counts, no data rows)
for s in "$SRC"/Data/intpathv2/release/*/stats.json; do
  [ -f "$s" ] || continue
  org=$(basename "$(dirname "$s")")
  python3 - "$s" "$DST/stats/${org}_stats.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
d.pop("files", None)  # local paths
json.dump(d, open(sys.argv[2], "w"), indent=2)
PY
done
echo "synced -> $DST"
git -C "$DST" status --short | head -50
