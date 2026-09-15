#!/usr/bin/env bash
# Build a self-contained WebAssembly build of the target-discovery app.
# Output is a static site that runs entirely in the browser - no server, no R
# on the host. Deployed to GitHub Pages under the lab website.
#
#   bash app/build_shinylive.sh [OUTDIR]     # run from the repository root
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/app/_site}"
BUNDLE="$(mktemp -d)/aml-atlas"
mkdir -p "$BUNDLE/data"

# app.R, with the data path pinned to the bundled copy
sed 's#for (p in c("source_data", file.path("..", "source_data"),#for (p in c("data", file.path("..", "data"),#; s#file.path(dirname(getwd()), "source_data")))#file.path(dirname(getwd()), "data")))#' \
    "$ROOT/app/app.R" > "$BUNDLE/app.R"

# copy only the CSVs the app actually reads
grep -oE '<- rd\([^)]*\)' "$ROOT/app/app.R" \
  | sed -E 's/<- rd\(|\)//g; s/"//g; s/, */\//g' \
  | sort -u | while IFS= read -r rel; do
      src="$ROOT/source_data/$rel"
      [ -f "$src" ] || { echo "MISSING: $rel" >&2; exit 1; }
      mkdir -p "$BUNDLE/data/$(dirname "$rel")"
      cp "$src" "$BUNDLE/data/$rel"
    done
echo "bundle: $(du -sh "$BUNDLE" | cut -f1) across $(find "$BUNDLE/data" -type f | wc -l | tr -d ' ') data files"

Rscript -e "shinylive::export('$BUNDLE', '$OUT')"
echo "built -> $OUT ($(du -sh "$OUT" | cut -f1))"
