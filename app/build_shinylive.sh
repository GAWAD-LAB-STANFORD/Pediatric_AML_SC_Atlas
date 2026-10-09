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
mkdir -p "$BUNDLE"

cp "$ROOT/app/app.R" "$BUNDLE/app.R"

# One binary bundle instead of 24 MB of CSV: see app/prep_bundle.R for why.
Rscript "$ROOT/app/prep_bundle.R" "$ROOT/source_data" "$BUNDLE/bundle.rds"

Rscript -e "shinylive::export('$BUNDLE', '$OUT')"
echo "app.json: $(du -h "$OUT/app.json" | cut -f1)"
echo "built -> $OUT ($(du -sh "$OUT" | cut -f1))"
