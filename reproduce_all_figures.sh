#!/usr/bin/env bash
# Reproduce every figure that can be rebuilt from the shipped data.
#
#   bash Code/reproduce_all_figures.sh            # everything
#   bash Code/reproduce_all_figures.sh Figure_5   # one figure
#
# Layout required (see README_REPRODUCIBILITY.md):
#   <parent>/Code/   <parent>/Source_Data/
#
# Renderers write into each figure's figures/ folder. Nothing outside the
# unpacked archives is written and no network access is used.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
ONLY="${1:-}"
LOG="$ROOT/reproduce_logs"; rm -rf "$LOG"; mkdir -p "$LOG"

# ---- dependency check ----------------------------------------------------
echo "checking dependencies"
MISS=0
command -v Rscript >/dev/null || { echo "  MISSING: R (Rscript)"; MISS=1; }
command -v python3 >/dev/null || { echo "  MISSING: python3"; MISS=1; }
for p in pandas numpy scipy matplotlib anndata seaborn scanpy wordcloud; do
  python3 -c "import $p" 2>/dev/null || { echo "  MISSING python package: $p"; MISS=1; }
done
for p in ggplot2 dplyr scales patchwork tidyr ggrepel readr png ggdendro flexmix reshape; do
  Rscript -e "q(status=!requireNamespace('$p',quietly=TRUE))" 2>/dev/null || { echo "  MISSING R package: $p"; MISS=1; }
done
[ $MISS -eq 0 ] && echo "  all core dependencies present"
echo "  install anything listed above, then re-run"
echo

PASS=(); FAIL=(); GATED=()

# Several renderers ggsave() into an output folder they never create. Pre-make
# every output directory the shipped scripts name, so that is not a failure.
echo "preparing output directories"
grep -rhoE '0[45]_(Main|Supplementary)_Figures/[A-Za-z0-9_./-]+' . --include="*.R" --include="*.py" 2>/dev/null \
  | sed 's|/[^/]*\.[A-Za-z0-9]*$||' | sort -u | while read -r d; do
      case "$d" in */*) mkdir -p "../$d" 2>/dev/null ;; esac
    done
mkdir -p ../04_Main_Figures ../05_Supplementary_Figures 2>/dev/null
echo "  $(find ../04_Main_Figures ../05_Supplementary_Figures -type d 2>/dev/null | wc -l | tr -d ' ') output directories ready"
echo
is_builder() {
  case "$(basename "$1")" in
    build_*|compute_*|explore_*|cellchat_*|0[0-7]_*) return 0 ;; *) return 1 ;;
  esac
}
is_composite() { case "$(basename "$1")" in composite*|annotate*) return 0 ;; *) return 1 ;; esac; }
is_renderer()  {
  is_builder "$1" && return 1
  case "$(basename "$1")" in fig*|panel*|composite*|annotate*|*Dendrogram*) return 0 ;; *) return 1 ;; esac
}

run_one() {
  local s="$1" figdir="$2" sdir="$3"
  local name="${figdir}/$(basename "$s")"
  local out="$LOG/$(echo "$name" | tr '/' '_').log"
  if [[ "$s" == *.R ]]; then ( cd "$sdir" && Rscript "$(basename "$s")" ) >"$out" 2>&1
  else                       ( cd "$sdir" && python3 "$(basename "$s")" ) >"$out" 2>&1; fi
  local rc=$?
  if [ $rc -eq 0 ]; then PASS+=("$name"); return 0; fi
  if grep -qiE "\.\./raw/|raw/ALSF|No such file.*raw/" "$out"; then GATED+=("$name"); return 0; fi
  FAIL+=("$name"); return 1
}

# Two sweeps. Some panels consume PNGs that other panels write, so a second
# sweep picks up whatever became available during the first. Within each sweep
# panels run before composites, because composites stitch panel outputs.
for sweep in 1 2; do
 [ $sweep -eq 2 ] && { PASS=(); FAIL=(); GATED=(); echo "second sweep"; }
 for phase in panels composites; do
  for figdir in $(ls -d Figure_* 2>/dev/null | sort -V); do
    [ -n "$ONLY" ] && [ "$figdir" != "$ONLY" ] && continue
    sdir="$figdir/scripts"; [ -d "$sdir" ] || sdir="$figdir"
    mkdir -p "$figdir/figures"
    for s in "$sdir"/*.R "$sdir"/*.py; do
      [ -e "$s" ] || continue
      is_renderer "$s" || continue
      if [ "$phase" = panels ]; then is_composite "$s" && continue
      else is_composite "$s" || continue; fi
      run_one "$s" "$figdir" "$sdir"
    done
  done
 done
done

echo "================ reproduction summary ================"
printf '  reproduced        : %d\n' "${#PASS[@]}"
printf '  needs raw/ inputs : %d   (45 GB atlas, see README_REPRODUCIBILITY.md)\n' "${#GATED[@]}"
printf '  failed            : %d\n' "${#FAIL[@]}"
echo
[ ${#GATED[@]} -gt 0 ] && { echo "  gated on raw data:"; printf '    %s\n' "${GATED[@]}"; echo; }
[ ${#FAIL[@]}  -gt 0 ] && { echo "  failed (see reproduce_logs/):"; printf '    %s\n' "${FAIL[@]}"; echo; }
echo "  outputs are in each Figure_*/figures/ folder"
[ ${#FAIL[@]} -eq 0 ] && exit 0 || exit 1
