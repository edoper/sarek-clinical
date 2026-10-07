#!/usr/bin/env bash
#
# run_bge_annotate_filter.sh: feed the local BGE consensus VCFs into VEP + candidate-filtering.
#   1. VEP-annotate each <repo>/consensus-cohort/<sample>.consensus.vcf.gz
#      -> <WD>/<sample>.germline.vep.vcf.gz   (resumable; one isolated workdir so the
#         BGE samples don't collide with other VEP VCFs in candidate-filtering).
#   2. run_filtering.sh in that workdir -> <proband>.<panel>.candidatos (auto-discovers trios/duos).
#   3. Copy the .candidatos to $WIN.
# Progress: tail this script's log, or `watch <repo>/bge_filter_progress.sh`.
# FAIL-CLOSED: any VEP failure stops the run before filtering (a partial cohort
# changes the cohort-artifact denominator), a filtering failure stops before the
# copy-out, and only tables written by THIS run are copied.
set -uo pipefail

# Same completeness test consensus_from_results.sh uses: a killed VEP leaves a file
# whose header still reads, so "header OK" is not "done". Whole = BGZF EOF block + .tbi.
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
complete_vcf() {  # <path.vcf.gz>
    [[ -s "$1" && -s "$1.tbi" ]] || return 1
    [[ "$(tail -c 28 "$1" | od -An -tx1 | tr -d ' \n')" == "$BGZF_EOF" ]]
}

. "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/site.sh"
CONS_DIR="${CONS_DIR:-$SAREK_REPO/consensus-cohort}"
WD="${WD:-$CF/bge-cohort}"   # CF + WIN come from site.sh (WIN empty => no copy-out)
VEP="$CF/vep_annotate.sh"
[ -d "$CF" ] || { echo "ERROR: candidate-filtering repo not found at '$CF'."; echo "       Set CF=/path/to/candidate-filtering (or put it beside this repo)."; exit 1; }
[ -f "$VEP" ] || { echo "ERROR: $VEP not found; is '$CF' really the candidate-filtering repo?"; exit 1; }
mkdir -p "$WD"

# Isolated workdir needs the code + reference config (filtering_r.pl reads them from cwd).
# Link whatever reference files candidate-filtering currently ships, rather than a fixed
# list: the panel is versioned (g4e-2025 -> g4e-2026 -> ...) and new tables get added
# (gnomad-mis-constraint.txt drives ACMG PP2). A hardcoded list silently degrades: a
# dangling panel link, or a missing table that just switches a criterion off.
for f in filtering_r.pl parse_pangolin.pl site.sh; do
    [ -e "$CF/$f" ] && ln -sf "$CF/$f" "$WD/$f"
done
for f in "$CF"/*.txt; do
    [ -e "$f" ] && ln -sf "$f" "$WD/$(basename "$f")"
done
[ -e "$WD/filtering_r.pl" ] || { echo "ERROR: filtering_r.pl not found in $CF"; exit 1; }
ls "$WD"/g4e-*.txt >/dev/null 2>&1 || echo "  WARN: no g4e-*.txt panel found in $CF"
# Drop any dangling links left by an earlier run against a different panel version.
find "$WD" -maxdepth 1 -xtype l -delete 2>/dev/null || true

mapfile -t VCFS < <(ls "$CONS_DIR"/*.consensus.vcf.gz 2>/dev/null)
total=${#VCFS[@]}
[ "$total" -gt 0 ] || { echo "ERROR: no consensus VCFs in $CONS_DIR"; exit 1; }
echo "[annotate] $total consensus VCFs -> VEP (workdir $WD)"

# ── Step 1: VEP (resumable, per-sample tolerant) ──
done=0; failed=()
for v in "${VCFS[@]}"; do
    s=$(basename "$v" .consensus.vcf.gz)
    out="$WD/$s.germline.vep.vcf.gz"
    if complete_vcf "$out"; then
        :                                                  # already annotated, whole; skip
    else
        rm -f "$out" "$out.tbi"
        if ! bash "$VEP" "$v" "$out" > "$WD/vep.$s.log" 2>&1; then
            echo "  ERROR: VEP failed for $s (see $WD/vep.$s.log)"; failed+=("$s")
            rm -f "$out" "$out.tbi"
        fi
    fi
    done=$((done+1))
    printf "[annotate] %d/%d done | last: %-14s | failed: %d\n" "$done" "$total" "$s" "${#failed[@]}"
done
echo "[annotate] complete: $((total-${#failed[@]}))/$total ok${failed:+; failed: ${failed[*]}}"
if [ "${#failed[@]}" -gt 0 ]; then
    echo "ERROR: ${#failed[@]} sample(s) failed VEP; NOT filtering a partial cohort. Fix and re-run (resumable)."
    exit 1
fi

# ── Step 2: candidate-filtering (Pangolin + filter, all probands) ──
echo "[filter] running candidate-filtering in $WD ..."
STAMP="$WD/.filter_started"; : > "$STAMP"
if ! WORKDIR="$WD" bash "$CF/run_filtering.sh"; then
    echo "ERROR: run_filtering.sh failed; nothing copied out (tables in $WD may be stale)."
    exit 1
fi

# ── Step 3: collect candidatos to the deliverable folder ──
# WIN is a WSL convenience (a Windows-side folder). Unset: the normal case off WSL,
# means there is nowhere to copy to, so the results simply stay in $WD.
OUT_NAME="${OUT_NAME:-bge-candidatos}"          # override for other cohorts (e.g. epigen-candidatos)
# Only tables written by this run: an older table in a reused $WD (another panel,
# a removed sample) must not ride along into the deliverable.
mapfile -t NEW < <(find "$WD" -maxdepth 1 -name '*.candidatos' -newer "$STAMP" | sort)
n=${#NEW[@]}
[ "$n" -gt 0 ] || { echo "ERROR: filtering reported success but wrote no .candidatos"; exit 1; }
if [ -n "$WIN" ]; then
    mkdir -p "$WIN/$OUT_NAME"
    cp -- "${NEW[@]}" "$WIN/$OUT_NAME/"
    echo "[done] $n candidatos -> $WIN/$OUT_NAME/"
else
    echo "[done] $n candidatos in $WD/  (set WIN=/path/to/folder to copy them out)"
fi
