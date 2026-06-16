#!/usr/bin/env bash
#
# run_bge_annotate_filter.sh — feed the local BGE consensus VCFs into VEP + candidate-filtering.
#   1. VEP-annotate each ~/sarek-clinical/consensus-cohort/<sample>.consensus.vcf.gz
#      -> <WD>/<sample>.germline.vep.vcf.gz   (resumable; one isolated workdir so the
#         63 BGE samples don't collide with the DRAGEN VEP VCFs in ~/candidate-filtering).
#   2. run_filtering.sh in that workdir -> <proband>.<panel>.candidatos (auto-discovers trios/duos).
#   3. Copy the .candidatos to $WIN.
# Progress: tail this script's log, or `watch ~/sarek-clinical/bge_filter_progress.sh`.
set -uo pipefail

CONS_DIR="${CONS_DIR:-$HOME/sarek-clinical/consensus-cohort}"
CF="${CF:-$HOME/candidate-filtering}"
WD="${WD:-$CF/bge-cohort}"
WIN="${WIN:-/mnt/c/Users/epere/Documents}"
VEP="$CF/vep_annotate.sh"
mkdir -p "$WD"

# Isolated workdir needs the code + reference config (filtering_r.pl reads them from cwd).
for f in filtering_r.pl parse_pangolin.pl mane-plus-clinical-names.txt g4e-2025.txt typevar.txt acmg_sf_v3.2.txt; do
    ln -sf "$CF/$f" "$WD/$f"
done

mapfile -t VCFS < <(ls "$CONS_DIR"/*.consensus.vcf.gz 2>/dev/null)
total=${#VCFS[@]}
[ "$total" -gt 0 ] || { echo "ERROR: no consensus VCFs in $CONS_DIR"; exit 1; }
echo "[annotate] $total consensus VCFs -> VEP (workdir $WD)"

# ── Step 1: VEP (resumable, per-sample tolerant) ──
done=0; failed=()
for v in "${VCFS[@]}"; do
    s=$(basename "$v" .consensus.vcf.gz)
    out="$WD/$s.germline.vep.vcf.gz"
    if [ -s "$out" ] && bcftools view -h "$out" >/dev/null 2>&1; then
        :                                                  # already annotated — skip
    else
        rm -f "$out" "$out.tbi"
        if ! bash "$VEP" "$v" "$out" > "$WD/vep.$s.log" 2>&1; then
            echo "  WARN: VEP failed for $s (see $WD/vep.$s.log)"; failed+=("$s")
        fi
    fi
    done=$((done+1))
    printf "[annotate] %d/%d done | last: %-14s | failed: %d\n" "$done" "$total" "$s" "${#failed[@]}"
done
echo "[annotate] complete: $((total-${#failed[@]}))/$total ok${failed:+; failed: ${failed[*]}}"

# ── Step 2: candidate-filtering (Pangolin + filter, all probands) ──
echo "[filter] running candidate-filtering in $WD ..."
WORKDIR="$WD" bash "$CF/run_filtering.sh"

# ── Step 3: collect candidatos to Windows ──
mkdir -p "$WIN/bge-candidatos"
cp "$WD"/*.candidatos "$WIN/bge-candidatos/" 2>/dev/null
echo "[done] $(ls "$WD"/*.candidatos 2>/dev/null | wc -l) candidatos -> $WIN/bge-candidatos/"
