#!/usr/bin/env bash
#
# consensus.sh — Clinical germline consensus annotation (DeepVariant backbone)
# -----------------------------------------------------------------------------
# Takes per-caller germline VCFs from one sample (as produced by nf-core/sarek:
# DeepVariant, Strelka2, FreeBayes, HaplotypeCaller) and produces a single
# clinically-reportable VCF that:
#
#   * uses the DeepVariant VCF as the GENOTYPED BACKBONE  (keeps GT/DP/AD/GQ/VAF)
#   * annotates every variant with cross-caller concordance:
#         CALLERS=<comma list>   NCALLERS=<int>   CONF=HIGH|MEDIUM|LOW
#   * DROPS NOTHING — the keep/reject decision (e.g. "NCALLERS>=2 OR DeepVariant")
#     is deferred to the downstream candidate-filtering stage, fully auditable.
#
# Why this shape: Sarek's built-in --snv_consensus_calling intersects with
# `bcftools isec -n+2` and writes a SITES-ONLY consensus VCF (QUAL/FILTER='.',
# no genotype columns) — not reportable. We reuse the same bcftools building
# blocks but keep DeepVariant's genotypes and tier rather than hard-filter.
#
# CONF tiers reflect CROSS-CALLER CONCORDANCE, not absolute quality:
#   NCALLERS>=3 -> HIGH    NCALLERS==2 -> MEDIUM    NCALLERS==1 -> LOW
# A DeepVariant-only call is CONF=LOW = "unconfirmed by other callers", NOT "bad".
#
# Requires: bcftools, bgzip, tabix (htslib) on PATH.  Tested with bcftools >=1.18.
# -----------------------------------------------------------------------------
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  consensus.sh -r REF.fasta -d DEEPVARIANT.vcf.gz -o OUT_PREFIX \
               -c strelka=STRELKA.vcf.gz \
               -c freebayes=FREEBAYES.vcf.gz \
               -c haplotypecaller=HAPLOTYPECALLER.vcf.gz \
               [-f 'PASS,.']   [-w WORKDIR]

  -r  Reference FASTA (indexed .fai alongside) — used for indel left-alignment.
  -d  DeepVariant VCF (bgzipped) — the genotyped backbone. Required.
  -o  Output prefix. Writes <prefix>.consensus.vcf.gz (+ .tbi) and
      <prefix>.consensus.log.
  -c  name=path for each OTHER caller (repeatable). 'deepvariant' is implicit via -d.
  -f  FILTER values to keep (comma list passed to `bcftools view -f`). Default: PASS,.
  -w  Working dir for intermediates (default: a temp dir, auto-cleaned).
  -h  This help.
EOF
}

FILTER_KEEP="PASS,."
WORKDIR=""
REF=""
DV_VCF=""
OUT_PREFIX=""
declare -a OTHER_NAMES=()
declare -a OTHER_PATHS=()

while getopts ":r:d:o:c:f:w:h" opt; do
  case "$opt" in
    r) REF="$OPTARG" ;;
    d) DV_VCF="$OPTARG" ;;
    o) OUT_PREFIX="$OPTARG" ;;
    c) name="${OPTARG%%=*}"; path="${OPTARG#*=}"
       if [[ "$name" == "$OPTARG" || -z "$name" || -z "$path" ]]; then
         echo "ERROR: -c expects name=path, got '$OPTARG'" >&2; exit 2; fi
       OTHER_NAMES+=("$name"); OTHER_PATHS+=("$path") ;;
    f) FILTER_KEEP="$OPTARG" ;;
    w) WORKDIR="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    :) echo "ERROR: -$OPTARG needs an argument" >&2; exit 2 ;;
  esac
done

# ---- validation -------------------------------------------------------------
[[ -n "$REF" && -n "$DV_VCF" && -n "$OUT_PREFIX" ]] || { echo "ERROR: -r, -d, -o are required." >&2; usage >&2; exit 2; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH." >&2; exit 3; }
command -v bgzip   >/dev/null || { echo "ERROR: bgzip not on PATH." >&2; exit 3; }
command -v tabix   >/dev/null || { echo "ERROR: tabix not on PATH." >&2; exit 3; }
[[ -f "$REF" ]]      || { echo "ERROR: reference not found: $REF" >&2; exit 3; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: reference index missing: ${REF}.fai (run: samtools faidx $REF)" >&2; exit 3; }
[[ -f "$DV_VCF" ]]   || { echo "ERROR: DeepVariant VCF not found: $DV_VCF" >&2; exit 3; }
for p in "${OTHER_PATHS[@]:-}"; do [[ -z "$p" || -f "$p" ]] || { echo "ERROR: caller VCF not found: $p" >&2; exit 3; }; done

if [[ -z "$WORKDIR" ]]; then
  WORKDIR="$(mktemp -d)"; trap 'rm -rf "$WORKDIR"' EXIT
else
  mkdir -p "$WORKDIR"
fi

LOG="${OUT_PREFIX}.consensus.log"
mkdir -p "$(dirname "$OUT_PREFIX")"
: > "$LOG"
log() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG" >&2; }

log "consensus.sh start"
log "bcftools: $(bcftools --version | head -1)"
log "reference: $REF"
log "backbone (DeepVariant): $DV_VCF"
log "filter-keep: $FILTER_KEEP"

# ---- normalize one VCF: keep PASS, left-align + split multiallelics, dedup ---
# DeepVariant keeps its FORMAT (genotypes); others are normalized identically so
# their alleles match the backbone exactly for intersection.
norm_one() {
  local in="$1" out="$2"
  # Two passes: (1) keep wanted FILTERs, split multiallelics + left-align against
  # the reference; (2) drop ONLY exact-duplicate records. Split intentionally —
  # older bcftools (<=1.13) refuses -m and --rm-dup in a single `norm` call.
  # NB: --rm-dup exact (NOT all): `all` collapses by POSITION and would silently
  # delete alternate alleles of a just-split multiallelic site. `exact` removes
  # only records identical in CHROM/POS/REF/ALT, preserving distinct ALTs.
  { bcftools view -f "$FILTER_KEEP" "$in" \
      | bcftools norm -f "$REF" -m -both \
      | bcftools norm --rm-dup exact -Oz -o "$out" ; } 2>>"$LOG"
  tabix -f -p vcf "$out" 2>>"$LOG"
}

# Order is fixed and recorded: DeepVariant FIRST, then the others as given.
# isec's presence string has one char per input in THIS order.
declare -a CALLER_ORDER=("deepvariant")
declare -a NORM_VCFS=()

DV_NORM="${WORKDIR}/00.deepvariant.norm.vcf.gz"
log "normalizing deepvariant ..."
norm_one "$DV_VCF" "$DV_NORM"
NORM_VCFS+=("$DV_NORM")

idx=1
for i in "${!OTHER_NAMES[@]}"; do
  nm="${OTHER_NAMES[$i]}"; pth="${OTHER_PATHS[$i]}"
  outv="$(printf '%s/%02d.%s.norm.vcf.gz' "$WORKDIR" "$idx" "$nm")"
  log "normalizing ${nm} ..."
  norm_one "$pth" "$outv"
  CALLER_ORDER+=("$nm"); NORM_VCFS+=("$outv"); idx=$((idx+1))
done
log "callers (presence-string order): ${CALLER_ORDER[*]}"

# ---- intersect ALL callers, -n+1 = keep every site, record presence ---------
ISEC_DIR="${WORKDIR}/isec"
bcftools isec -n +1 -p "$ISEC_DIR" "${NORM_VCFS[@]}" 2>>"$LOG"
[[ -s "${ISEC_DIR}/sites.txt" ]] || { log "ERROR: isec produced no sites.txt"; exit 4; }
log "isec sites: $(wc -l < "${ISEC_DIR}/sites.txt")"

# ---- build a presence-annotation VCF (CALLERS/NCALLERS/CONF for every site) --
# sites.txt cols: CHROM POS REF ALT <binary-presence-string>
PRESENCE="${WORKDIR}/presence.vcf.gz"
awk -v callers="$(IFS=,; echo "${CALLER_ORDER[*]}")" '
  BEGIN{
    OFS="\t"; n=split(callers,c,",")
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"Callers that found this variant\">"
    print "##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description=\"Number of callers concurring\">"
    print "##INFO=<ID=CONF,Number=1,Type=String,Description=\"Cross-caller concordance tier: HIGH(>=3)/MEDIUM(2)/LOW(1)\">"
    print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"
  }
  {
    m=split($5,bits,""); list=""; k=0
    for(i=1;i<=n;i++) if(bits[i]=="1"){ list=list (list?",":"") c[i]; k++ }
    conf = (k>=3?"HIGH":(k==2?"MEDIUM":"LOW"))
    print $1,$2,".",$3,$4,".",".","CALLERS=" list ";NCALLERS=" k ";CONF=" conf
  }' "${ISEC_DIR}/sites.txt" | bgzip -c > "$PRESENCE"
tabix -f -p vcf "$PRESENCE" 2>>"$LOG"

# ---- PRIMARY: annotate the DeepVariant backbone (genotypes preserved) --------
OUT_MAIN="${OUT_PREFIX}.consensus.vcf.gz"
bcftools annotate -a "$PRESENCE" -c INFO/CALLERS,INFO/NCALLERS,INFO/CONF \
  -Oz -o "$OUT_MAIN" "$DV_NORM" 2>>"$LOG"
tabix -f -p vcf "$OUT_MAIN" 2>>"$LOG"
log "consensus (DeepVariant backbone): $OUT_MAIN  [$(bcftools view -H "$OUT_MAIN" | wc -l) variants]"

log "consensus.sh done"
echo "$OUT_MAIN"
