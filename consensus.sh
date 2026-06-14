#!/usr/bin/env bash
#
# consensus.sh — Clinical germline UNION consensus (DeepVariant + rescued concordant)
# -----------------------------------------------------------------------------
# Builds one clinically-reportable VCF for a sample from its per-caller germline
# VCFs (nf-core/sarek: DeepVariant, Strelka2, FreeBayes, HaplotypeCaller):
#
#   1. BACKBONE  — every DeepVariant call is kept, with DeepVariant's genotype
#                  fields (GT/GQ/DP/AD/VAF). GT_SOURCE=deepvariant.
#   2. RESCUE    — any variant DeepVariant did NOT call but >=2 of the OTHER
#                  callers did, is added back. Its genotype is borrowed from
#                  Strelka2 if Strelka2 called it, otherwise HaplotypeCaller
#                  (FreeBayes counts toward concordance but is never a genotype
#                  donor). GT_SOURCE=strelka|haplotypecaller. Rescued records
#                  carry GT/GQ/DP/AD (no VAF — those callers don't emit it; AD
#                  is kept so allele fraction stays derivable).
#
# Every variant is annotated with:
#   CALLERS=<list>   NCALLERS=<int>   CONF=HIGH(>=3)/MEDIUM(2)/LOW(1)   GT_SOURCE=<caller>
#
# Nothing else is dropped/tiered here — downstream candidate-filtering decides
# strictness (e.g. on NCALLERS/CONF). CONF reflects CROSS-CALLER CONCORDANCE,
# not absolute quality (a DeepVariant-only LOW call is unconfirmed, not wrong).
#
# Requires: bcftools, bgzip, tabix (htslib). Tested on bcftools 1.13 + newer.
# -----------------------------------------------------------------------------
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  consensus.sh -r REF.fasta -d DEEPVARIANT.vcf.gz -o OUT_PREFIX \
               -c strelka=STRELKA.vcf.gz \
               -c freebayes=FREEBAYES.vcf.gz \
               -c haplotypecaller=HAPLOTYPECALLER.vcf.gz \
               [-s SAMPLE_ID] [-f 'PASS,.'] [-w WORKDIR]

  -r  Reference FASTA (indexed .fai alongside) — used for indel left-alignment.
  -d  DeepVariant VCF (bgzipped) — the genotyped backbone. Required.
  -o  Output prefix. Writes <prefix>.consensus.vcf.gz (+ .tbi) and <prefix>.consensus.log.
  -c  name=path for each OTHER caller (repeatable). Recognized genotype donors,
      in priority order: 'strelka' then 'haplotypecaller'. 'freebayes' counts
      toward concordance only.
  -s  Sample name to use in the output (default: basename of OUT_PREFIX).
  -f  FILTER values to keep (passed to `bcftools view -f`). Default: PASS,.
  -w  Working dir for intermediates (default: a temp dir, auto-cleaned).
  -h  This help.
EOF
}

FILTER_KEEP="PASS,."
WORKDIR=""
REF=""
DV_VCF=""
OUT_PREFIX=""
SAMPLE_OVERRIDE=""
declare -a OTHER_NAMES=()
declare -a OTHER_PATHS=()
# Genotype-donor priority for rescued (non-DeepVariant) variants:
declare -a DONOR_PRIORITY=("strelka" "haplotypecaller")

while getopts ":r:d:o:c:s:f:w:h" opt; do
  case "$opt" in
    r) REF="$OPTARG" ;;
    d) DV_VCF="$OPTARG" ;;
    o) OUT_PREFIX="$OPTARG" ;;
    c) name="${OPTARG%%=*}"; path="${OPTARG#*=}"
       if [[ "$name" == "$OPTARG" || -z "$name" || -z "$path" ]]; then
         echo "ERROR: -c expects name=path, got '$OPTARG'" >&2; exit 2; fi
       OTHER_NAMES+=("$name"); OTHER_PATHS+=("$path") ;;
    s) SAMPLE_OVERRIDE="$OPTARG" ;;
    f) FILTER_KEEP="$OPTARG" ;;
    w) WORKDIR="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    :) echo "ERROR: -$OPTARG needs an argument" >&2; exit 2 ;;
  esac
done

# ---- validation -------------------------------------------------------------
[[ -n "$REF" && -n "$DV_VCF" && -n "$OUT_PREFIX" ]] || { echo "ERROR: -r, -d, -o are required." >&2; usage >&2; exit 2; }
for tool in bcftools bgzip tabix; do command -v "$tool" >/dev/null || { echo "ERROR: $tool not on PATH." >&2; exit 3; }; done
[[ -f "$REF" ]]       || { echo "ERROR: reference not found: $REF" >&2; exit 3; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: reference index missing: ${REF}.fai (run: samtools faidx $REF)" >&2; exit 3; }
[[ -f "$DV_VCF" ]]    || { echo "ERROR: DeepVariant VCF not found: $DV_VCF" >&2; exit 3; }
for p in "${OTHER_PATHS[@]:-}"; do [[ -z "$p" || -f "$p" ]] || { echo "ERROR: caller VCF not found: $p" >&2; exit 3; }; done

SAMPLE_ID="${SAMPLE_OVERRIDE:-$(basename "$OUT_PREFIX")}"

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
log "reference: $REF   sample: $SAMPLE_ID   filter-keep: $FILTER_KEEP"
log "backbone (DeepVariant): $DV_VCF"

# Fields kept on every output record:
KEEP_FMT_TAGS=" GT GQ DP AD VAF "                       # harmonized genotype schema
KEEP_INFO_TAGS=" CALLERS NCALLERS CONF GT_SOURCE "      # our consensus annotations
# Harmonize a VCF to the kept fields by DROPPING every other INFO/FORMAT tag.
# We compute the drop-list from the header rather than use `-x ^keep,list`:
# bcftools 1.13's caret form errors ("No matching tag") when nothing needs
# removing, and also mishandles mixed INFO+FORMAT caret lists.
harmonize() {
  local in="$1" out="$2" drop="" tag
  for tag in $(bcftools view -h "$in" | sed -n 's/^##FORMAT=<ID=\([^,]*\),.*/\1/p' | sort -u); do
    [[ "$KEEP_FMT_TAGS" == *" $tag "* ]] || drop="${drop:+$drop,}FORMAT/$tag"
  done
  for tag in $(bcftools view -h "$in" | sed -n 's/^##INFO=<ID=\([^,]*\),.*/\1/p' | sort -u); do
    [[ "$KEEP_INFO_TAGS" == *" $tag "* ]] || drop="${drop:+$drop,}INFO/$tag"
  done
  if [[ -n "$drop" ]]; then bcftools annotate -x "$drop" "$in" -Oz -o "$out" 2>>"$LOG"
  else cp "$in" "$out"; fi
  tabix -f -p vcf "$out" 2>>"$LOG"
}

# ---- normalize one VCF: keep PASS, split multiallelics + left-align, dedup,
#      and harmonize the single sample name so files concat cleanly later -----
norm_one() {
  local in="$1" out="$2" tmp="$2.tmp.vcf.gz"
  # --rm-dup exact (NOT all): `all` collapses by POSITION, silently deleting
  # alternate alleles of a just-split multiallelic. Two passes because older
  # bcftools (<=1.13) refuses -m and --rm-dup in one `norm` call.
  { bcftools view -f "$FILTER_KEEP" "$in" \
      | bcftools norm -f "$REF" -m -both \
      | bcftools norm --rm-dup exact -Oz -o "$tmp" ; } 2>>"$LOG"
  local ns; ns=$(bcftools query -l "$tmp" | wc -l)
  if [[ "$ns" -ge 1 ]]; then
    printf '%s\n' "$SAMPLE_ID" | bcftools reheader -s - "$tmp" -o "$out" 2>>"$LOG"; rm -f "$tmp"
  else
    mv "$tmp" "$out"
  fi
  tabix -f -p vcf "$out" 2>>"$LOG"
}

# Fixed order (DeepVariant first); presence-string bits follow this order.
declare -a CALLER_ORDER=("deepvariant")
declare -a NORM_VCFS=()
declare -A NORM_BY_NAME=()

DV_NORM="${WORKDIR}/00.deepvariant.norm.vcf.gz"
log "normalizing deepvariant ..."
norm_one "$DV_VCF" "$DV_NORM"; NORM_VCFS+=("$DV_NORM"); NORM_BY_NAME[deepvariant]="$DV_NORM"

idx=1
for i in "${!OTHER_NAMES[@]}"; do
  nm="${OTHER_NAMES[$i]}"; outv="$(printf '%s/%02d.%s.norm.vcf.gz' "$WORKDIR" "$idx" "$nm")"
  log "normalizing ${nm} ..."
  norm_one "${OTHER_PATHS[$i]}" "$outv"
  CALLER_ORDER+=("$nm"); NORM_VCFS+=("$outv"); NORM_BY_NAME["$nm"]="$outv"; idx=$((idx+1))
done
log "callers (presence-string order): ${CALLER_ORDER[*]}"

# ---- intersect ALL callers, -n+1 = every site, record presence --------------
ISEC_DIR="${WORKDIR}/isec"; bcftools isec -n +1 -p "$ISEC_DIR" "${NORM_VCFS[@]}" 2>>"$LOG"
[[ -s "${ISEC_DIR}/sites.txt" ]] || { log "ERROR: isec produced no sites.txt"; exit 4; }
log "isec sites: $(wc -l < "${ISEC_DIR}/sites.txt")"

# ---- presence VCF: CALLERS/NCALLERS/CONF + the chosen GT_SOURCE per site -----
PRESENCE="${WORKDIR}/presence.vcf.gz"
awk -v callers="$(IFS=,; echo "${CALLER_ORDER[*]}")" -v donors="$(IFS=,; echo "${DONOR_PRIORITY[*]}")" '
  BEGIN{
    OFS="\t"; n=split(callers,c,","); nd=split(donors,d,",")
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"Callers that found this variant\">"
    print "##INFO=<ID=NCALLERS,Number=1,Type=Integer,Description=\"Number of callers concurring\">"
    print "##INFO=<ID=CONF,Number=1,Type=String,Description=\"Concordance tier: HIGH(>=3)/MEDIUM(2)/LOW(1)\">"
    print "##INFO=<ID=GT_SOURCE,Number=1,Type=String,Description=\"Caller the output genotype was taken from\">"
    print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"
  }
  {
    m=split($5,bits,""); list=""; k=0; delete pres
    for(i=1;i<=n;i++) if(bits[i]=="1"){ list=list (list?",":"") c[i]; k++; pres[c[i]]=1 }
    conf=(k>=3?"HIGH":(k==2?"MEDIUM":"LOW"))
    if("deepvariant" in pres) src="deepvariant"
    else { src="."; for(j=1;j<=nd;j++) if(d[j] in pres){ src=d[j]; break } }
    print $1,$2,".",$3,$4,".",".","CALLERS=" list ";NCALLERS=" k ";CONF=" conf ";GT_SOURCE=" src
  }' "${ISEC_DIR}/sites.txt" | bgzip -c > "$PRESENCE"
tabix -f -p vcf "$PRESENCE" 2>>"$LOG"

# ---- BACKBONE: annotate DeepVariant records, harmonize to kept fields -------
PRIMARY="${WORKDIR}/primary.vcf.gz"; PRIMARY_ANN="${WORKDIR}/primary.ann.vcf.gz"
bcftools annotate -a "$PRESENCE" -c INFO/CALLERS,INFO/NCALLERS,INFO/CONF,INFO/GT_SOURCE "$DV_NORM" -Oz -o "$PRIMARY_ANN" 2>>"$LOG"
harmonize "$PRIMARY_ANN" "$PRIMARY"
log "backbone (DeepVariant): $(bcftools view -H "$PRIMARY" | wc -l) variants"

declare -a CONCAT=("$PRIMARY")

# ---- RESCUE: non-DeepVariant sites with NCALLERS>=2, genotype from donor ----
for donor in "${DONOR_PRIORITY[@]}"; do
  [[ -n "${NORM_BY_NAME[$donor]:-}" ]] || continue
  sites="${WORKDIR}/rescue_sites.${donor}.vcf.gz"
  bcftools view -i "NCALLERS>=2 && GT_SOURCE=\"$donor\"" "$PRESENCE" -Oz -o "$sites" 2>>"$LOG"
  tabix -f -p vcf "$sites" 2>>"$LOG"
  nsite=$(bcftools view -H "$sites" | wc -l)
  [[ "$nsite" -gt 0 ]] || { log "rescue via ${donor}: 0 sites"; continue; }
  # donor's genotyped records at exactly those sites
  rdir="${WORKDIR}/isec_rescue_${donor}"
  bcftools isec -n=2 -w1 -p "$rdir" -Oz "${NORM_BY_NAME[$donor]}" "$sites" 2>>"$LOG"
  rescann="${WORKDIR}/rescue.${donor}.ann.vcf.gz"; resc="${WORKDIR}/rescue.${donor}.vcf.gz"
  bcftools annotate -a "$PRESENCE" -c INFO/CALLERS,INFO/NCALLERS,INFO/CONF,INFO/GT_SOURCE "${rdir}/0000.vcf.gz" -Oz -o "$rescann" 2>>"$LOG"
  harmonize "$rescann" "$resc"
  log "rescue via ${donor}: $(bcftools view -H "$resc" | wc -l) variants"
  CONCAT+=("$resc")
done

# ---- assemble: concat backbone + rescued, sort, index -----------------------
OUT_MAIN="${OUT_PREFIX}.consensus.vcf.gz"
{ bcftools concat -a "${CONCAT[@]}" | bcftools sort -Oz -o "$OUT_MAIN" ; } 2>>"$LOG"
tabix -f -p vcf "$OUT_MAIN" 2>>"$LOG"
log "FINAL consensus: $OUT_MAIN  [$(bcftools view -H "$OUT_MAIN" | wc -l) variants total]"
log "consensus.sh done"
echo "$OUT_MAIN"
