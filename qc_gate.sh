#!/usr/bin/env bash
#
# qc_gate.sh — per-SAMPLE quality gate. Decides whether a sample is fit to report.
#
# WHY THIS EXISTS
#   The variant-level flags in candidate-filtering (lowDP, lowGQ, homopolymer) judge individual
#   calls. Nothing judged the SAMPLE. A contaminated library, a sample swap, or an exome that
#   simply failed capture would flow through the whole pipeline and produce a confident-looking
#   candidate list. That is the failure mode most likely to cause harm in routine operation —
#   far more than a small difference in caller accuracy.
#
#   This script fails such a sample loudly, BEFORE anyone reads its candidates.
#
# WHAT IT CHECKS (all from the sample's own VCF + optional CRAM/BAM — no extra cloud cost)
#   1. CALL COUNT      — a collapsed library yields far too few variants; a contaminated or
#                        badly-mapped one yields far too many.
#   2. TRANSITION/TRANSVERSION — Ti/Tv is the classic sanity metric. Genome-wide ~2.0-2.1,
#                        exome ~2.8-3.3 (exons are enriched for transitions). A value near 0.5
#                        means random noise, i.e. the calls are mostly artefact.
#   3. HETEROZYGOSITY  — het/hom-alt ratio. Contamination inflates it (foreign alleles appear
#                        heterozygous); a wrong reference or inbreeding deflates it.
#   4. CONTAMINATION   — proxy: the fraction of het calls with extreme allele balance. A clean
#                        sample's hets sit near AB 0.5; contamination drags a tail toward 0.
#   5. SEX CONCORDANCE — chrX het rate + chrY call count give observed sex; compared with the
#                        expected sex when you supply one. Catches sample swaps.
#   6. COVERAGE        — mean depth over called sites (and over a target BED if given).
#
# EXIT CODES:  0 = PASS   1 = FAIL (do not report)   2 = WARN (review before reporting)
#
# USAGE
#   ./qc_gate.sh <sample.vcf.gz> [--sex M|F] [--targets targets.bed] [--json out.json]
#
#   Thresholds are the DEFAULTS below and are deliberately conservative. A clinical lab must set
#   its own from its own validation data and record them in the SOP — see validation/SOP.md.
set -uo pipefail

VCF="${1:?Usage: $0 <sample.vcf.gz> [--sex M|F] [--targets t.bed] [--json out.json]}"; shift || true
EXPECTED_SEX=""; TARGETS=""; JSON=""
while [ $# -gt 0 ]; do
  case "$1" in
    --sex)     EXPECTED_SEX="${2:-}"; shift 2 ;;
    --targets) TARGETS="${2:-}";      shift 2 ;;
    --json)    JSON="${2:-}";         shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
command -v bcftools >/dev/null || { echo "ERROR: bcftools not on PATH" >&2; exit 1; }
[ -s "$VCF" ] || { echo "ERROR: VCF not found or empty: $VCF" >&2; exit 1; }

# ── thresholds (override by exporting; record whatever you use in the SOP) ────────────
: "${QC_MIN_VARIANTS:=20000}"      # exome floor; a WGS sample should be far above this
: "${QC_MAX_VARIANTS:=8000000}"    # above this = almost certainly artefact-laden
: "${QC_MIN_TITV:=1.8}"            # below ~1.8 => noise contaminating the call set
: "${QC_MAX_TITV:=3.6}"
: "${QC_MIN_HETHOM:=1.0}"
: "${QC_MAX_HETHOM:=3.0}"          # contamination inflates het/hom
: "${QC_MAX_SKEWED_HET_FRAC:=0.15}" # frac of hets with AB<0.25 or >0.75
: "${QC_MIN_MEAN_DP:=20}"          # mean depth at called sites

fails=(); warns=(); notes=()
add_fail(){ fails+=("$1"); }
add_warn(){ warns+=("$1"); }

SAMPLE=$(bcftools query -l "$VCF" | head -1)
[ -n "$SAMPLE" ] || { echo "ERROR: no sample column in $VCF" >&2; exit 1; }

REGION_ARG=(); [ -n "$TARGETS" ] && [ -s "$TARGETS" ] && REGION_ARG=(-T "$TARGETS")

# ── 1. call count ─────────────────────────────────────────────────────────────────────
NVAR=$(bcftools view -H "${REGION_ARG[@]}" "$VCF" 2>/dev/null | wc -l)
(( NVAR < QC_MIN_VARIANTS )) && add_fail "variant count $NVAR < $QC_MIN_VARIANTS (library/capture failure?)"
(( NVAR > QC_MAX_VARIANTS )) && add_fail "variant count $NVAR > $QC_MAX_VARIANTS (artefact-laden call set?)"

# ── 2. Ti/Tv ──────────────────────────────────────────────────────────────────────────
read -r TI TV <<<"$(bcftools query "${REGION_ARG[@]}" -f '%REF\t%ALT\n' "$VCF" 2>/dev/null | awk '
  BEGIN{ti=0;tv=0}
  length($1)==1 && length($2)==1 {
    p=$1$2
    if(p=="AG"||p=="GA"||p=="CT"||p=="TC") ti++; else if($1!=$2) tv++
  } END{print ti, tv}')"
TITV=$(awk -v a="$TI" -v b="$TV" 'BEGIN{printf "%.3f", (b>0)? a/b : 0}')
awk -v v="$TITV" -v lo="$QC_MIN_TITV" 'BEGIN{exit !(v+0 < lo+0)}' && add_fail "Ti/Tv $TITV < $QC_MIN_TITV (call set dominated by noise)"
awk -v v="$TITV" -v hi="$QC_MAX_TITV" 'BEGIN{exit !(v+0 > hi+0)}' && add_warn "Ti/Tv $TITV > $QC_MAX_TITV (unusual; check target region)"

# ── 3-4. zygosity, allele balance, depth ──────────────────────────────────────────────
read -r NHET NHOM SKEWED MEANDP <<<"$(bcftools query "${REGION_ARG[@]}" -f '[%GT\t%AD\t%DP\n]' "$VCF" 2>/dev/null | awk -F'\t' '
  BEGIN{het=0;hom=0;skew=0;dpsum=0;dpn=0}
  {
    gt=$1; ad=$2; dp=$3
    gsub(/\|/,"/",gt)
    if(gt=="0/1"||gt=="1/0") {
      het++
      if(ad ~ /,/){ split(ad,a,","); t=a[1]+a[2]; if(t>0){ ab=a[2]/t; if(ab<0.25||ab>0.75) skew++ } }
    } else if(gt=="1/1") hom++
    if(dp ~ /^[0-9]+$/){ dpsum+=dp; dpn++ }
  }
  END{ printf "%d %d %d %.1f", het, hom, skew, (dpn>0)? dpsum/dpn : 0 }')"
HETHOM=$(awk -v a="$NHET" -v b="$NHOM" 'BEGIN{printf "%.3f", (b>0)? a/b : 0}')
SKEWFRAC=$(awk -v s="$SKEWED" -v h="$NHET" 'BEGIN{printf "%.4f", (h>0)? s/h : 0}')

awk -v v="$HETHOM" -v lo="$QC_MIN_HETHOM" 'BEGIN{exit !(v+0 < lo+0)}' && add_warn "het/hom $HETHOM < $QC_MIN_HETHOM (consanguinity, or wrong reference)"
awk -v v="$HETHOM" -v hi="$QC_MAX_HETHOM" 'BEGIN{exit !(v+0 > hi+0)}' && add_fail "het/hom $HETHOM > $QC_MAX_HETHOM (possible CONTAMINATION or sample mixture)"
awk -v v="$SKEWFRAC" -v hi="$QC_MAX_SKEWED_HET_FRAC" 'BEGIN{exit !(v+0 > hi+0)}' && add_fail "skewed-AB het fraction $SKEWFRAC > $QC_MAX_SKEWED_HET_FRAC (possible CONTAMINATION)"
awk -v v="$MEANDP" -v lo="$QC_MIN_MEAN_DP" 'BEGIN{exit !(v+0 < lo+0)}' && add_fail "mean depth ${MEANDP}x < ${QC_MIN_MEAN_DP}x at called sites"

# ── 5. sex concordance ────────────────────────────────────────────────────────────────
# TWO assay-independence fixes, both learned from failing a known-male EXOME (2026-07-28):
#
#   (a) EXCLUDE THE PSEUDOAUTOSOMAL REGIONS. PAR1/PAR2 are present on both X and Y, so a
#       male is genuinely DIPLOID there and heterozygous calls are expected. Including them
#       inflated a male exome's chrX het rate from 0.090 to 0.169 and pushed it into the
#       "female" band.
#   (b) chrY CALL COUNT IS NOT COMPARABLE ACROSS ASSAYS. The same HG002 gave 11,375 chrY
#       calls by WGS and 57 by exome — a ~200x difference, because capture kits barely
#       target chrY. An absolute chrY threshold silently means "WGS only". chrY is therefore
#       NORMALISED per 1,000 autosomal calls and used only to CORROBORATE.
#
# Primary signal = chrX non-PAR heterozygosity, which is assay-independent: a male has one
# X and is near-hemizygous there; a female is heterozygous at roughly half her variant sites.
PAR1_END=2781479; PAR2_START=155701383            # GRCh38
XHET=$(bcftools query -r "chrX:$((PAR1_END+1))-$((PAR2_START-1)),X:$((PAR1_END+1))-$((PAR2_START-1))" \
        -f '[%GT\n]' "$VCF" 2>/dev/null | sed 's/|/\//' | awk '
   /^0\/1$|^1\/0$/{h++} /^1\/1$|^0\/1$|^1\/0$/{t++} END{printf "%.3f", (t>0)? h/t : -1}')
YN=$(bcftools view -H -r chrY,Y "$VCF" 2>/dev/null | wc -l)
AUTO=$(bcftools view -H -r chr1,chr2,chr3,1,2,3 "$VCF" 2>/dev/null | wc -l)
YNORM=$(awk -v y="$YN" -v a="$AUTO" 'BEGIN{printf "%.2f", (a>0)? y*1000/a : 0}')

: "${QC_XHET_MALE_MAX:=0.15}"    # non-PAR chrX het rate at/below this => male
: "${QC_XHET_FEMALE_MIN:=0.30}"  # at/above this => female; between the two => undetermined
OBSERVED_SEX="undetermined"
if   awk -v x="$XHET" -v m="$QC_XHET_MALE_MAX"   'BEGIN{exit !(x>=0 && x<=m)}'; then OBSERVED_SEX="M"
elif awk -v x="$XHET" -v f="$QC_XHET_FEMALE_MIN" 'BEGIN{exit !(x>=f)}';         then OBSERVED_SEX="F"
fi
notes+=("chrX non-PAR het-rate $XHET -> sex $OBSERVED_SEX (chrY $YN calls = $YNORM per 1000 autosomal, corroborating only)")
if [ -n "$EXPECTED_SEX" ]; then
    if [ "$OBSERVED_SEX" = "undetermined" ]; then
        add_warn "sex indeterminate (chrX non-PAR het $XHET falls between $QC_XHET_MALE_MAX and $QC_XHET_FEMALE_MIN; expected $EXPECTED_SEX)"
    elif [ "$OBSERVED_SEX" != "$EXPECTED_SEX" ]; then
        add_fail "SEX MISMATCH: expected $EXPECTED_SEX, observed $OBSERVED_SEX — possible SAMPLE SWAP"
    elif [ "$OBSERVED_SEX" = "M" ] && awk -v y="$YNORM" 'BEGIN{exit !(y+0 < 0.5)}'; then
        notes+=("chrY yield is low ($YNORM/1000) — normal for capture kits, which barely target chrY")
    fi
fi

# ── verdict ───────────────────────────────────────────────────────────────────────────
STATUS=PASS; CODE=0
(( ${#warns[@]} > 0 )) && { STATUS=WARN; CODE=2; }
(( ${#fails[@]} > 0 )) && { STATUS=FAIL; CODE=1; }

printf '== QC gate: %s ==\n' "$SAMPLE"
printf '  variants        %s\n  Ti/Tv           %s\n  het/hom         %s\n  skewed-AB hets  %s\n  mean depth      %sx\n' \
       "$NVAR" "$TITV" "$HETHOM" "$SKEWFRAC" "$MEANDP"
for n in "${notes[@]}"; do printf '  %s\n' "$n"; done
for w in "${warns[@]}"; do printf '  WARN  %s\n' "$w"; done
for f in "${fails[@]}"; do printf '  FAIL  %s\n' "$f"; done
printf '  ---> %s\n' "$STATUS"

if [ -n "$JSON" ]; then
  { printf '{"sample":"%s","status":"%s","variants":%s,"titv":%s,"het_hom":%s,' "$SAMPLE" "$STATUS" "$NVAR" "$TITV" "$HETHOM"
    printf '"skewed_ab_het_frac":%s,"mean_dp":%s,"chrx_het_rate":%s,"chry_calls":%s,' "$SKEWFRAC" "$MEANDP" "$XHET" "$YN"
    printf '"observed_sex":"%s","expected_sex":"%s","failures":%d,"warnings":%d}\n' \
           "$OBSERVED_SEX" "$EXPECTED_SEX" "${#fails[@]}" "${#warns[@]}"; } > "$JSON"
fi
exit $CODE
