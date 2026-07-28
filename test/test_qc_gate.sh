#!/usr/bin/env bash
# test_qc_gate.sh — regression test for the per-sample QC gate.
#
# A gate that never fails is worse than no gate, and a gate that fails GOOD samples is
# worse still: it trains people to ignore it. Both directions are tested here.
#
# The sex checks exist because the first version of qc_gate.sh FAILED a known-male exome
# (2026-07-28). Two WGS-calibrated assumptions were at fault:
#   (a) pseudoautosomal regions were included, inflating a male's chrX het rate 0.090 -> 0.169
#   (b) an ABSOLUTE chrY call threshold was used; the same sample gives 11,375 chrY calls by
#       WGS and 57 by exome, so any absolute cutoff silently means "WGS only"
#
# Needs a real consensus VCF to work from; skips cleanly if none is present.
set -uo pipefail
REPO="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)"
GATE="$REPO/qc_gate.sh"
fails=0
ok(){ echo "  PASS  $*"; }
bad(){ echo "  FAIL  $*"; fails=$((fails+1)); }

# find any male sample VCF produced by the validation runs
VCF=""
for c in "$REPO"/validation/exome/consensus/*.consensus.vcf.gz \
         "$REPO"/giab-validation/consensus/*.consensus.vcf.gz; do
    [ -s "$c" ] && { VCF="$c"; break; }
done
[ -n "$VCF" ] || { echo "  SKIP  no consensus VCF available (run validation/ first)"; exit 0; }
echo "== QC gate regression ($(basename "$VCF")) =="

"$GATE" "$VCF" --sex M >/dev/null 2>&1; rc=$?
[ $rc -eq 0 ] && ok "known-male sample declared M -> PASS" || bad "known-male declared M returned exit $rc"
"$GATE" "$VCF" --sex F >/dev/null 2>&1; rc=$?
[ $rc -eq 1 ] && ok "known-male declared F -> FAIL (swap detection)" || bad "sex swap not caught (exit $rc)"

TD="$(mktemp -d)"; trap 'rm -rf "$TD"' EXIT
# synthetic female: force chrX non-PAR genotypes heterozygous
bcftools view -h "$VCF" > "$TD/f.vcf"
bcftools view -H "$VCF" | awk -F'\t' 'BEGIN{OFS="\t"}
  $1=="chrX" && $2>2781479 && $2<155701383 {sub(/^[01][\/|][01]/,"0/1",$10)} {print}' >> "$TD/f.vcf"
bgzip -f "$TD/f.vcf" && tabix -f -p vcf "$TD/f.vcf.gz"
"$GATE" "$TD/f.vcf.gz" --sex F >/dev/null 2>&1; rc=$?
[ $rc -eq 0 ] && ok "synthetic female declared F -> PASS" || bad "female not recognised (exit $rc)"
"$GATE" "$TD/f.vcf.gz" --sex M >/dev/null 2>&1; rc=$?
[ $rc -eq 1 ] && ok "synthetic female declared M -> FAIL" || bad "reverse swap not caught (exit $rc)"

# collapsed library
bcftools view -h "$VCF" > "$TD/tiny.vcf"; bcftools view -H "$VCF" | head -5000 >> "$TD/tiny.vcf"
bgzip -f "$TD/tiny.vcf" && tabix -f -p vcf "$TD/tiny.vcf.gz"
"$GATE" "$TD/tiny.vcf.gz" >/dev/null 2>&1; rc=$?
[ $rc -eq 1 ] && ok "collapsed library (5k variants) -> FAIL" || bad "low variant count not caught (exit $rc)"

# noise-dominated call set (Ti/Tv collapses)
bcftools view -h "$VCF" > "$TD/noise.vcf"
bcftools view -H "$VCF" | head -60000 | awk 'BEGIN{OFS="\t";b["A"]="C";b["C"]="A";b["G"]="T";b["T"]="G"}
  length($4)==1&&length($5)==1{$5=b[$4]; print}' >> "$TD/noise.vcf"
bgzip -f "$TD/noise.vcf" && tabix -f -p vcf "$TD/noise.vcf.gz"
"$GATE" "$TD/noise.vcf.gz" >/dev/null 2>&1; rc=$?
[ $rc -eq 1 ] && ok "noise-dominated call set (Ti/Tv ~0) -> FAIL" || bad "noise not caught (exit $rc)"

echo
[ $fails -eq 0 ] && echo "ALL TESTS PASSED" || { echo "$fails TEST(S) FAILED"; exit 1; }
