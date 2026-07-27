#!/usr/bin/env bash
# benchmark_giab.sh — compare the consensus VCF against the GIAB HG002 v4.2.1 truth set
# with RTG vcfeval, restricted to GIAB's high-confidence regions.
#
# Evaluates several CONFIDENCE TIERS of the same consensus file, because that is the
# decision the pipeline defers downstream: the union is most sensitive, NCALLERS>=2 and
# CONF=HIGH trade recall for precision. The point of the table is to make that trade
# measurable instead of assumed.
#
# Also reports an EXOME-RESTRICTED row, since the clinical arms are exome/BGE: the same
# calls intersected with the capture BED say more about production than a genome-wide
# number does.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

RTG=./rtg-tools-3.12.1/rtg
SDF=GRCh38.sdf
TRUTH=truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
CONF_BED=truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
CONS="${CONS:-consensus/HG002.consensus.vcf.gz}"
EXOME_BED="${EXOME_BED:-../refs/S07604514_V6r2_Padded.GRCh38.bed}"
OUT=results-eval; mkdir -p "$OUT"

for f in "$RTG" "$TRUTH" "$CONF_BED" "$CONS"; do
    [ -e "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done
[ -d "$SDF" ] || { echo "ERROR: missing SDF $SDF (run: $RTG format -o $SDF <ref.fasta>)"; exit 1; }

# vcfeval needs a single sample and PASS-only comparison semantics we control ourselves.
run_eval() { # <label> <query.vcf.gz> [extra-bed]
    # NOTE: separate statements — `local` expands ALL its arguments before assigning
    # any of them, so `dir="$OUT/$label"` on the same line sees $label unbound under set -u.
    local label="$1" q="$2" bed="${3:-$CONF_BED}"
    local dir="$OUT/$label"
    rm -rf "$dir"
    $RTG vcfeval -b "$TRUTH" -c "$q" -t "$SDF" -e "$bed" -o "$dir" \
         --ref-overlap --all-records >"$OUT/$label.log" 2>&1 || {
        echo "  WARN: vcfeval failed for $label (see $OUT/$label.log)"; return 1; }
    # rtg writes summary.txt with SNP/INDEL/total rows
    echo "  [$label]"; sed 's/^/     /' "$dir/summary.txt" 2>/dev/null | head -6
}

echo "== building confidence tiers from the consensus VCF =="
bcftools view -i 'NCALLERS>=2' "$CONS" -Oz -o "$OUT/tier.ncallers2.vcf.gz" && tabix -f -p vcf "$OUT/tier.ncallers2.vcf.gz"
bcftools view -i 'CONF="HIGH"'  "$CONS" -Oz -o "$OUT/tier.confhigh.vcf.gz"  && tabix -f -p vcf "$OUT/tier.confhigh.vcf.gz"
printf "  union: %s | NCALLERS>=2: %s | CONF=HIGH: %s\n" \
  "$(bcftools index -n "$CONS")" "$(bcftools index -n "$OUT/tier.ncallers2.vcf.gz")" "$(bcftools index -n "$OUT/tier.confhigh.vcf.gz")"

echo
echo "== genome-wide, inside GIAB high-confidence regions =="
run_eval union-genomewide  "$CONS"                      || true
run_eval ncallers2         "$OUT/tier.ncallers2.vcf.gz" || true
run_eval confhigh          "$OUT/tier.confhigh.vcf.gz"  || true

if [ -s "$EXOME_BED" ]; then
    echo
    echo "== exome-restricted (capture BED ∩ GIAB high-confidence) =="
    # intersect the two BEDs so vcfeval evaluates only the exome's confident part
    sort -k1,1 -k2,2n "$CONF_BED" > "$OUT/conf.sorted.bed"
    sort -k1,1 -k2,2n "$EXOME_BED" > "$OUT/exome.sorted.bed"
    if command -v bedtools >/dev/null; then
        bedtools intersect -a "$OUT/conf.sorted.bed" -b "$OUT/exome.sorted.bed" > "$OUT/conf_exome.bed"
    else
        # bedtools-free intersection via awk on sorted interval lists
        awk 'NR==FNR{s[FNR]=$1"\t"$2"\t"$3; n=FNR; next}
             {for(i=1;i<=n;i++){split(s[i],a,"\t");
                if(a[1]==$1 && a[2]<$3 && $2<a[3]){st=(a[2]>$2?a[2]:$2); en=(a[3]<$3?a[3]:$3); print $1"\t"st"\t"en}}}' \
            "$OUT/exome.sorted.bed" "$OUT/conf.sorted.bed" > "$OUT/conf_exome.bed"
    fi
    echo "  exome ∩ high-conf regions: $(wc -l < "$OUT/conf_exome.bed")"
    run_eval union-exome "$CONS" "$OUT/conf_exome.bed" || true
else
    echo "  (no exome BED at $EXOME_BED — skipping the exome-restricted row)"
fi

echo
echo "Full per-tier output under $OUT/<tier>/ (summary.txt, roc curves, fp/fn VCFs)."
