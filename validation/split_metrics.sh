#!/usr/bin/env bash
# split_metrics.sh — per-class (SNV / INDEL) precision, recall and F1 for each vcfeval tier.
#
# RTG's summary.txt reports one combined row. Indel accuracy is always the lower of the two
# and is the number that actually discriminates between pipelines, so split it out from the
# tp / fp / fn VCFs vcfeval already wrote. Counting is exact — no re-running of vcfeval.
#
#   precision = TP_call     / (TP_call     + FP)      # of what we called, how much is real
#   recall    = TP_baseline / (TP_baseline + FN)      # of what is true, how much we found
set -uo pipefail
cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
OUT=${OUT:-results-eval}

n() { # <vcf.gz> <snps|indels>  -> variant count of that class
    [ -s "$1" ] || { echo 0; return; }
    bcftools view -v "$2" "$1" 2>/dev/null | grep -cv '^#' || echo 0
}
row() { # <tier> <class>
    local d="$OUT/$1" c="$2"
    local tpc tpb fp fn
    tpc=$(n "$d/tp.vcf.gz" "$c"); tpb=$(n "$d/tp-baseline.vcf.gz" "$c")
    fp=$(n "$d/fp.vcf.gz" "$c");  fn=$(n "$d/fn.vcf.gz" "$c")
    awk -v t="$1" -v c="$c" -v tpc="$tpc" -v tpb="$tpb" -v fp="$fp" -v fn="$fn" 'BEGIN{
        p = (tpc+fp) ? tpc/(tpc+fp) : 0;
        r = (tpb+fn) ? tpb/(tpb+fn) : 0;
        f = (p+r)    ? 2*p*r/(p+r)  : 0;
        printf "| %-18s | %-5s | %9d | %7d | %7d | %.4f | %.4f | %.4f |\n", t, c, tpb, fp, fn, p, r, f;
    }'
}

echo "| tier | class | TP | FP | FN | precision | recall | F1 |"
echo "|---|---|---|---|---|---|---|---|"
for tier in dvonly union-genomewide ncallers2 confhigh union-exome; do
    [ -d "$OUT/$tier" ] || continue
    row "$tier" snps
    row "$tier" indels
done
