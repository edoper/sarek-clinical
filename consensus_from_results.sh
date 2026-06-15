#!/usr/bin/env bash
#
# consensus_from_results.sh — bridge step between cloud calling and local interpretation.
# For each sample it pulls the four small per-caller VCFs from the bucket and runs the
# canonical consensus.sh (reused from ~/sarek-clinical — NOT duplicated here), producing
# <sample>.consensus.vcf.gz ready for vep_annotate.sh + candidate-filtering.
#
# Consensus is a ~seconds bcftools step on tiny exome VCFs, so it runs locally on the
# downloaded per-caller VCFs (the heavy calling already happened in the cloud). Sample
# names carry the -P/-M/-F convention, so the consensus outputs feed candidate-filtering
# trio/duo auto-discovery directly.
set -euo pipefail
cd "$(dirname "$0")"

SAMPLESHEET="${SAMPLESHEET:-samplesheet.csv}"
OUTDIR="${OUTDIR:-gs://intergenica-sarek-clinical/bge-wes/results}"      # Sarek --outdir (cloud)
LOCAL_OUT="${LOCAL_OUT:-./consensus}"                                    # where consensus VCFs land
CONSENSUS_SH="${CONSENSUS_SH:-$HOME/sarek-clinical/consensus.sh}"        # canonical script (reused)
# Reuse the chr-named GRCh38 already on disk (matches GATK chr naming on the main
# chromosomes, where exome candidates live). Override REF if you prefer the full GATK ref.
REF="${REF:-$HOME/vep_refs/pangolin/GRCh38.primary_assembly.genome.fa}"

[[ -f "$SAMPLESHEET"  ]] || { echo "ERROR: $SAMPLESHEET not found"; exit 1; }
[[ -x "$CONSENSUS_SH" ]] || { echo "ERROR: consensus.sh not found/executable: $CONSENSUS_SH"; exit 1; }
[[ -f "$REF.fai"      ]] || { echo "ERROR: reference index missing: $REF.fai"; exit 1; }
mkdir -p "$LOCAL_OUT"

# Find the published VCF for a caller/sample under the Sarek outdir (filenames vary
# slightly by caller, e.g. strelka.variants / haplotypecaller.filtered).
find_vcf() {  # <caller> <sample> — prefer filtered > variants > plain; never genome/gvcf
    local all f
    all=$(gcloud storage ls "$OUTDIR/variant_calling/$1/$2/**.vcf.gz" 2>/dev/null | grep -viE '\.(g|genome)\.vcf\.gz$')
    for pat in 'filtered\.vcf\.gz$' 'variants\.vcf\.gz$' '\.vcf\.gz$'; do
        f=$(echo "$all" | grep -iE "$pat" | head -1)
        [ -n "$f" ] && { echo "$f"; return; }
    done
}

# Samples = column 2 of the samplesheet (skip header)
samples=$(awk -F',' 'NR>1{print $2}' "$SAMPLESHEET" | sort -u)
[[ -n "$samples" ]] || { echo "ERROR: no samples in $SAMPLESHEET"; exit 1; }

for s in $samples; do
    echo "=== $s ==="
    tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' RETURN
    declare -A vcf=()
    ok=1
    for caller in deepvariant strelka freebayes haplotypecaller; do
        uri=$(find_vcf "$caller" "$s")
        if [[ -z "$uri" ]]; then echo "  WARN: no $caller VCF for $s — skipping sample"; ok=0; break; fi
        gcloud storage cp "$uri" "$uri.tbi" "$tmp/" 2>/dev/null || true
        vcf[$caller]="$tmp/$(basename "$uri")"
    done
    [[ "$ok" == 1 ]] || { rm -rf "$tmp"; continue; }

    "$CONSENSUS_SH" -r "$REF" \
        -d "${vcf[deepvariant]}" \
        -o "$LOCAL_OUT/$s" \
        -c strelka="${vcf[strelka]}" \
        -c freebayes="${vcf[freebayes]}" \
        -c haplotypecaller="${vcf[haplotypecaller]}"
    rm -rf "$tmp"
done

echo
echo "Consensus VCFs in: $LOCAL_OUT/<sample>.consensus.vcf.gz"
echo "Next (in candidate-filtering): vep_annotate.sh <sample>.consensus.vcf.gz <sample>.germline.vep.vcf.gz  ->  filtering_r.pl"
