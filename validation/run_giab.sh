#!/usr/bin/env bash
# GIAB validation of the sarek-clinical pipeline — HG002 (NA24385), GRCh38.
#
# Runs the pipeline EXACTLY as it is run clinically (same config, same four callers,
# same --skip_tools) so the accuracy numbers describe the real pipeline, not a
# special-cased variant of it. Then: consensus.sh -> RTG vcfeval vs GIAB v4.2.1.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
. ../site.sh 2>/dev/null || true
: "${SAREK_BUCKET:=gs://intergenica-sarek-clinical}"
B="$SAREK_BUCKET/giab"

source ../env.sh >/dev/null 2>&1

echo "===== Stage 1: sarek WGS from FASTQ (4 callers, Spot) ====="
nextflow run nf-core/sarek -r 3.8.1 -profile docker -c "$(readlink -f ../gcb.config)" \
  --input      "$PWD/samplesheet-giab.csv" \
  --outdir     "$B/results" \
  --genome     GATK.GRCh38 \
  --tools      deepvariant,strelka,freebayes,haplotypecaller \
  --skip_tools baserecalibrator,vcftools,multiqc \
  -work-dir    "$B/work" \
  -ansi-log false -resume

echo "===== Stage 2: consensus ====="
SAMPLESHEET="$PWD/samplesheet-giab.csv" OUTDIR="$B/results" \
  LOCAL_OUT="$PWD/consensus" ../consensus_from_results.sh

echo "===== Stage 3: benchmark vs GIAB v4.2.1 ====="
./benchmark_giab.sh
echo "===== GIAB validation complete ====="
