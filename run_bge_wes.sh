#!/usr/bin/env bash
#
# run_bge_wes.sh — launch nf-core/sarek WES variant-calling (from CRAM) on Google Batch.
# The heavy step (4 callers over the exome target). Per-caller VCFs land in the bucket;
# build the consensus with consensus_from_results.sh afterwards.
#
# Prereqs (one-time): see README. You must set INTERVALS to your exome target BED/
# interval_list (gs:// path); see "Inputs you must provide".
set -euo pipefail
cd "$(dirname "$0")"

source ~/sarek-clinical/env.sh          # nextflow + JDK + NXF_SYNTAX_PARSER=v1

SAMPLESHEET="${SAMPLESHEET:-samplesheet.csv}"
OUTDIR="${OUTDIR:-gs://intergenica-sarek-clinical/bge-wes/results}"
INTERVALS="${INTERVALS:?set INTERVALS to your exome targets, e.g. gs://intergenica-sarek-clinical/bge-wes/targets/bge_calling_regions.bed}"

[[ -f "$SAMPLESHEET" ]] || { echo "ERROR: $SAMPLESHEET not found (run make_samplesheet.sh first)"; exit 1; }

echo "[run_bge_wes] samplesheet=$SAMPLESHEET  intervals=$INTERVALS  outdir=$OUTDIR"
nextflow run nf-core/sarek -r 3.8.1 \
  -profile docker -c gcb-bge-wes.config \
  --step variant_calling \
  --input "$SAMPLESHEET" \
  --intervals "$INTERVALS" \
  --outdir "$OUTDIR" \
  --tools deepvariant,strelka,freebayes,haplotypecaller \
  --genome GATK.GRCh38 \
  -ansi-log false -resume

echo "[run_bge_wes] done. Per-caller VCFs in: $OUTDIR/variant_calling/<caller>/<sample>/"
echo "[run_bge_wes] next: ./consensus_from_results.sh"
