#!/usr/bin/env bash
#
# run_epigen_wes.sh — launch nf-core/sarek for the 20 EPIGEN exomes FROM FASTQ on
# Google Batch (Spot). Full path: align (BWA) -> markdup -> [BQSR skipped] -> 4-caller
# variant calling over the exome target. Per-caller VCFs land in the bucket; build the
# consensus with consensus_from_results.sh afterwards, then candidate-filtering locally.
#
# Cost lever notes (cheapest efficient path, ~$25-40 for 20 exomes):
#   * Spot VMs (gcb.config), exome-scoped (--wes --intervals) -> minimal compute
#   * --skip_tools baserecalibrator : BQSR adds cost/time and is unnecessary with a
#     DeepVariant consensus backbone (DeepVariant guidance = no BQSR). Drop this flag
#     if you want strict GATK-best-practice parity with the original LNGC calls.
#   * genome=GATK.GRCh38 via igenomes provides the BWA index + known sites for mapping.
#   * Delete $WORKDIR + the uploaded FASTQ after results are in, to stop storage cost.
set -euo pipefail
cd "$(dirname "$0")"

source ~/sarek-clinical/env.sh          # nextflow + JDK + NXF_SYNTAX_PARSER=v1

SAMPLESHEET="${SAMPLESHEET:-samplesheet-epigen.csv}"
OUTDIR="${OUTDIR:-gs://intergenica-sarek-clinical/epigen-wes/results}"
WORKDIR="${WORKDIR:-gs://intergenica-sarek-clinical/epigen-wes/work}"
# GRCh38, chr-prefixed. Default = the exact Agilent SureSelect V6 r2 (S07604514, hg38)
# Padded BED — kit-exact for the EPIGEN capture. 187k intervals / 100.8 Mb (padded).
# Fallback: gs://intergenica-sarek-clinical/bge-wes/targets/twist_coding_targets.bed (Twist coding).
INTERVALS="${INTERVALS:-gs://intergenica-sarek-clinical/epigen-wes/targets/S07604514_V6r2_Padded.GRCh38.bed}"  # Agilent SureSelect V6 r2 (kit-exact)

[[ -f "$SAMPLESHEET" ]] || { echo "ERROR: $SAMPLESHEET not found (run the prep first)"; exit 1; }

echo "[run_epigen_wes] samplesheet=$SAMPLESHEET"
echo "[run_epigen_wes] intervals=$INTERVALS"
echo "[run_epigen_wes] outdir=$OUTDIR  workdir=$WORKDIR"
nextflow run nf-core/sarek -r 3.8.1 \
  -profile docker -c gcb.config -c epigen-throttle.config \
  --step mapping \
  --wes \
  --input "$SAMPLESHEET" \
  --intervals "$INTERVALS" \
  --skip_tools baserecalibrator,vcftools,multiqc \
  --outdir "$OUTDIR" \
  -work-dir "$WORKDIR" \
  --tools deepvariant,strelka,freebayes,haplotypecaller \
  --genome GATK.GRCh38 \
  -ansi-log false -resume

echo "[run_epigen_wes] done. Per-caller VCFs in: $OUTDIR/variant_calling/<caller>/<sample>/"
echo "[run_epigen_wes] next: SAMPLESHEET=$SAMPLESHEET OUTDIR=$OUTDIR ./consensus_from_results.sh"
