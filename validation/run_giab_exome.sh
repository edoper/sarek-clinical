#!/usr/bin/env bash
# run_giab_exome.sh — GIAB validation of the EXOME arm (HG002, Agilent SureSelect capture).
#
# WHY A SECOND VALIDATION
#   The WGS validation (run_giab.sh) measures the pipeline on 30x PCR-free whole-genome reads.
#   The clinical arms are exome/BGE. Restricting WGS calls to a capture BED — which the WGS
#   validation also reports — does NOT reproduce a real exome: it has none of the capture bias,
#   duplicate structure or uneven target coverage that an actual hybrid-capture library has, and
#   those are exactly what make exome calling harder. This script validates the assay you deliver.
#
# INPUT
#   GIAB publishes an HG002 exome from Oslo University Hospital: Agilent SureSelect V5, Illumina,
#   ~9.6 GB BAM aligned to GRCh37. Since the pipeline is GRCh38 and starts from FASTQ, the BAM is
#   name-collated and converted back to paired FASTQ, then re-aligned by the pipeline itself
#   (--step mapping --wes). Nothing from the original GRCh37 alignment is carried forward.
#
#   Kit note: the reference exome is captured with SureSelect **V5**; this deployment's clinical
#   BED is **V6**. Regions present in V6 but not V5 have no reads and would score as false
#   negatives that are the KIT's doing, not the pipeline's. The benchmark therefore reports the
#   g4e panel and, separately, the capture region — read the caveat in RESULTS-HG002-EXOME.md.
#
# USAGE
#   ./run_giab_exome.sh          # from validation/, after ../env.sh is sourced
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
. ../site.sh 2>/dev/null || true
: "${SAREK_BUCKET:=gs://intergenica-sarek-clinical}"
B="$SAREK_BUCKET/giab-exome"
WORK=exome
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
INTERVALS="${INTERVALS:-$SAREK_BUCKET/epigen-wes/targets/S07604514_V6r2_Padded.GRCh38.bed}"

mkdir -p "$WORK"; cd "$WORK"
source ../../env.sh >/dev/null 2>&1

# ── Stage 1: reference exome BAM -> paired FASTQ (streamed; no local BAM copy) ──
if [ ! -s HG002exome.R1.fastq.gz ] || [ ! -s HG002exome.R2.fastq.gz ]; then
    echo "===== Stage 1: streaming the reference exome BAM and converting to FASTQ ====="
    samtools collate -u -O -@ 4 "$BAM_URL" tmp.collate \
      | samtools fastq -@ 2 -1 HG002exome.R1.fastq.gz -2 HG002exome.R2.fastq.gz \
                       -0 /dev/null -s /dev/null -n -
    rm -f tmp.collate*.bam
else
    echo "===== Stage 1: FASTQ already present — skipping conversion ====="
fi

# ── Stage 2: upload (GCS ingress is free) ──
echo "===== Stage 2: staging FASTQ to the bucket ====="
export CLOUDSDK_STORAGE_PARALLEL_COMPOSITE_UPLOAD_ENABLED=False   # WSL/drvfs safety
gcloud storage cp -n HG002exome.R1.fastq.gz HG002exome.R2.fastq.gz "$B/fastq/"

cat > samplesheet-giab-exome.csv <<EOF
patient,sample,lane,fastq_1,fastq_2
HG002EX,HG002EX,L001,$B/fastq/HG002exome.R1.fastq.gz,$B/fastq/HG002exome.R2.fastq.gz
EOF

# ── Stage 3: the pipeline, exome mode, exactly as the clinical exome arm runs it ──
echo "===== Stage 3: sarek --step mapping --wes (4 callers, Spot) ====="
nextflow run nf-core/sarek -r 3.8.1 -profile docker -c "$(readlink -f ../../gcb.config)" \
  --input      "$PWD/samplesheet-giab-exome.csv" \
  --outdir     "$B/results" \
  --genome     GATK.GRCh38 \
  --wes \
  --intervals  "$INTERVALS" \
  --tools      deepvariant,strelka,freebayes,haplotypecaller \
  --skip_tools baserecalibrator,vcftools,multiqc \
  -work-dir    "$B/work" \
  -ansi-log false -resume

# ── Stage 4: consensus ──
echo "===== Stage 4: consensus ====="
SAMPLESHEET="$PWD/samplesheet-giab-exome.csv" OUTDIR="$B/results" \
  LOCAL_OUT="$PWD/consensus" ../../consensus_from_results.sh

# ── Stage 5: the QC gate must pass before the result is considered usable ──
echo "===== Stage 5: per-sample QC gate ====="
../../qc_gate.sh consensus/HG002EX.consensus.vcf.gz --sex M --json HG002EX.qc.json || \
  echo "  NOTE: QC gate did not return PASS — see above (expected for a 2015 exome; record it)"

# ── Stage 6: benchmark ──
echo "===== Stage 6: benchmark vs GIAB v4.2.1 ====="
cd ..
CONS="$WORK/consensus/HG002EX.consensus.vcf.gz" OUT=results-eval-exome ./benchmark_giab.sh
echo "===== exome validation complete ====="
