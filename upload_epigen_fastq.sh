#!/usr/bin/env bash
# Upload the 20 EPIGEN exome FASTQ (local -> bucket) for the Sarek from-FASTQ run.
# GCS ingress is free; gcloud storage cp parallelises automatically. Idempotent
# (skips files already present with matching size via --no-clobber off / rsync-style).
#
#   bash ~/sarek-clinical/upload_epigen_fastq.sh          # upload
#   CHECK=1 bash ~/sarek-clinical/upload_epigen_fastq.sh  # just compare local vs bucket counts
set -euo pipefail
SRC="${SRC:-/mnt/c/Users/epere/Documents/UDD/Acceso-UNICAMP/Data_LNGC/Fastq_files_EPIGEN}"
DST="${DST:-gs://intergenica-sarek-clinical/epigen-wes/fastq}"

local_n=$(find "$SRC" -iname '*.fastq.gz' | wc -l)
echo "local fastq.gz: $local_n   (expect 56)   size: $(du -sh "$SRC" | cut -f1)"

if [ "${CHECK:-0}" = "1" ]; then
  echo "bucket fastq.gz: $(gcloud storage ls "$DST/**.fastq.gz" 2>/dev/null | wc -l)"
  exit 0
fi

# Disable parallel composite uploads: splitting big files into parallel chunks read over
# the WSL /mnt/c (drvfs) bridge corrupts components ("Temporary components were not uploaded
# correctly"). Single-stream per file is reliable. Scoped via env var (no global config change).
export CLOUDSDK_STORAGE_PARALLEL_COMPOSITE_UPLOAD_ENABLED=False

# rsync = resumable + skips already-uploaded identical files (safe to re-run after an interruption)
gcloud storage rsync -r \
  --exclude '.*\.(txt|md5)$' \
  "$SRC" "$DST"

echo "--- verify ---"
echo "bucket fastq.gz: $(gcloud storage ls "$DST/**.fastq.gz" 2>/dev/null | wc -l) / $local_n"
