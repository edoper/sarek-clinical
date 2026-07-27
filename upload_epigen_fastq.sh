#!/usr/bin/env bash
# Upload exome FASTQ (local -> bucket) for a Sarek from-FASTQ run.
# GCS ingress is free; gcloud storage cp parallelises automatically. Idempotent
# (skips files already present with matching size via --no-clobber off / rsync-style).
#
#   SRC=/path/to/fastq bash upload_epigen_fastq.sh           # upload
#   SRC=/path/to/fastq CHECK=1 bash upload_epigen_fastq.sh   # just compare local vs bucket counts
set -euo pipefail
. "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/site.sh"
# SRC has no default on purpose: it is your own local FASTQ folder.
SRC="${SRC:?set SRC to the local folder holding your *.fastq.gz (e.g. SRC=/data/run01 $0)}"
DST="${DST:-$SAREK_BUCKET/epigen-wes/fastq}"

local_n=$(find "$SRC" -iname '*.fastq.gz' | wc -l)
echo "local fastq.gz: $local_n   size: $(du -sh "$SRC" | cut -f1)"

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
