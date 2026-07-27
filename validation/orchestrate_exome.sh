#!/usr/bin/env bash
# Wait for the streaming BAM->FASTQ conversion, then run the exome validation end to end.
# Fail-closed: any stage failure stops the chain and is recorded in the log.
set -uo pipefail
cd "$(dirname "$(readlink -f "$0")")"
echo "[$(date -u +%H:%M:%SZ)] waiting for BAM->FASTQ conversion"
while pgrep -f 'samtools collate' >/dev/null; do sleep 60; done
grep -q CONVERT_OK ../giab-validation/exome/convert.log 2>/dev/null \
  || { echo "[$(date -u +%H:%M:%SZ)] FATAL: conversion did not report success"; exit 1; }
ls -la ../giab-validation/exome/*.fastq.gz
echo "[$(date -u +%H:%M:%SZ)] conversion done -> launching exome validation"
exec ./run_giab_exome.sh
