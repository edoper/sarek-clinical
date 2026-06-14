#!/usr/bin/env bash
#
# make_samplesheet.sh — build an nf-core/sarek (--step variant_calling) samplesheet
# from a simple family/CRAM table, naming samples with the -P/-M/-F convention so the
# downstream candidate-filtering pipeline auto-discovers trios/duos.
#
# INPUT  (TSV, 1 row per individual; '#' comments and blanks ignored):
#     family   role   sex    cram_uri                         [crai_uri]
#   - family : family/case id (becomes Sarek `patient`)
#   - role   : P | M | F   (proband / mother / father)
#   - sex    : XX | XY | NA
#   - cram_uri : gs:// path to the CRAM (Terra workspace bucket or intergenica bucket)
#   - crai_uri : optional; defaults to <cram_uri>.crai
#
# OUTPUT (CSV, Sarek variant_calling schema):
#     patient,sample,sex,status,cram,crai
#   sample = <family>-<role>  (e.g. FAM01-P), status = 0 (germline normal).
#
# Usage:  ./make_samplesheet.sh families.tsv > samplesheet.csv
set -euo pipefail

IN="${1:?usage: make_samplesheet.sh families.tsv > samplesheet.csv}"
[[ -f "$IN" ]] || { echo "ERROR: not found: $IN" >&2; exit 1; }

echo "patient,sample,sex,status,cram,crai"
awk -F'\t' '
  /^[[:space:]]*#/ || NF==0 { next }
  { gsub(/^[[:space:]]+|[[:space:]]+$/,"",$1) }
  $1=="" { next }
  {
    fam=$1; role=$2; sex=($3==""?"NA":$3); cram=$4; crai=$5
    if (role!="P" && role!="M" && role!="F") { print "ERROR: bad role \""role"\" (need P/M/F) for "fam > "/dev/stderr"; exit 1 }
    if (cram !~ /^gs:\/\//) { print "ERROR: cram must be a gs:// URI: "cram > "/dev/stderr"; exit 1 }
    if (crai=="") crai=cram".crai"
    printf "%s,%s-%s,%s,0,%s,%s\n", fam, fam, role, sex, cram, crai
  }
' "$IN"
