#!/usr/bin/env bash
# Zero-cost progress for the BGE Sarek cohort run: reads the local Nextflow log and
# lists the result bucket (object listing only — no compute/egress). Run once, or live:
#     watch -n 60 ~/sarek-clinical/bge_progress.sh
LOG="${LOG:-$HOME/sarek-clinical/bge-cohort.log}"
RESULTS="${RESULTS:-gs://intergenica-sarek-clinical/bge-wes/results-cohort}"
SHEET="${SHEET:-$HOME/sarek-clinical/samplesheet-cohort.csv}"
CALLERS=(deepvariant strelka freebayes haplotypecaller)
mapfile -t SAMPLES < <(awk -F',' 'NR>1{print $2}' "$SHEET" | sort -u)
TOTAL=$(( ${#SAMPLES[@]} * ${#CALLERS[@]} ))

alive=$(pgrep -f "java.*nextflow" >/dev/null && echo "RUNNING" || echo "stopped")
listing=$(gcloud storage ls "$RESULTS/variant_calling/**" 2>/dev/null | grep -iE '\.vcf\.gz$' | grep -viE '\.(g|genome)\.vcf\.gz$')
done=$(echo "$listing" | grep -ciE '/(deepvariant|strelka|freebayes|haplotypecaller)/' )

pct=$(( TOTAL ? done*100/TOTAL : 0 )); filled=$(( pct/5 ))
bar=$(printf '%*s' "$filled" '' | tr ' ' '#')$(printf '%*s' $((20-filled)) '' | tr ' ' '.')
last=$(grep -E 'Submitted process|Pipeline completed|ERROR|Staging foreign' "$LOG" 2>/dev/null | tail -1 | sed 's/.*SAREK://; s/ (.*//' | cut -c1-46)
[ -z "$last" ] && last="(starting / staging)"
grep -q "Pipeline completed successfully" "$LOG" 2>/dev/null && last="COMPLETED ✅"
printf "nextflow:%-8s | %d samples | caller VCFs [%s] %d/%d (%d%%) | %s\n" "$alive" "${#SAMPLES[@]}" "$bar" "$done" "$TOTAL" "$pct" "$last"
