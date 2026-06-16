#!/usr/bin/env bash
# Zero-cost progress bar for the BGE annotate+filter stage. Live:  watch -n 30 ~/sarek-clinical/bge_filter_progress.sh
CONS="${CONS_DIR:-$HOME/sarek-clinical/consensus-cohort}"
WD="${WD:-$HOME/candidate-filtering/bge-cohort}"
total=$(ls "$CONS"/*.consensus.vcf.gz 2>/dev/null | wc -l); [ "$total" -gt 0 ] || total=63
vep=$(ls "$WD"/*.germline.vep.vcf.gz 2>/dev/null | wc -l)
cand=$(ls "$WD"/*.candidatos 2>/dev/null | wc -l)
running=$(pgrep -f "vep|run_bge_annotate_filter|run_filtering" >/dev/null && echo yes || echo no)
pct=$(( total ? vep*100/total : 0 )); fill=$(( pct/5 ))
bar=$(printf '%*s' "$fill" '' | tr ' ' '#')$(printf '%*s' $((20-fill)) '' | tr ' ' '.')
printf "VEP [%s] %d/%d (%d%%) | candidatos: %d | active: %s\n" "$bar" "$vep" "$total" "$pct" "$cand" "$running"
