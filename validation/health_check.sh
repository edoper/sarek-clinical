#!/usr/bin/env bash
# health_check.sh — one hourly status snapshot of the GIAB validation run.
# Zero cost: only lists Batch/Compute control-plane state (no compute, no egress).
#
# Writes STATUS.md locally and cp -f's it to $WIN (never append — drvfs append-caching
# hides live progress from the Windows side).
#
# HARD BUDGET GUARD: the user authorised <= $15 total. If reconstructed spend crosses
# KILL_AT, this cancels the Nextflow driver and every live Batch job rather than letting
# an unattended run walk past the ceiling.
set -uo pipefail
cd "$(dirname "$(readlink -f "$0")")"
. ../site.sh 2>/dev/null || true

SINCE="$(cat run_since.txt 2>/dev/null || date -u -d '6 hours ago' +%Y-%m-%dT%H:%M:%SZ)"
BUDGET="${BUDGET:-15}"
KILL_AT="${KILL_AT:-12.00}"
WINDIR="${WIN:-}/GIAB-validation"
STATUS=STATUS.md

vms=$(timeout 90 gcloud compute instances list --format='value(name)' 2>/dev/null | wc -l)
# Count CLIENT-side: Google Batch's server-side status filter is unreliable (it has
# returned 0 then 5 a minute apart when the truth was 37). One listing, counted locally.
states=$(timeout 120 gcloud batch jobs list --location=us-central1 \
   --filter="createTime>=\"$SINCE\"" --format='value(status.state)' 2>/dev/null)
jobs_running=$(printf '%s\n' "$states" | grep -c '^RUNNING$'   || true)
jobs_done=$(printf '%s\n'    "$states" | grep -c '^SUCCEEDED$' || true)
jobs_failed=$(printf '%s\n'  "$states" | grep -c '^FAILED$'    || true)

cost_out=$(SINCE="$SINCE" BUDGET="$BUDGET" timeout 120 ../bge_cost.sh 2>/dev/null | tail -6)
spend=$(echo "$cost_out" | grep -oE '\$[0-9]+\.[0-9]+' | head -1 | tr -d '$')
spend="${spend:-0.00}"

driver=$(pgrep -f 'nextflow.*sarek' >/dev/null 2>&1 && echo RUNNING || echo "not running")
completed=$(grep -cE '^\[' sarek.log 2>/dev/null || true)
errs=$(grep -ciE 'ERROR ~|Execution aborted|Pipeline completed with errors' sarek.log 2>/dev/null || true)
finished=$(grep -c 'Pipeline completed successfully' sarek.log 2>/dev/null || true)

{
  echo "# GIAB validation — live status"
  echo
  echo "_updated $(date -u '+%Y-%m-%d %H:%M:%SZ') (checks run hourly)_"
  echo
  echo "| | |"
  echo "|---|---|"
  echo "| Nextflow driver | $driver |"
  echo "| Tasks submitted | $completed |"
  echo "| Batch jobs | $jobs_done succeeded / $jobs_running running / $jobs_failed failed |"
  echo "| VMs up now | $vms |"
  echo "| Spend so far | \$$spend of \$$BUDGET authorised |"
  echo "| Driver errors | $errs |"
  echo "| Pipeline finished | $([ "$finished" -gt 0 ] && echo yes || echo 'not yet') |"
  echo
  echo '```'
  echo "$cost_out"
  echo '```'
  echo
  echo "Recent pipeline lines:"
  echo '```'
  grep -E '^\[|Pipeline completed|ERROR' sarek.log 2>/dev/null | tail -8
  echo '```'
} > "$STATUS"

[ -n "${WIN:-}" ] && { mkdir -p "$WINDIR" 2>/dev/null; cp -f "$STATUS" "$WINDIR/STATUS.md" 2>/dev/null; }

# --- hard budget guard -------------------------------------------------------
over=$(awk -v s="$spend" -v k="$KILL_AT" 'BEGIN{print (s+0 > k+0) ? 1 : 0}')
if [ "$over" = "1" ]; then
    echo "BUDGET GUARD TRIPPED: \$$spend > \$$KILL_AT — cancelling run" | tee -a guard.log
    for p in $(pgrep -f 'nextflow.*sarek' 2>/dev/null); do kill "$p" 2>/dev/null; done
    for j in $(timeout 90 gcloud batch jobs list --location=us-central1 \
                 --filter="createTime>=\"$SINCE\" AND status.state=RUNNING" \
                 --format='value(name)' 2>/dev/null); do
        timeout 60 gcloud batch jobs delete "$j" --location=us-central1 --quiet 2>/dev/null
    done
    echo "GUARD: cancelled at \$$spend" >> "$STATUS"
    [ -n "${WIN:-}" ] && cp -f "$STATUS" "$WINDIR/STATUS.md" 2>/dev/null
    exit 2
fi
printf '[%s] driver=%s tasks=%s jobs=%s/%s/%s vms=%s spend=$%s\n' \
   "$(date -u +%H:%M:%SZ)" "$driver" "$completed" "$jobs_done" "$jobs_running" "$jobs_failed" "$vms" "$spend" \
   | tee -a health.log
