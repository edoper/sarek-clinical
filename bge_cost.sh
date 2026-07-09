#!/usr/bin/env bash
# Live COST bar for a BGE Sarek Batch run — reconstructs accrued Google Batch Spot
# spend from the job records (control-plane listing only: NO compute, NO egress).
# Live:   watch -n 60 ~/sarek-clinical/bge_cost.sh
#
# Env:
#   SINCE   ISO8601 lower bound on job createTime (default: 36h ago, i.e. this run).
#           For an exact run set e.g. SINCE=2026-06-14T00:00:00Z
#   BUDGET  budget ceiling in USD for the bar (default 20)
#   REGION  Batch location (default us-central1)
#   FRAC    optional completion fraction 0<f<=1 (from the progress bar) → prints a
#           projected final cost = accrued / FRAC. The dashboard passes this in.
set -euo pipefail
REGION="${REGION:-us-central1}"
BUDGET="${BUDGET:-20}"
SINCE="${SINCE:-$(date -u -d '36 hours ago' +%Y-%m-%dT%H:%M:%SZ)}"
NOW=$(date -u +%s)

# machineType -> "vCPU memGB"  (types seen in this project's runs)
declare -A SPEC=(
  [c2-standard-30]="30 120" [c2-standard-8]="8 32" [c2-standard-4]="4 16"
  [c2d-highcpu-2]="2 4"
  [n1-standard-1]="1 3.75"
  [n2-standard-2]="2 8"
  [n2-highcpu-4]="4 4" [n2-highcpu-8]="8 8" [n2-highcpu-16]="16 16" [n2-highcpu-32]="32 32"
  [n2-highmem-2]="2 16" [n2-highmem-4]="4 32" [n2-highmem-8]="8 64" [n2-highmem-16]="16 128"
  [n2d-highcpu-2]="2 2" [n2d-highcpu-4]="4 4" [n2d-highcpu-8]="8 8" [n2d-highcpu-16]="16 16" [n2d-highcpu-32]="32 32"
  [n2d-standard-2]="2 8" [n2d-standard-4]="4 16" [n2d-standard-8]="8 32" [n2d-standard-16]="16 64" [n2d-standard-32]="32 128"
  [n2d-highmem-2]="2 16" [n2d-highmem-4]="4 32" [n2d-highmem-8]="8 64" [n2d-highmem-16]="16 128"
)
# family -> SPOT USD "perVCPUhr perGBhr" (us-central1, ~early-2026)
declare -A RATE=(
  [n1]="0.006655 0.000892" [n2]="0.007540 0.001010" [n2d]="0.006554 0.000878"
  [c2]="0.007820 0.001047" [c2d]="0.006810 0.000912"
)

csv=$(gcloud batch jobs list --location="$REGION" \
        --filter="createTime>=\"$SINCE\"" \
        --format="csv[no-heading](createTime,allocationPolicy.instances[0].policy.machineType,allocationPolicy.instances[0].policy.provisioningModel,status.state,status.runDuration)" \
        2>/dev/null || true)

cost=0; njobs=0; running=0; failed=0; unknown=""
while IFS=, read -r ctime mt model state dur; do
  [ -z "${mt:-}" ] && continue
  njobs=$((njobs+1))
  case "$state" in RUNNING|SCHEDULED|QUEUED|*RUNNING*) running=$((running+1));; FAILED) failed=$((failed+1));; esac
  spec="${SPEC[$mt]:-}"
  if [ -z "$spec" ]; then unknown="$unknown $mt"; continue; fi
  read -r vcpu memgb <<<"$spec"
  fam="${mt%%-*}"; read -r rc rm <<<"${RATE[$fam]:-0 0}"
  # seconds: prefer reported runDuration; for still-running jobs use elapsed wall time
  sec=0
  if [ -n "${dur:-}" ]; then sec="${dur%s}"; else
    case "$state" in *RUNNING*|SCHEDULED|QUEUED) sec=$(( NOW - $(date -u -d "$ctime" +%s 2>/dev/null || echo "$NOW") ));; esac
  fi
  cost=$(awk -v c="$cost" -v s="$sec" -v v="$vcpu" -v m="$memgb" -v rc="$rc" -v rm="$rm" \
             'BEGIN{printf "%.4f", c + (s/3600.0)*(v*rc + m*rm)}')
done <<< "$csv"

# budget bar (20 cells)
pct=$(awk -v c="$cost" -v b="$BUDGET" 'BEGIN{p=b>0?c*100/b:0; if(p>100)p=100; printf "%d", p}')
fill=$(( pct/5 ))
bar=$(printf '%*s' "$fill" '' | tr ' ' '#')$(printf '%*s' $((20-fill)) '' | tr ' ' '.')

proj=""
if [ -n "${FRAC:-}" ]; then
  proj=$(awk -v c="$cost" -v f="$FRAC" 'BEGIN{ if(f>0 && f<=1) printf "  proj ~$%.2f @100%%", c/f }')
fi
printf "cost Spot [%s] \$%.2f / \$%.0f budget (%d%%) | %d jobs, %d run, %d fail%s\n" \
  "$bar" "$cost" "$BUDGET" "$pct" "$njobs" "$running" "$failed" "$proj"
[ -n "$unknown" ] && printf "  (unpriced machine types:%s — extend SPEC[] in bge_cost.sh)\n" "$unknown"
exit 0
