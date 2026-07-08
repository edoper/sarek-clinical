#!/usr/bin/env bash
# One-glance BGE cohort dashboard: cloud calling progress + local VEP/filter progress
# + live Spot cost/budget bar with a projected final cost. All monitoring is
# control-plane/listing only (no compute, no egress).
#
# Live:   watch -n 30 ~/sarek-clinical/bge_dashboard.sh
# Env:    BUDGET (USD, default 20), SINCE (cost window, default 36h ago),
#         plus anything bge_progress.sh / bge_filter_progress.sh / bge_cost.sh accept.
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "── BGE cohort ─ $(date -u '+%Y-%m-%d %H:%M UTC') ──────────────────────────"

# 1) cloud calling — capture so we can extract completion % for the cost projection
calling="$("$HERE/bge_progress.sh" 2>/dev/null || echo 'calling: (progress unavailable)')"
printf '  %s\n' "$calling"

# 2) local annotate/filter stage
printf '  %s\n' "$("$HERE/bge_filter_progress.sh" 2>/dev/null || echo 'filter: (progress unavailable)')"

# 3) cost bar — feed calling % (…(NN%)…) in as the completion fraction for projection
pct=$(sed -n 's/.*(\([0-9]\+\)%).*/\1/p' <<<"$calling" | head -1)
if [ -n "${pct:-}" ] && [ "$pct" -gt 0 ] 2>/dev/null; then
  FRAC=$(awk -v p="$pct" 'BEGIN{printf "%.4f", p/100}')
  printf '  %s\n' "$(FRAC="$FRAC" "$HERE/bge_cost.sh" 2>/dev/null || echo 'cost: (unavailable)')"
else
  printf '  %s\n' "$("$HERE/bge_cost.sh" 2>/dev/null || echo 'cost: (unavailable)')"
fi
