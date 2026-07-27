# Shared site settings for every script in this repo. Sourced (not executed), silent.
#   . "$(dirname "${BASH_SOURCE[0]}")/site.sh"
#
# Nothing here is tied to one machine or one person. Every value can be overridden by
# exporting it beforehand, or by creating an untracked `site.env` next to this file:
#
#   # site.env
#   SAREK_PROJECT=my-gcp-project
#   SAREK_BUCKET=gs://my-bucket
#   CF=$HOME/code/candidate-filtering
#   WIN=/mnt/c/Users/me/Documents      # WSL only; leave unset elsewhere
#
# The defaults are simply the deployment this repo was developed against — they are a
# starting point, not a requirement. gcb*.config read the same SAREK_* variables, so
# setting them once retargets both the shell scripts and Nextflow.

# Repo root = wherever this file lives, so any checkout location works.
SAREK_REPO="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export SAREK_REPO

# Untracked per-site overrides, if present.
# Precedence: an explicitly exported value  >  site.env  >  the defaults below. So a
# one-off `SAREK_BUCKET=gs://other ./run_bge_wes.sh` still wins over a committed-to-disk
# site.env, which is what a person (or an agent, or CI) exporting a variable expects.
if [ -f "$SAREK_REPO/site.env" ]; then
    _keys=(); _vals=()
    for _v in SAREK_PROJECT SAREK_REGION SAREK_BUCKET CF WIN; do
        [ -n "${!_v:-}" ] && { _keys+=("$_v"); _vals+=("${!_v}"); }
    done
    . "$SAREK_REPO/site.env"
    for _i in "${!_keys[@]}"; do printf -v "${_keys[$_i]}" '%s' "${_vals[$_i]}"; done
    unset _v _i _keys _vals
fi

# --- Google Cloud target ----------------------------------------------------
: "${SAREK_PROJECT:=intergenica}"                       # GCP project id
: "${SAREK_REGION:=us-central1}"                        # Batch location
: "${SAREK_BUCKET:=gs://intergenica-sarek-clinical}"    # bucket root: fastq/ work/ results/ refs/
export SAREK_PROJECT SAREK_REGION SAREK_BUCKET

# --- Local layout -----------------------------------------------------------
# candidate-filtering: the downstream repo. Defaults to a sibling of this checkout.
: "${CF:=$(dirname "$SAREK_REPO")/candidate-filtering}"
# WIN: Windows-side deliverable folder (WSL convenience only). Deliberately EMPTY by
# default — on Linux/macOS there is no /mnt/c, and scripts skip the copy-out step and
# print where the results are instead of failing.
: "${WIN:=}"
export CF WIN
