# Source this to load the Sarek/Nextflow toolchain plus this deployment's settings:
#   source /path/to/sarek-clinical/env.sh
#
# Site settings (GCP project, bucket, local paths) live in site.sh — override them by
# exporting first, or in an untracked site.env next to it. See site.sh for the list.
# The Nextflow configs read the same SAREK_* variables, so one place retargets everything.

. "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/site.sh"

# --- Toolchain (override JAVA_HOME / NXF_HOME if yours live elsewhere) ------
export JAVA_HOME="${JAVA_HOME:-$HOME/jdk21}"
export PATH="$JAVA_HOME/bin:$HOME/bin:$PATH"
# Required: sarek 3.8.1's configs use legacy syntax that Nextflow 26.x's new
# parser rejects ("Unexpected input: ':'"). Force the v1 config parser.
export NXF_SYNTAX_PARSER=v1
echo "Sarek toolchain loaded: java=$(java -version 2>&1 | head -1), nextflow=$(nextflow -version 2>&1 | grep -m1 version | tr -s ' ')"
echo "  repo=$SAREK_REPO  project=$SAREK_PROJECT  bucket=$SAREK_BUCKET${WIN:+  win=$WIN}"
