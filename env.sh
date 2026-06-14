# Source this to load the Sarek/Nextflow toolchain into your shell:
#   source ~/sarek-clinical/env.sh
export JAVA_HOME="$HOME/jdk21"
export PATH="$HOME/jdk21/bin:$HOME/bin:$PATH"
echo "Sarek toolchain loaded: java=$(java -version 2>&1 | head -1), nextflow=$(nextflow -version 2>&1 | grep -m1 version | tr -s ' ')"
