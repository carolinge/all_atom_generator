#!/bin/bash
# server_setup.sh
# ===============
# Run ONCE on the bio cluster (Linux) to set up the ver_2 server-side
# environment. Idempotent — safe to re-run.
#
# Cluster: newton.pks.mpg.de (login as `ssh bio`)
# Project root: /data/biophys/carolinge/clawork/37_OXR
# Repo location after sync: $PROJECT/repo (rsync target of d:/all_atom)
#
# Usage on the server:
#     cd ~/work/37_OXR/repo
#     bash project_RRM/ver2/pipeline/server_setup.sh

set -euo pipefail

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
ENV_NAME="allatom_v2"

# 1. Verify cluster modules
echo ">>> Loading CUDA module"
module purge
module load cuda/12.5

# 2. Make persistent project sub-dirs (idempotent)
echo ">>> Ensuring project directories"
mkdir -p "$PROJECT/replicas_v2"
mkdir -p "$PROJECT/force_fields"
mkdir -p "$PROJECT/structures/original"
mkdir -p "$PROJECT/structures/prepared"
mkdir -p "$PROJECT/structures/modified"
mkdir -p "$PROJECT/pipeline_logs"

# 3. Conda env
echo ">>> Conda env '$ENV_NAME'"
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found in PATH. Install miniconda first."
    exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
ENV_YAML="$REPO_ROOT/project_RRM/environment_v2.yml"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "    Env $ENV_NAME exists. Updating from $ENV_YAML"
    conda env update -n "$ENV_NAME" -f "$ENV_YAML" --prune
else
    echo "    Creating env $ENV_NAME from $ENV_YAML"
    conda env create -f "$ENV_YAML"
fi

conda activate "$ENV_NAME"

# 4. Verify the AmberTools components we need
echo ""
echo ">>> Verify AmberTools binaries"
for cmd in cpptraj tleap sander; do
    if command -v "$cmd" >/dev/null 2>&1; then
        echo "    OK: $cmd  ($(which "$cmd"))"
    else
        echo "    FAIL: $cmd not found in env. AmberTools install incomplete."
        exit 1
    fi
done

# 5. Verify OpenMM CUDA visibility
echo ""
echo ">>> Verify OpenMM CUDA platform"
python -c "
from openmm import Platform
print(f'OpenMM platforms: {[Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]}')
" || true

echo ""
echo "=== server_setup.sh DONE ==="
echo "Next: bash project_RRM/ver2/pipeline/build_modxna_residues.sh"
