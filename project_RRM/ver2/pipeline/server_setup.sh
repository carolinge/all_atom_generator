#!/bin/bash
# server_setup.sh
# ===============
# Run ONCE on the bio cluster (Linux) to set up ver_2 server-side env.
# Idempotent — safe to re-run.
#
# Strategy: CLONE the working ver_1 `allatom` env, then add ambertools as
# a single pinned package. Avoids the 30-min "conda solve forever" trap
# that hits when env-create tries to resolve 25 packages from scratch.
#
# Usage:
#     ssh bio
#     cd /data/biophys/carolinge/clawork/37_OXR/repo
#     bash project_RRM/ver2/pipeline/server_setup.sh

set -euo pipefail

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
SOURCE_ENV="allatom"      # ver_1's env, used as baseline
TARGET_ENV="allatom_v2"

# Pin the ambertools build to avoid solver thrashing (Python 3.11 + no GPU
# CUDA + OpenMPI variant). conda-forge build hashes change occasionally;
# the wildcard '=*py311*' lets the solver pick the latest matching build
# without exploring py310/312/313/cuda branches.
AMBERTOOLS_SPEC='ambertools=24.8=*py311*'

# 1. CUDA module (matches the cudatoolkit pinned in allatom)
echo ">>> module load cuda/12.5"
module purge
module load cuda/12.5

# 2. Persistent project sub-dirs
echo ">>> Ensuring $PROJECT/ subdirs"
mkdir -p "$PROJECT"/{replicas_v2,force_fields/lib_amber,structures/{original,prepared,modified},pipeline_logs}

# 3. Conda env
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not in PATH"; exit 1
fi
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "$TARGET_ENV"; then
    echo ">>> Env $TARGET_ENV already exists — skipping clone"
else
    if ! conda env list | awk '{print $1}' | grep -qx "$SOURCE_ENV"; then
        echo "ERROR: source env '$SOURCE_ENV' not found. ver_1 not installed?"
        exit 1
    fi
    echo ">>> Cloning $SOURCE_ENV -> $TARGET_ENV (~30 s, just a file copy)"
    conda create --name "$TARGET_ENV" --clone "$SOURCE_ENV" -y
fi

conda activate "$TARGET_ENV"

# 4. Add ambertools (idempotent — skip if already present)
if ! conda list -n "$TARGET_ENV" 2>/dev/null | grep -q '^ambertools'; then
    echo ">>> Installing $AMBERTOOLS_SPEC into $TARGET_ENV (~5 min)"
    conda install -n "$TARGET_ENV" -c conda-forge -y "$AMBERTOOLS_SPEC"
else
    echo ">>> ambertools already installed; skipping"
fi

# 5. Verify
echo ""
echo "=== Verify AmberTools binaries ==="
for cmd in cpptraj tleap sander; do
    if command -v "$cmd" >/dev/null 2>&1; then
        echo "    OK: $cmd  ($(which "$cmd"))"
    else
        echo "    FAIL: $cmd not found"
        exit 1
    fi
done

echo ""
echo "=== Verify OpenMM platform ==="
python -c "
from openmm import Platform
print('OpenMM platforms:', [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
" || true

echo ""
echo "=== server_setup.sh DONE ==="
echo "Next: bash project_RRM/ver2/pipeline/build_modxna_residues.sh"
