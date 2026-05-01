#!/bin/bash
# setup_modxna_server.sh
# ======================
# One-time setup on a Linux server (e.g. the bio cluster) to enable
# running modxna.sh. Installs AmberTools 24+ via conda into the
# allatom_v2 environment.
#
# AmberTools provides: cpptraj, tleap, sander — required by modxna.sh.
#
# Usage:
#     bash project_RRM/ver2/pipeline/setup_modxna_server.sh

set -euo pipefail

ENV_NAME="allatom_v2"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Install miniconda first."
    exit 1
fi

# Source conda for the rest of this script
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Create env if missing
if ! conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo ">>> Creating conda env $ENV_NAME from environment_v2.yml"
    conda env create -f "$(dirname "$0")/../../environment_v2.yml"
fi

conda activate "$ENV_NAME"

# Install ambertools (only on Linux — this is a Linux-only conda package)
echo ">>> Installing ambertools into $ENV_NAME (conda-forge channel)"
conda install -n "$ENV_NAME" -c conda-forge -y ambertools

# Verify
echo ""
echo "=== Verify ==="
which cpptraj && cpptraj --version | head -1
which tleap   && echo "tleap OK"
which sander  && echo "sander OK"

echo ""
echo "Done. Next:"
echo "  bash project_RRM/ver2/pipeline/build_modxna_residues.sh"
