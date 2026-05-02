#!/bin/bash
# server_setup.sh
# ===============
# Run ONCE on the bio cluster (Linux) to set up ver_2 server-side env.
# Idempotent — safe to re-run.
#
# Two conda envs (separation of concerns):
#   allatom_v2  -- OpenMM + ParmEd + PDBFixer + analysis. Cloned from
#                  ver_1's `allatom` env (already validated). Used by
#                  run_md.py and Windows-side prepare_system.py.
#   amber24     -- ambertools (tleap, cpptraj, sander) + python 3.11.
#                  Used only by build_amber_system.sh on the cluster.
#
# Why two envs?  Trying to install ambertools INTO allatom_v2 with the
# stock conda 2022.10 solver hits "phantom __glibc conflict" pseudo-
# failures and 30-minute solve thrash. micromamba (libsolv) on a
# narrow new env solves cleanly in ~1 min.
#
# Usage:
#     ssh bio
#     cd /data/biophys/carolinge/clawork/37_OXR/repo
#     bash project_RRM/ver2/pipeline/server_setup.sh

set -euo pipefail

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
SOURCE_ENV="allatom"
OPENMM_ENV="allatom_v2"
AMBER_ENV="amber24"
AMBERTOOLS_SPEC='ambertools=24.8=*nompi*py311*'

# 1. CUDA module (matches cudatoolkit pinned in allatom)
echo ">>> module load cuda/12.5"
module purge
module load cuda/12.5

# 2. Persistent project sub-dirs
echo ">>> Ensuring $PROJECT/ subdirs"
mkdir -p "$PROJECT"/{replicas_v2,force_fields/lib_amber,structures/{original,prepared,modified},pipeline_logs}

# 3. conda baseline (used for allatom_v2 clone)
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not in PATH"; exit 1
fi
CONDA_BASE="$(conda info --base)"
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

# 4. allatom_v2 = clone of ver_1's allatom (instant; no solve)
if conda env list | awk '{print $1}' | grep -qx "$OPENMM_ENV"; then
    echo ">>> $OPENMM_ENV already exists — skip clone"
else
    if ! conda env list | awk '{print $1}' | grep -qx "$SOURCE_ENV"; then
        echo "ERROR: source env '$SOURCE_ENV' not found. ver_1 not installed?"
        exit 1
    fi
    echo ">>> Cloning $SOURCE_ENV -> $OPENMM_ENV"
    conda create --name "$OPENMM_ENV" --clone "$SOURCE_ENV" -y
fi

# 5. micromamba (user-local, single static binary; required for amber24)
if [ ! -x "$HOME/bin/micromamba" ]; then
    echo ">>> Installing micromamba into ~/bin/"
    mkdir -p "$HOME/bin"
    cd /tmp
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xjf - bin/micromamba
    mv bin/micromamba "$HOME/bin/micromamba"
    chmod +x "$HOME/bin/micromamba"
    rm -rf bin
    cd -
fi
echo ">>> micromamba: $($HOME/bin/micromamba --version)"

# 6. amber24 env via micromamba (tleap + cpptraj + sander only)
export MAMBA_ROOT_PREFIX="$HOME/.conda"
if [ -d "$HOME/.conda/envs/$AMBER_ENV/bin" ] && \
   [ -x "$HOME/.conda/envs/$AMBER_ENV/bin/tleap" ]; then
    echo ">>> $AMBER_ENV already has tleap — skip create"
else
    echo ">>> Creating $AMBER_ENV via micromamba (~1-3 min, libsolv solver)"
    "$HOME/bin/micromamba" create -n "$AMBER_ENV" -c conda-forge -y \
        "$AMBERTOOLS_SPEC" python=3.11
fi

# 7. cpptraj from source (>= 6.26 required by modxna.sh; AmberTools 24.8
#    ships 6.24, which silently fails on `change crdset charge by`).
#    Build into ~/cpptraj_src and override amber24's stale 6.24 binary
#    with a wrapper that sets LD_LIBRARY_PATH.
CPPTRAJ_SRC="$HOME/cpptraj_src"
CPPTRAJ_BIN="$HOME/.conda/envs/$AMBER_ENV/bin/cpptraj"
if "$CPPTRAJ_BIN" --version 2>&1 | grep -qE 'V[6-9]\.(2[6-9]|[3-9][0-9]|[1-9][0-9]{2,})|V[7-9]'; then
    echo ">>> cpptraj already >= 6.26 — skip build"
else
    echo ">>> Building cpptraj from source (5-10 min) — amber24's 6.24 too old"
    if [ ! -d "$CPPTRAJ_SRC" ]; then
        git clone --depth 1 https://github.com/Amber-MD/cpptraj.git "$CPPTRAJ_SRC"
    fi
    (
        cd "$CPPTRAJ_SRC"
        export CONDA_PREFIX="$HOME/.conda/envs/$AMBER_ENV"
        export PATH="$CONDA_PREFIX/bin:$PATH"
        ./configure -shared --with-netcdf="$CONDA_PREFIX" --with-zlib="$CONDA_PREFIX" gnu
        make -j4
    )
    [ -x "$CPPTRAJ_SRC/bin/cpptraj" ] || { echo "ERROR: cpptraj build failed"; exit 1; }

    # Replace amber24/bin/cpptraj with a wrapper that sets LD_LIBRARY_PATH.
    if [ -f "$CPPTRAJ_BIN" ] && [ ! -f "${CPPTRAJ_BIN}.6.24.bak" ]; then
        mv "$CPPTRAJ_BIN" "${CPPTRAJ_BIN}.6.24.bak"
    fi
    cat > "$CPPTRAJ_BIN" <<'WRAPEOF'
#!/bin/bash
# Wrapper installed by server_setup.sh — points to source-built cpptraj
# and ensures amber24's libreadline/libnetcdf etc. are findable.
export LD_LIBRARY_PATH="$HOME/.conda/envs/amber24/lib:${LD_LIBRARY_PATH:-}"
exec "$HOME/cpptraj_src/bin/cpptraj" "$@"
WRAPEOF
    chmod +x "$CPPTRAJ_BIN"
fi

# 8. Verify
echo ""
echo "=== Verify amber24 binaries ==="
for cmd in cpptraj tleap sander; do
    if [ -x "$HOME/.conda/envs/$AMBER_ENV/bin/$cmd" ]; then
        echo "    OK: $HOME/.conda/envs/$AMBER_ENV/bin/$cmd"
    else
        echo "    FAIL: $cmd not in $AMBER_ENV"
        exit 1
    fi
done
echo "    cpptraj version: $($HOME/.conda/envs/$AMBER_ENV/bin/cpptraj --version 2>&1 | head -1)"

echo ""
echo "=== Verify allatom_v2 OpenMM ==="
conda activate "$OPENMM_ENV"
python -c "
from openmm import Platform
print('OpenMM platforms:', [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])
import openmm; print('OpenMM version:', openmm.__version__)
" || true

echo ""
echo "=== server_setup.sh DONE ==="
echo "Next: bash project_RRM/ver2/pipeline/build_modxna_residues.sh"
