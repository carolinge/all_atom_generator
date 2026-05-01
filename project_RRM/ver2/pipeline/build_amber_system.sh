#!/bin/bash
# build_amber_system.sh
# ======================
# Server-side: take a prepared 4BS2 PDB (with the modified residue already
# placed via NeRF) and run tleap to generate the AMBER prmtop + inpcrd.
#
# OpenMM reads the prmtop directly via AmberPrmtopFile — no OpenMM XML
# conversion needed.
#
# Usage:
#     bash build_amber_system.sh <input_pdb> <output_prefix>
# Example:
#     bash build_amber_system.sh \
#         project_RRM/ver2/structures/modified/4BS2_8OG_G3.pdb \
#         runs/4BS2_8OG_G3
#
# Prerequisites:
#   - allatom_v2 conda env with ambertools active
#   - .lib files for the modified residues already built (run
#     build_modxna_residues.sh first)

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_pdb> <output_prefix>"
    exit 1
fi

INPUT_PDB=$1
OUT_PREFIX=$2

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
LIB_DIR="$REPO_ROOT/project_RRM/ver2/force_fields/openmm_xml/lib_amber"
FRCMOD="$REPO_ROOT/project_RRM/ver2/force_fields/modxna/dat/frcmod.modxna"

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: input PDB $INPUT_PDB not found"
    exit 1
fi
if [ ! -f "$FRCMOD" ]; then
    echo "ERROR: $FRCMOD not found"
    exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")"

# Auto-detect which modxna libs to load by scanning lib_amber/
LIB_LOAD_LINES=""
for lib in "$LIB_DIR"/*.lib; do
    [ -f "$lib" ] || continue
    LIB_LOAD_LINES+="loadOff $lib"$'\n'
done

cat > "${OUT_PREFIX}.tleap.in" <<EOF
# Auto-generated tleap script for $INPUT_PDB
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.water.tip3p

# Load modXNA shared parameters (frcmod) and assembled residue libraries.
loadAmberParams $FRCMOD
$LIB_LOAD_LINES
# Read the prepared, modified structure
mol = loadPdb $INPUT_PDB

# Solvation: 1.0 nm padding, 0.15 M NaCl
solvateBox mol TIP3PBOX 10.0
addIonsRand mol Na+ 0
addIonsRand mol Na+ 30
addIonsRand mol Cl- 30

# Output
saveAmberParm mol ${OUT_PREFIX}.prmtop ${OUT_PREFIX}.inpcrd
savePdb mol ${OUT_PREFIX}.solvated.pdb

quit
EOF

echo ">>> Running tleap"
tleap -f "${OUT_PREFIX}.tleap.in" 2>&1 | tee "${OUT_PREFIX}.tleap.log"

# Validate
if [ ! -f "${OUT_PREFIX}.prmtop" ] || [ ! -f "${OUT_PREFIX}.inpcrd" ]; then
    echo "ERROR: prmtop or inpcrd not generated. Check ${OUT_PREFIX}.tleap.log"
    exit 1
fi

echo ""
echo "=== Built ==="
ls -la "${OUT_PREFIX}.prmtop" "${OUT_PREFIX}.inpcrd" "${OUT_PREFIX}.solvated.pdb"
