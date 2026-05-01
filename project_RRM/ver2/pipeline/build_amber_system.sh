#!/bin/bash
# build_amber_system.sh
# ======================
# Server-side: take a prepared 4BS2-class PDB (with the modified residue
# already placed) and run tleap to produce a solvated AMBER prmtop+inpcrd.
#
# Usage:
#     bash build_amber_system.sh <name> <input_pdb>
#
# Example:
#     bash build_amber_system.sh 4BS2_8OG_G3 \
#         project_RRM/ver2/structures/modified/4BS2_8OG_G3.pdb
#
# Output (under $PROJECT/replicas_v2/<name>/system/):
#     <name>.prmtop
#     <name>.inpcrd
#     <name>.solvated.pdb
#     <name>.tleap.in
#     <name>.tleap.log
#
# Each replica directory under replicas_v2/<name>/r{1,2,3}/ then symlinks
# or copies prmtop+inpcrd from system/.

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <name> <input_pdb>"
    echo "  e.g.  $0 4BS2_8OG_G3 project_RRM/ver2/structures/modified/4BS2_8OG_G3.pdb"
    exit 1
fi

NAME=$1
INPUT_PDB=$2

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
LIB_DIR="$PROJECT/force_fields/lib_amber"
FRCMOD="$PROJECT/repo/project_RRM/ver2/force_fields/modxna/dat/frcmod.modxna"
SYSDIR="$PROJECT/replicas_v2/$NAME/system"

if [ ! -f "$INPUT_PDB" ]; then
    # Try resolving relative to repo
    if [ -f "$PROJECT/repo/$INPUT_PDB" ]; then
        INPUT_PDB="$PROJECT/repo/$INPUT_PDB"
    else
        echo "ERROR: input PDB not found: $INPUT_PDB"
        exit 1
    fi
fi
if [ ! -f "$FRCMOD" ]; then
    echo "ERROR: $FRCMOD not found. Sync repo first."
    exit 1
fi

if ! command -v tleap >/dev/null 2>&1; then
    echo "ERROR: tleap not in PATH. conda activate allatom_v2 first."
    exit 1
fi

mkdir -p "$SYSDIR"

# Auto-detect which modxna libs to load from $LIB_DIR
LIB_LOAD_LINES=""
for lib in "$LIB_DIR"/*.lib; do
    [ -f "$lib" ] || continue
    LIB_LOAD_LINES+="loadOff $lib"$'\n'
done
if [ -z "$LIB_LOAD_LINES" ]; then
    echo "WARN: no .lib files in $LIB_DIR. Did you run build_modxna_residues.sh?"
fi

cat > "$SYSDIR/$NAME.tleap.in" <<EOF
# Auto-generated tleap script for $NAME
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.water.tip3p

# modXNA shared parameters + assembled residue libraries
loadAmberParams $FRCMOD
$LIB_LOAD_LINES
mol = loadPdb $INPUT_PDB

# Solvation: 1.0 nm padding, 0.15 M NaCl
solvateBox mol TIP3PBOX 10.0
addIonsRand mol Na+ 0
addIonsRand mol Na+ 30
addIonsRand mol Cl- 30

saveAmberParm mol $SYSDIR/$NAME.prmtop $SYSDIR/$NAME.inpcrd
savePdb mol $SYSDIR/$NAME.solvated.pdb

quit
EOF

echo ">>> Running tleap"
tleap -f "$SYSDIR/$NAME.tleap.in" 2>&1 | tee "$SYSDIR/$NAME.tleap.log"

if [ ! -f "$SYSDIR/$NAME.prmtop" ] || [ ! -f "$SYSDIR/$NAME.inpcrd" ]; then
    echo "ERROR: prmtop or inpcrd not produced. Check $SYSDIR/$NAME.tleap.log"
    exit 1
fi

echo ""
echo "=== Built ==="
ls -la "$SYSDIR/$NAME".{prmtop,inpcrd,solvated.pdb}
echo ""
echo "Next: copy run_md.py + sub.sh into $PROJECT/replicas_v2/$NAME/r{1,2,3}/"
echo "      (use pipeline/spawn_replicas.sh)"
