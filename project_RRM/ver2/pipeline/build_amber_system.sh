#!/bin/bash
# build_amber_system.sh
# ======================
# Server-side: take a prepared PDB (with the modified residue placed)
# and run tleap to produce a solvated AMBER prmtop+inpcrd.
#
# Usage:
#     bash build_amber_system.sh <name> <input_pdb>
# Example:
#     bash build_amber_system.sh 4BS2_8OG_G3 \
#         project_RRM/ver2/structures/modified/4BS2_8OG_G3.pdb
#
# Output: $PROJECT/replicas_v2/<name>/system/{prmtop,inpcrd,solvated.pdb}
# tleap runs in the amber24 env; OpenMM later reads the prmtop in
# allatom_v2 via AmberPrmtopFile.

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <name> <input_pdb>"
    exit 1
fi

NAME=$1
INPUT_PDB=$2

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
REPO="$PROJECT/repo"
LIB_DIR="$PROJECT/force_fields/lib_amber"
FRCMOD="$REPO/project_RRM/ver2/force_fields/modxna/dat/frcmod.modxna"
SYSDIR="$PROJECT/replicas_v2/$NAME/system"
AMBER_ENV="amber24"

# Resolve input PDB
if [ ! -f "$INPUT_PDB" ]; then
    if [ -f "$REPO/$INPUT_PDB" ]; then
        INPUT_PDB="$REPO/$INPUT_PDB"
    else
        echo "ERROR: input PDB not found: $INPUT_PDB"
        exit 1
    fi
fi
[ -f "$FRCMOD" ] || { echo "ERROR: $FRCMOD not found"; exit 1; }

# Activate amber24 for tleap
if [ -d "$HOME/.conda/envs/$AMBER_ENV/bin" ]; then
    export PATH="$HOME/.conda/envs/$AMBER_ENV/bin:$PATH"
fi
command -v tleap >/dev/null 2>&1 || { echo "ERROR: tleap not in PATH"; exit 1; }

mkdir -p "$SYSDIR"

# Auto-detect modxna libs to load
LIB_LOAD_LINES=""
for lib in "$LIB_DIR"/*.lib; do
    [ -f "$lib" ] || continue
    LIB_LOAD_LINES+="loadOff $lib"$'\n'
done
[ -n "$LIB_LOAD_LINES" ] || \
    echo "WARN: no .lib in $LIB_DIR — did build_modxna_residues.sh run?"

cat > "$SYSDIR/$NAME.tleap.in" <<EOF
# Auto-generated tleap script for $NAME
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.water.tip3p

loadAmberParams $FRCMOD
$LIB_LOAD_LINES
mol = loadPdb $INPUT_PDB

solvateBox mol TIP3PBOX 10.0
addIonsRand mol Na+ 0
addIonsRand mol Na+ 30
addIonsRand mol Cl- 30

saveAmberParm mol $SYSDIR/$NAME.prmtop $SYSDIR/$NAME.inpcrd
savePdb mol $SYSDIR/$NAME.solvated.pdb

quit
EOF

echo ">>> tleap -f $SYSDIR/$NAME.tleap.in"
tleap -f "$SYSDIR/$NAME.tleap.in" 2>&1 | tee "$SYSDIR/$NAME.tleap.log"

[ -f "$SYSDIR/$NAME.prmtop" ] && [ -f "$SYSDIR/$NAME.inpcrd" ] || {
    echo "ERROR: prmtop/inpcrd not produced. Check $SYSDIR/$NAME.tleap.log"
    exit 1
}

echo ""
echo "=== Built ==="
ls -la "$SYSDIR/$NAME".{prmtop,inpcrd,solvated.pdb}
echo ""
echo "Next: bash project_RRM/ver2/pipeline/spawn_replicas.sh $NAME"
