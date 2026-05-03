#!/bin/bash
# build_free_rna.sh
# =================
# Build a small free-RNA validation system de novo using tleap's
# sequence command. No input PDB needed — tleap places residues in
# idealized linear geometry, then solvates.
#
# Usage:
#     bash build_free_rna.sh <name> "<sequence>"
# Example (m1A in the middle of a tetranucleotide):
#     bash build_free_rna.sh m1A_test "A5 A M1A A3"
#
# Convention: sequence is space-separated AMBER residue names. Use
# 5'/3' caps for terminal residues:
#     standard:    A5 ... A3   C5 ... C3   G5 ... G3   U5 ... U3
#     modxna:      M1A5 (5cap) M1A (internal) M1A3 (3cap)
#                  PUU5  PUU  PUU3
#                  8OG5  8OG  8OG3
#
# Output: $PROJECT/replicas_v2/<name>/system/{prmtop,inpcrd,solvated.pdb}

set -euo pipefail

if [ $# -ne 2 ]; then
    echo 'Usage: $0 <name> "<sequence>"'
    echo 'Example: $0 m1A_test "A5 A M1A A3"'
    exit 1
fi

NAME=$1
SEQ=$2

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
REPO="$PROJECT/repo"
LIB_DIR="$PROJECT/force_fields/lib_amber"
FRCMOD="$REPO/project_RRM/ver2/force_fields/modxna/dat/frcmod.modxna"
SYSDIR="$PROJECT/replicas_v2/$NAME/system"
AMBER_ENV="amber24"

if [ -d "$HOME/.conda/envs/$AMBER_ENV/bin" ]; then
    export PATH="$HOME/.conda/envs/$AMBER_ENV/bin:$PATH"
fi
command -v tleap >/dev/null 2>&1 || { echo "ERROR: tleap not in PATH"; exit 1; }

mkdir -p "$SYSDIR"

LIB_LOAD_LINES=""
for lib in "$LIB_DIR"/*.lib; do
    [ -f "$lib" ] || continue
    LIB_LOAD_LINES+="loadOff $lib"$'\n'
done

cat > "$SYSDIR/$NAME.tleap.in" <<EOF
# Auto-generated tleap script for free-RNA system "$NAME"
# Sequence: $SEQ
source leaprc.RNA.OL3
source leaprc.water.tip3p

loadAmberParams $FRCMOD
$LIB_LOAD_LINES
mol = sequence { $SEQ }

# Save unsolvated topology (small) for inspection
saveAmberParm mol $SYSDIR/$NAME.dry.prmtop $SYSDIR/$NAME.dry.inpcrd
savePdb mol $SYSDIR/$NAME.dry.pdb

# Solvate with TIP3P, 1.0 nm padding, neutralize
solvateBox mol TIP3PBOX 10.0
addIonsRand mol Na+ 0
addIonsRand mol Cl- 0

saveAmberParm mol $SYSDIR/$NAME.prmtop $SYSDIR/$NAME.inpcrd
savePdb mol $SYSDIR/$NAME.solvated.pdb
quit
EOF

echo ">>> tleap"
tleap -f "$SYSDIR/$NAME.tleap.in" 2>&1 | tee "$SYSDIR/$NAME.tleap.log"

[ -f "$SYSDIR/$NAME.prmtop" ] && [ -f "$SYSDIR/$NAME.inpcrd" ] || {
    echo "ERROR: prmtop/inpcrd not produced"
    exit 1
}

echo ""
echo "=== Built ==="
ls -la "$SYSDIR/$NAME".{prmtop,inpcrd,solvated.pdb}
echo ""
echo "Next: bash project_RRM/ver2/pipeline/spawn_replicas.sh $NAME [--gpu a100|h200]"
