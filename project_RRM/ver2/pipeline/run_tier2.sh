#!/bin/bash
# Run analyze_tier2.py on the 3 8OG replicas (Tier 2 only relevant for
# the modified system — WT G has no N7-H7 / O8 to analyze, but we run
# the WT residue too so the backbone dihedral comparison has both sides).
set -euo pipefail
PROJECT="/data/biophys/carolinge/clawork/37_OXR"
RESULTS="$PROJECT/results"
ENV_BIN="$HOME/.conda/envs/allatom_v2/bin"
SCRIPT="$(dirname "$0")/analyze_tier2.py"

# Resolve mod resid (177)
MOD_RESID=$($ENV_BIN/python -c "
import MDAnalysis as mda
u = mda.Universe('$PROJECT/replicas_v2/4BS2_8OG_G3/system/4BS2_8OG_G3.solvated.pdb')
print(u.select_atoms('resname OGI or resname 8OG').residues[0].resid)
")
echo ">>> Modified residue resid = $MOD_RESID"

for sys_label in 4BS2_WT 4BS2_8OG_G3; do
    TOPO="$PROJECT/replicas_v2/$sys_label/system/${sys_label}.solvated.pdb"
    PRMTOP="$PROJECT/replicas_v2/$sys_label/system/${sys_label}.prmtop"
    if [ "$sys_label" = "4BS2_8OG_G3" ]; then MOD_NAME="8OG"; else MOD_NAME="G"; fi
    for r in $PROJECT/replicas_v2/$sys_label/r*; do
        rep=$(basename "$r")
        traj="$r/traj.dcd"
        [ -f "$traj" ] || { echo "skip $sys_label/$rep (no traj)"; continue; }
        if [ -f "$RESULTS/$sys_label/$rep/tier2_summary.json" ]; then
            echo "skip $sys_label/$rep (cached)"; continue
        fi
        echo ""
        echo ">>> Tier-2 $sys_label / $rep"
        $ENV_BIN/python "$SCRIPT" \
            --topology "$TOPO" --prmtop "$PRMTOP" --trajectory "$traj" \
            --label "$sys_label" --replica "$rep" \
            --mod-resid "$MOD_RESID" --mod-resname "$MOD_NAME" \
            --out "$RESULTS"
    done
done
echo ""
echo "=== Tier-2 DONE ==="
echo "Per-replica: $RESULTS/4BS2_*/r*/tier2_*.csv"
