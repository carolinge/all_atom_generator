#!/bin/bash
# run_analysis.sh
# ===============
# Run analyze_v2.py on all 6 4BS2 trajectories (WT and 8OG, 3 replicas
# each), then compare_v2.py to aggregate + compare + plot. Run on the
# bio cluster from the repo root.

set -euo pipefail

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
RESULTS="$PROJECT/results"
ENV_BIN="$HOME/.conda/envs/allatom_v2/bin"
ANALYZE="$(dirname "$0")/analyze_v2.py"
COMPARE="$(dirname "$0")/compare_v2.py"

# Resolve modified residue id from solvated PDB (look up the 8OG / OGI residue).
# Use a one-shot python on the 8OG_G3 prmtop to find its resid in the solvated
# system (protein offset + RNA position 3).
MOD_RESID=$($ENV_BIN/python -c "
import MDAnalysis as mda
u = mda.Universe('$PROJECT/replicas_v2/4BS2_8OG_G3/system/4BS2_8OG_G3.solvated.pdb')
mods = u.select_atoms('resname OGI or resname 8OG').residues
if not len(mods):
    raise SystemExit('no OGI/8OG residue found')
print(mods[0].resid)
")
echo ">>> Modified residue resid in solvated system: $MOD_RESID"

mkdir -p "$RESULTS"

# Per-replica analysis
for sys_label in 4BS2_WT 4BS2_8OG_G3; do
    TOPO="$PROJECT/replicas_v2/$sys_label/system/${sys_label}.solvated.pdb"
    if [ "$sys_label" = "4BS2_8OG_G3" ]; then
        MOD_NAME="8OG"
    else
        MOD_NAME="G"
    fi
    PRMTOP="$PROJECT/replicas_v2/$sys_label/system/${sys_label}.prmtop"
    for r in $PROJECT/replicas_v2/$sys_label/r*; do
        rep=$(basename "$r")
        traj="$r/traj.dcd"
        if [ ! -f "$traj" ]; then
            echo "  SKIP $sys_label/$rep (traj.dcd missing)"
            continue
        fi
        out_check="$RESULTS/$sys_label/$rep/summary.json"
        if [ -f "$out_check" ]; then
            echo "  exists $sys_label/$rep — skip (rerun: rm $out_check)"
            continue
        fi
        echo ""
        echo ">>> Analyzing $sys_label / $rep"
        $ENV_BIN/python "$ANALYZE" \
            --topology "$TOPO" --prmtop "$PRMTOP" --trajectory "$traj" \
            --label "$sys_label" --replica "$rep" \
            --mod-resid "$MOD_RESID" --mod-resname "$MOD_NAME" \
            --out "$RESULTS"
    done
done

# Aggregate + compare
echo ""
echo ">>> compare_v2.py — aggregate WT vs 8OG"
$ENV_BIN/python "$COMPARE" \
    --label-replicas WT  4BS2_WT/r1 4BS2_WT/r2 4BS2_WT/r3 \
    --label-replicas 8OG 4BS2_8OG_G3/r1 4BS2_8OG_G3/r2 4BS2_8OG_G3/r3 \
    --mod-resid "$MOD_RESID" \
    --in "$RESULTS" \
    --out "$RESULTS/_aggregate"

echo ""
echo "=== DONE ==="
echo "Per-replica:    $RESULTS/4BS2_*/r*/summary.json"
echo "Aggregate:      $RESULTS/_aggregate/"
echo "Report:         $RESULTS/_aggregate/report.md"
