#!/bin/bash
# spawn_replicas.sh
# =================
# Server-side: create r1/r2/r3 replica directories under
# replicas_v2/<name>/, each with run_md.py + sub.sh + replica.txt (seed),
# all sharing the prmtop+inpcrd in ../system/.
#
# Usage:
#     bash spawn_replicas.sh <name> [--n 3] [--prefix-seed 0]
#
# Example:
#     bash spawn_replicas.sh 4BS2_8OG_G3                # seeds 1,2,3
#     bash spawn_replicas.sh 4BS2_8OG_G3 --n 5          # seeds 1,2,3,4,5
#     bash spawn_replicas.sh 4BS2_PUU_U2 --prefix-seed 10  # seeds 11,12,13
#
# After this, submit:
#     for r in $PROJECT/replicas_v2/<name>/r*; do
#         (cd $r && dos2unix sub.sh run_md.py 2>/dev/null; sbatch sub.sh)
#     done

set -euo pipefail

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
REPO="$PROJECT/repo"
PIPELINE="$REPO/project_RRM/ver2/pipeline"

NAME=""
N=3
PREFIX_SEED=0
GPU="a100"   # default; can override with --gpu h200

while [ $# -gt 0 ]; do
    case "$1" in
        --n) N=$2; shift 2 ;;
        --prefix-seed) PREFIX_SEED=$2; shift 2 ;;
        --gpu) GPU=$2; shift 2 ;;
        --*) echo "Unknown option $1"; exit 1 ;;
        *)
            if [ -z "$NAME" ]; then NAME=$1; shift
            else echo "Extra arg: $1"; exit 1
            fi ;;
    esac
done

if [ -z "$NAME" ]; then
    echo "Usage: $0 <name> [--n 3] [--prefix-seed 0] [--gpu a100|h200]"
    exit 1
fi

case "$GPU" in
    a100)  TEMPLATE="$PIPELINE/sub.sh.template" ;;
    h200)  TEMPLATE="$PIPELINE/sub_h200.sh.template" ;;
    *)     echo "Unknown --gpu $GPU (use a100 or h200)"; exit 1 ;;
esac

SYSDIR="$PROJECT/replicas_v2/$NAME/system"
if [ ! -f "$SYSDIR"/*.prmtop ] 2>/dev/null; then
    if ! ls "$SYSDIR"/*.prmtop >/dev/null 2>&1; then
        echo "ERROR: no prmtop in $SYSDIR — run build_amber_system.sh first."
        exit 1
    fi
fi

for i in $(seq 1 $N); do
    SEED=$((PREFIX_SEED + i))
    REP_DIR="$PROJECT/replicas_v2/$NAME/r$i"
    mkdir -p "$REP_DIR"
    cp "$PIPELINE/run_md.py" "$REP_DIR/run_md.py"
    sed "s/__JOBNAME__/${NAME}_r${i}/" "$TEMPLATE" \
        > "$REP_DIR/sub.sh"
    chmod +x "$REP_DIR/sub.sh"
    echo "$SEED" > "$REP_DIR/replica.txt"
    # Defensive: strip CRLF on copied scripts
    sed -i 's/\r$//' "$REP_DIR/run_md.py" "$REP_DIR/sub.sh"
    echo "    spawned $REP_DIR  (seed=$SEED)"
done

echo ""
echo "=== Spawned $N replicas under $PROJECT/replicas_v2/$NAME/ ==="
echo "Submit them with:"
echo "  for r in $PROJECT/replicas_v2/$NAME/r*; do (cd \$r && sbatch sub.sh); done"
