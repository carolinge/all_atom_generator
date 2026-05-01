#!/bin/bash
# build_modxna_residues.sh
# =========================
# Run modxna.sh on a Linux host with AmberTools installed, producing
# assembled .lib files for our 3 modXNA-covered modifications:
#
#   8OG  = 8-oxoguanosine          (canonical)
#   PUU  = pseudouridine           (canonical, NOT modXNA's "PSU" code)
#   M1A  = N1-methyladenosine      (canonical)
#
# Note: m6A is NOT in modXNA as canonical. Use Bussi 2022 separately.
#
# Each mod is built in three variants:
#   internal  -- standard middle-of-chain residue (uses RPO backbone)
#   5cap      -- 5'-terminal (uses 5PO backbone + --5cap flag)
#   3cap      -- 3'-terminal (uses RPO backbone + --3cap flag)
#
# Output: project_RRM/ver2/force_fields/openmm_xml/lib_amber/<CODE>.lib (etc.)
# Bring those back to the Windows side and feed them to convert_lib_to_openmm.py.
#
# Prerequisites on the server (run setup_modxna_server.sh first):
#   - conda env with ambertools (cpptraj, tleap, sander)
#   - modxna.sh in PATH or invoked by relative path

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
MODXNA="$REPO_ROOT/project_RRM/ver2/force_fields/modxna/modxna.sh"
INPUTS="$REPO_ROOT/project_RRM/ver2/pipeline/build_modxna_inputs"
OUTDIR="$REPO_ROOT/project_RRM/ver2/force_fields/openmm_xml/lib_amber"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

if ! command -v cpptraj >/dev/null 2>&1; then
    echo "ERROR: cpptraj not found in PATH. Activate the AmberTools conda env first."
    exit 1
fi
if ! command -v tleap >/dev/null 2>&1; then
    echo "ERROR: tleap not found. Same env."
    exit 1
fi
if [ ! -x "$MODXNA" ]; then
    echo "ERROR: $MODXNA not executable."
    exit 1
fi

build_one () {
    local code=$1     # e.g. 8OG
    local variant=$2  # internal | 5cap | 3cap
    local resname    # name of the assembled residue
    local infile="$INPUTS/$code.in"
    local logfile="$OUTDIR/${code}_${variant}.modxna.log"

    if [ ! -f "$infile" ]; then
        echo "  SKIP $code $variant (input file $infile missing)"
        return
    fi

    # Construct per-variant input + flags
    local extra=""
    case "$variant" in
        internal) resname="${code}I"; extra="" ;;
        5cap)     resname="${code}5"; extra="--5cap"
                  # For 5cap, override the backbone to 5PO
                  infile="$OUTDIR/${code}_5cap.in"
                  awk '{$1="5PO"; print}' "$INPUTS/$code.in" > "$infile" ;;
        3cap)     resname="${code}3"; extra="--3cap" ;;
    esac

    echo ">>> Building $code [$variant] -> ${resname}.lib"
    "$MODXNA" -i "$infile" -m "$resname" $extra --clean --nomin > "$logfile" 2>&1 \
        || { echo "  FAILED. Log: $logfile"; return 1; }
    echo "  OK: $(ls -la "$resname.lib" 2>/dev/null | awk '{print $5,$9}')"
}

for code in 8OG PUU M1A; do
    build_one "$code" internal || true
    build_one "$code" 3cap     || true
    build_one "$code" 5cap     || true
done

echo ""
echo "=== Assembled libs in $OUTDIR ==="
ls -la "$OUTDIR"/*.lib 2>/dev/null || echo "(none)"
