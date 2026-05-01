#!/bin/bash
# build_modxna_residues.sh
# ========================
# Server-side: assemble the 3 canonical-RNA modifications we use
# into AMBER .lib files, via modXNA's modxna.sh.
#
#   8OG  = 8-oxoguanosine          (canonical, modXNA 8OG.mol2)
#   PUU  = pseudouridine           (canonical psi; NOT modXNA's "PSU")
#   M1A  = N1-methyladenosine      (canonical, modXNA M1A.mol2)
#
# Note: m6A is NOT in modXNA as canonical natural form. Use Bussi 2022
# externally if you need it.
#
# Each modification is built in three variants:
#   <CODE>I  internal   (RPO + RC3 + base)        — middle of chain
#   <CODE>3  3'-cap     (RPO + RC3 + base, --3cap)
#   <CODE>5  5'-cap     (5PO + RC3 + base, --5cap)
#
# Output: $PROJECT/force_fields/lib_amber/<CODE>{I,3,5}.lib
# These libs are SHARED across all replicas in replicas_v2/.
#
# Prerequisites:
#   - server_setup.sh has been run (allatom_v2 env + AmberTools)
#   - PROJECT root exists at /data/biophys/carolinge/clawork/37_OXR

set -euo pipefail

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
REPO="$PROJECT/repo"
MODXNA_DIR="$REPO/project_RRM/ver2/force_fields/modxna"
MODXNA="$MODXNA_DIR/modxna.sh"
INPUTS_DIR="$REPO/project_RRM/ver2/pipeline/build_modxna_inputs"
OUTDIR="$PROJECT/force_fields/lib_amber"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Ensure ambertools env active
if ! command -v cpptraj >/dev/null 2>&1; then
    echo "ERROR: cpptraj not found. conda activate allatom_v2 first."
    exit 1
fi
if [ ! -x "$MODXNA" ]; then
    echo "ERROR: $MODXNA not found or not executable."
    echo "Did you sync the repo to $REPO ?"
    exit 1
fi

build_one () {
    local code=$1     # 8OG, PUU, M1A
    local variant=$2  # internal | 5cap | 3cap

    local resname infile extra logfile
    case "$variant" in
        internal) resname="${code}I"; infile="$INPUTS_DIR/$code.in"; extra="" ;;
        3cap)     resname="${code}3"; infile="$INPUTS_DIR/$code.in"; extra="--3cap" ;;
        5cap)
            # 5'-cap requires a different backbone (5PO) — generate a
            # one-shot input file in OUTDIR.
            resname="${code}5"
            infile="$OUTDIR/${code}_5cap.in"
            awk '{$1="5PO"; print}' "$INPUTS_DIR/$code.in" > "$infile"
            extra="--5cap"
            ;;
        *) echo "Bad variant $variant"; return 1 ;;
    esac

    if [ ! -f "$infile" ]; then
        echo "  SKIP $code $variant (input $infile missing)"
        return
    fi

    logfile="$OUTDIR/${resname}.modxna.log"
    echo ">>> Building $resname  (input: $(basename "$infile"))"

    if "$MODXNA" -i "$infile" -m "$resname" $extra --clean --nomin > "$logfile" 2>&1; then
        # modxna.sh writes ${resname}.lib to CWD
        if [ -f "$OUTDIR/$resname.lib" ]; then
            echo "    OK -> $OUTDIR/$resname.lib"
        else
            echo "    WARN: $resname.lib not produced; check $logfile"
        fi
    else
        echo "    FAILED. Last 20 lines of $logfile:"
        tail -20 "$logfile" | sed 's/^/        /'
    fi
}

# Run all 3 mods × 3 variants = 9 builds
for code in 8OG PUU M1A; do
    for variant in internal 3cap 5cap; do
        build_one "$code" "$variant" || true
    done
done

echo ""
echo "=== Assembled libs in $OUTDIR ==="
ls -la "$OUTDIR"/*.lib 2>/dev/null || echo "(none yet — see *.modxna.log for errors)"
