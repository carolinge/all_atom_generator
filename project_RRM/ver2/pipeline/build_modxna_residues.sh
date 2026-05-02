#!/bin/bash
# build_modxna_residues.sh
# ========================
# Server-side: assemble 3 canonical-RNA modifications via modxna.sh.
#
#   8OG  = 8-oxoguanosine
#   PUU  = pseudouridine        (NOT modXNA's "PSU" — that's a designer ψ)
#   M1A  = N1-methyladenosine
#
# Each in 3 variants: internal | 5'-cap | 3'-cap.
#
# Output: $PROJECT/force_fields/lib_amber/<CODE>{I,3,5}.lib
#
# Prereq: server_setup.sh has installed amber24 env + repo synced.

set -euo pipefail

PROJECT="/data/biophys/carolinge/clawork/37_OXR"
REPO="$PROJECT/repo"
MODXNA_DIR="$REPO/project_RRM/ver2/force_fields/modxna"
MODXNA="$MODXNA_DIR/modxna.sh"
INPUTS_DIR="$REPO/project_RRM/ver2/pipeline/build_modxna_inputs"
OUTDIR="$PROJECT/force_fields/lib_amber"
AMBER_ENV="amber24"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Activate amber24 (where tleap/cpptraj/sander live)
if [ -d "$HOME/.conda/envs/$AMBER_ENV/bin" ]; then
    export PATH="$HOME/.conda/envs/$AMBER_ENV/bin:$PATH"
fi

if ! command -v cpptraj >/dev/null 2>&1; then
    echo "ERROR: cpptraj not in PATH. Did server_setup.sh complete?"
    exit 1
fi
if ! command -v tleap >/dev/null 2>&1; then
    echo "ERROR: tleap not in PATH"
    exit 1
fi
if [ ! -x "$MODXNA" ]; then
    echo "ERROR: $MODXNA not found. Did sync_to_server.sh run?"
    exit 1
fi

echo ">>> Using:"
echo "    cpptraj: $(which cpptraj)"
echo "    tleap:   $(which tleap)"
echo "    modxna:  $MODXNA"
echo ""

build_one () {
    local code=$1
    local variant=$2
    local resname infile extra logfile
    case "$variant" in
        internal) resname="${code}I"; infile="$INPUTS_DIR/$code.in"; extra="" ;;
        3cap)     resname="${code}3"; infile="$INPUTS_DIR/$code.in"; extra="--3cap" ;;
        5cap)
            resname="${code}5"
            infile="$OUTDIR/${code}_5cap.in"
            awk '{$1="5PO"; print}' "$INPUTS_DIR/$code.in" > "$infile"
            extra="--5cap"
            ;;
        *) echo "Bad variant $variant"; return 1 ;;
    esac

    [ -f "$infile" ] || { echo "  SKIP $code $variant (input $infile missing)"; return; }

    logfile="$OUTDIR/${resname}.modxna.log"
    echo ">>> Building $resname  (input: $(basename "$infile"))"

    if "$MODXNA" -i "$infile" -m "$resname" $extra --clean --nomin > "$logfile" 2>&1; then
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

for code in 8OG PUU M1A; do
    for variant in internal 3cap 5cap; do
        build_one "$code" "$variant" || true
    done
done

echo ""
echo "=== Assembled libs in $OUTDIR ==="
ls -la "$OUTDIR"/*.lib 2>/dev/null || echo "(none yet — see *.modxna.log)"
