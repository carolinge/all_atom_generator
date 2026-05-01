#!/bin/bash
# sync_to_server.sh
# =================
# Push local d:/all_atom (project_RRM/) to the bio cluster via rsync over
# the `bio` SSH alias. Handles CRLF -> LF conversion for shell + python
# scripts (Newton's sbatch chokes on DOS line endings — see manual).
#
# Usage (from local repo root):
#     bash project_RRM/ver2/pipeline/sync_to_server.sh
#
# Excludes:
#   - .git/, .last_push, output/ (large), runs/ (large)
#   - Node/Python caches
#   - LaTeX build artifacts
#
# Destination: bio:/data/biophys/carolinge/clawork/37_OXR/repo/

set -euo pipefail

# Resolve repo root robustly
REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
DEST_HOST="bio"
DEST_ROOT="/data/biophys/carolinge/clawork/37_OXR/repo"

cd "$REPO_ROOT"

# 1. Convert CRLF -> LF on shell + python scripts (in place, idempotent)
echo ">>> Normalizing line endings (CRLF -> LF) for *.sh and *.py under project_RRM/ver2/"
find project_RRM/ver2/ \( -name '*.sh' -o -name '*.py' \) -print0 \
    | while IFS= read -r -d '' f; do
        # Skip if no CR present (avoid touching mtime)
        if grep -lq $'\r' "$f" 2>/dev/null; then
            sed -i 's/\r$//' "$f"
            echo "    fixed: $f"
        fi
    done

# 2. rsync (note: trailing slash on source = "contents of"; on dest = "into")
echo ""
echo ">>> rsync -> $DEST_HOST:$DEST_ROOT/"
rsync -avz --delete-after \
    --exclude='.git/' \
    --exclude='.last_push' \
    --exclude='__pycache__/' \
    --exclude='*.pyc' \
    --exclude='.vscode/' \
    --exclude='.claude/' \
    --exclude='*.aux' \
    --exclude='*.log' \
    --exclude='*.out' \
    --exclude='*.toc' \
    --exclude='project_RRM/ver1/output/' \
    --exclude='project_RRM/ver1/uploaded/' \
    --exclude='project_RRM/ver1/prepared/' \
    --exclude='project_RRM/ver1/example_2RRM/' \
    --exclude='project_RRM/ver2/runs/' \
    --exclude='project_RRM/ver2/results/' \
    "$REPO_ROOT/" \
    "$DEST_HOST:$DEST_ROOT/"

echo ""
echo "=== sync_to_server.sh DONE ==="
echo "Next steps on the cluster:"
echo "  ssh $DEST_HOST"
echo "  cd $DEST_ROOT"
echo "  bash project_RRM/ver2/pipeline/server_setup.sh         # one-time"
echo "  bash project_RRM/ver2/pipeline/build_modxna_residues.sh # one-time per modXNA update"
