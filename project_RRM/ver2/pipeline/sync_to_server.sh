#!/bin/bash
# sync_to_server.sh
# =================
# Push local d:/all_atom (project_RRM/ subset) to the bio cluster.
#
# Uses `tar | ssh tar` instead of rsync, because Git Bash on Windows
# doesn't ship rsync by default. tar + ssh is single-pass, fast, and
# works with stock Git for Windows + OpenSSH.
#
# Run from EITHER:
#   - Git Bash:  bash project_RRM/ver2/pipeline/sync_to_server.sh
#   - Or via Claude's Bash tool (same environment)
#
# DO NOT run from `bash` in cmd.exe on Windows 10+ — that defaults to
# WSL bash which has no Linux distro on this system.
#
# Destination: bio:/data/biophys/carolinge/clawork/37_OXR/repo/

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
DEST_HOST="bio"
DEST_ROOT="/data/biophys/carolinge/clawork/37_OXR/repo"

cd "$REPO_ROOT"

# 1. Verify ssh works
echo ">>> Testing ssh to $DEST_HOST"
if ! ssh -o BatchMode=yes -o ConnectTimeout=10 "$DEST_HOST" "true"; then
    echo "ERROR: ssh to $DEST_HOST failed. Check ~/.ssh/config and key."
    exit 1
fi

# 2. Normalize CRLF -> LF on shell + python + tleap scripts
# IMPORTANT: also covers project_RRM/ver2/force_fields/modxna/ — git on
# Windows can convert these on checkout, and modxna.sh's shebang then
# breaks on Linux ("/bin/bash^M: bad interpreter").
echo ">>> Normalizing CRLF -> LF (project_RRM + modxna scripts)"
find project_RRM/ \( -name '*.sh' -o -name '*.py' -o -name '*.template' \
                   -o -name '*.in' -o -name '*.tleap' \) -print0 \
    | while IFS= read -r -d '' f; do
        if grep -lq $'\r' "$f" 2>/dev/null; then
            sed -i 's/\r$//' "$f"
            echo "    fixed: $f"
        fi
    done

# 3. Build the tar exclude list
EXCLUDES=(
    --exclude='.git'
    --exclude='.last_push'
    --exclude='__pycache__'
    --exclude='*.pyc'
    --exclude='.vscode'
    --exclude='.claude'
    --exclude='*.aux'
    --exclude='*.log'
    --exclude='*.toc'
    --exclude='project_RRM/ver1/output'
    --exclude='project_RRM/ver1/uploaded'
    --exclude='project_RRM/ver1/prepared'
    --exclude='project_RRM/ver1/example_2RRM'
    --exclude='project_RRM/ver2/runs'
    --exclude='project_RRM/ver2/results'
)

# 4. Tar over SSH
echo ""
echo ">>> Streaming repo to $DEST_HOST:$DEST_ROOT/"
echo "    (excluding ver_1 output/uploaded/prepared and ver_2 runs/results)"

# Pre-create dest dir
ssh "$DEST_HOST" "mkdir -p $DEST_ROOT"

# Stream tar
tar -cf - "${EXCLUDES[@]}" \
    app/ \
    auto_push.sh \
    docs/ \
    DEVELOPER.md \
    UPDATE_LOG.md \
    environment.yml \
    environment_server.yml \
    setup_server.sh \
    launch_ui.bat \
    project_RRM/ \
    .gitignore \
    | ssh "$DEST_HOST" "tar -xf - -C $DEST_ROOT"

echo ""
echo ">>> Verify on remote"
ssh "$DEST_HOST" "ls -la $DEST_ROOT/project_RRM/ver2/pipeline/ | head -15"

echo ""
echo "=== sync_to_server.sh DONE ==="
echo "Next: ssh $DEST_HOST"
echo "      cd $DEST_ROOT"
echo "      bash project_RRM/ver2/pipeline/server_setup.sh"
