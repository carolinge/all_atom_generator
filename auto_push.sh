#!/bin/bash
# auto_push.sh — commit and push if >4 hours since last push
# Called by Claude at session start and optionally via Task Scheduler

REPO="d:/all_atom"
STAMP="$REPO/.last_push"

NOW=$(date +%s)
LAST=0
[ -f "$STAMP" ] && LAST=$(cat "$STAMP")
DIFF=$((NOW - LAST))

if [ "$DIFF" -lt 14400 ]; then
  echo "[auto_push] Last push was $((DIFF/60)) min ago, skipping."
  exit 0
fi

cd "$REPO" || exit 1

# Stage tracked files that may have changed
git add app/ DEVELOPER.md UPDATE_LOG.md environment.yml launch_ui.bat .gitignore

# Commit only if there are staged changes
if ! git diff --cached --quiet; then
  MSG="Auto-update: $(date '+%Y-%m-%d %H:%M')"
  git commit -m "$MSG"
  echo "[auto_push] Committed: $MSG"
fi

# Push if there are unpushed commits
UNPUSHED=$(git log origin/main..HEAD --oneline 2>/dev/null)
if [ -n "$UNPUSHED" ]; then
  git push origin main && echo "$NOW" > "$STAMP" && echo "[auto_push] Pushed OK."
else
  echo "[auto_push] Nothing to push."
  echo "$NOW" > "$STAMP"
fi
