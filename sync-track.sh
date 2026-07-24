#!/usr/bin/env bash
# sync-track.sh <track-name> "<commit message>"
#
# Run from inside a track worktree (~/projects/riemann-tracks/<name>) after every real finding
# (every LEDGER append). Commits pending work, merges it into main LOCALLY in the canonical
# ~/projects/riemann checkout (same .git object store, no network needed — fast, proven
# reliable), rebases main back into the track branch so it sees what other tracks have landed,
# and best-effort pushes main to origin (GitHub) with a speed-timeout guard so a flaky network
# can't hang the agent for minutes.
#
# IMPORTANT PRECONDITION: this assumes ~/projects/riemann's working directory is NOT being
# actively edited by another live session right now. If a second live session is ever working
# directly in that directory again, this script's local-merge step is unsafe (it can collide
# with uncommitted edits there) — switch back to an origin-only fetch/rebase/push design in
# that case instead (see git history of this file for that version).
#
# On any real conflict it aborts and does NOT guess a resolution — it logs to SESSIONS.md
# NEEDS JACOB and exits nonzero. See CLAUDE.md "Multi-session protocol".
set -uo pipefail

NAME="${1:?usage: sync-track.sh <track-name> \"<commit message>\"}"
MSG="${2:-sync}"
CANON="$HOME/projects/riemann"
TRACK_DIR="$HOME/projects/riemann-tracks/$NAME"
BRANCH="track/$NAME"
TS="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
SESSIONS="$CANON/SESSIONS.md"

if [ "$(pwd -P)" != "$(cd "$TRACK_DIR" && pwd -P)" ]; then
  echo "sync-track.sh: run this from $TRACK_DIR (current dir: $(pwd))" >&2
  exit 2
fi

log_sync() {
  printf -- '- %s track:%s %s %s\n' "$TS" "$NAME" "$1" "$2" >> "$SESSIONS"
}

flag_conflict() {
  {
    echo ""
    echo "### CONFLICT $TS — track:$NAME"
    echo "$1"
    echo ""
  } >> "$SESSIONS"
  echo "sync-track.sh: CONFLICT — logged to SESSIONS.md NEEDS JACOB, stopping. Do not guess a resolution." >&2
}

# 1. commit any pending work on the track branch
if [ -n "$(git status --porcelain)" ]; then
  git add -A
  git commit -m "[track:$NAME] $MSG" >/dev/null
fi
COMMIT_SHA="$(git rev-parse --short HEAD)"

# 2. merge track branch into main, in the canonical checkout — local, no network, near-instant
if ! git -C "$CANON" merge --no-ff "$BRANCH" -m "merge $BRANCH ($TS)" >/tmp/sync-track-merge-$NAME.log 2>&1; then
  git -C "$CANON" merge --abort >/dev/null 2>&1
  flag_conflict "merge $BRANCH -> main failed. Your work IS still committed locally at $COMMIT_SHA. Log: /tmp/sync-track-merge-$NAME.log"
  log_sync "FAILED" "commit:$COMMIT_SHA merge-conflict-see-tmp-log"
  exit 1
fi

# 3. best-effort push to origin — speed-timeout guard so a flaky network can't hang this for
#    minutes; non-fatal either way, this is a backup, not part of the coordination mechanism
if git -C "$CANON" -c http.lowSpeedLimit=1000 -c http.lowSpeedTime=15 push origin main \
      >/tmp/sync-track-push-$NAME.log 2>&1; then
  PUSH_STATE="pushed"
else
  PUSH_STATE="push-failed(non-fatal)"
fi

# 4. rebase main back into the track branch so this track sees everyone else's landed work
if ! git rebase main >/tmp/sync-track-rebase-$NAME.log 2>&1; then
  git rebase --abort >/dev/null 2>&1
  flag_conflict "rebase main -> track/$NAME failed after a successful merge to main. Your merge DID land on main; only your local branch's rebase failed. Log: /tmp/sync-track-rebase-$NAME.log"
  log_sync "PARTIAL" "commit:$COMMIT_SHA merged-to-main rebase-conflict-see-tmp-log $PUSH_STATE"
  exit 1
fi

log_sync "ok" "commit:$COMMIT_SHA $PUSH_STATE"
echo "sync-track.sh: OK — $NAME synced ($COMMIT_SHA, $PUSH_STATE)"
