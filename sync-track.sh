#!/usr/bin/env bash
# sync-track.sh <track-name> "<commit message>"
#
# Run from inside a track worktree (~/projects/riemann-tracks/<name>) after every real finding
# (every LEDGER append). Commits pending work, then syncs PURELY THROUGH origin (GitHub) —
# fetch origin/main, rebase onto it, push. This deliberately never touches the canonical
# ~/projects/riemann checkout's working directory: that directory may be actively edited,
# uncommitted, by another live session at any time, and merging into a dirty working tree is
# unsafe. Going through origin means this track and any other session only ever interact via
# normal git push/pull, never by touching each other's working tree.
#
# On any real conflict (rebase hits actually-overlapping content) it aborts and does NOT guess
# a resolution — it logs to SESSIONS.md NEEDS JACOB (in THIS worktree's own copy, which travels
# with the next commit) and exits nonzero. See CLAUDE.md "Multi-session protocol".
set -uo pipefail

NAME="${1:?usage: sync-track.sh <track-name> \"<commit message>\"}"
MSG="${2:-sync}"
TRACK_DIR="$HOME/projects/riemann-tracks/$NAME"
TS="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
SESSIONS="./SESSIONS.md"
MAX_RETRIES=3

if [ "$(pwd -P)" != "$(cd "$TRACK_DIR" && pwd -P)" ]; then
  echo "sync-track.sh: run this from $TRACK_DIR (current dir: $(pwd))" >&2
  exit 2
fi

log_sync() {
  # $1 = status word, $2 = detail — appended locally; travels with the next commit.
  printf -- '- %s track:%s %s %s\n' "$TS" "$NAME" "$1" "$2" >> "$SESSIONS"
}

flag_conflict() {
  # $1 = short description
  {
    echo ""
    echo "### CONFLICT $TS — track:$NAME"
    echo "$1"
    echo ""
  } >> "$SESSIONS"
  # best-effort: get this visible on origin immediately, even if nothing else follows
  git add "$SESSIONS" >/dev/null 2>&1
  git commit -m "[track:$NAME] flag conflict $TS" >/dev/null 2>&1
  git push origin "HEAD:track/$NAME" >/dev/null 2>&1
  echo "sync-track.sh: CONFLICT — logged to SESSIONS.md NEEDS JACOB (local + best-effort pushed to track/$NAME), stopping. Do not guess a resolution." >&2
}

# 1. commit any pending work on the track branch
if [ -n "$(git status --porcelain)" ]; then
  git add -A
  git commit -m "[track:$NAME] $MSG" >/dev/null
fi
COMMIT_SHA="$(git rev-parse --short HEAD)"

# 2. fetch + rebase + push loop against origin/main (retries on a benign non-fast-forward race;
#    a REAL conflict during rebase is not retried — it's flagged immediately)
attempt=0
while :; do
  attempt=$((attempt+1))
  if ! git fetch origin main >/tmp/sync-track-fetch-$NAME.log 2>&1; then
    log_sync "FAILED" "commit:$COMMIT_SHA fetch-failed-see-tmp-log(network?)"
    echo "sync-track.sh: fetch origin failed (network?) — work is committed locally at $COMMIT_SHA, not synced yet. Log: /tmp/sync-track-fetch-$NAME.log" >&2
    exit 1
  fi

  if ! git rebase origin/main >/tmp/sync-track-rebase-$NAME.log 2>&1; then
    git rebase --abort >/dev/null 2>&1
    flag_conflict "rebase origin/main -> track/$NAME failed (real conflicting content, not just a race). Your prior work IS still committed locally at $COMMIT_SHA. Log: /tmp/sync-track-rebase-$NAME.log"
    log_sync "FAILED" "commit:$COMMIT_SHA rebase-conflict-see-tmp-log"
    exit 1
  fi

  if git push origin "HEAD:main" >/tmp/sync-track-push-$NAME.log 2>&1; then
    log_sync "ok" "commit:$COMMIT_SHA pushed(attempt:$attempt)"
    echo "sync-track.sh: OK — $NAME synced to origin/main ($COMMIT_SHA, attempt $attempt)"
    exit 0
  fi

  if [ "$attempt" -ge "$MAX_RETRIES" ]; then
    log_sync "PARTIAL" "commit:$COMMIT_SHA rebased-locally push-failed-after-$MAX_RETRIES-tries-see-tmp-log"
    echo "sync-track.sh: rebased onto origin/main locally but push failed $MAX_RETRIES times (network, or push-protected). Work is committed+rebased at $COMMIT_SHA — retry sync later. Log: /tmp/sync-track-push-$NAME.log" >&2
    exit 1
  fi
  # otherwise: someone else likely pushed between our fetch and our push — loop and retry
done
