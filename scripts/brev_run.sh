#!/usr/bin/env bash
# Run ON the instance. Fits the named models one after another, detached, with a log each.
#
#   bash scripts/brev_run.sh tune  mf lightgcn lightgcn_ctx lightgcn_ctx_pruned   # validation only
#   bash scripts/brev_run.sh final mf lightgcn                                    # scores on the test links: once per model,
#                                                                                 # after its setting was pulled back and committed
set -euo pipefail
STAGE="${1:?usage: scripts/brev_run.sh tune|final <model>...}"; shift
[ "$STAGE" = tune ] || [ "$STAGE" = final ] || { echo "stage is tune or final"; exit 1; }
[ "$#" -gt 0 ] || { echo "name at least one model: mf lightgcn lightgcn_ctx lightgcn_ctx_pruned"; exit 1; }
cd "$(dirname "$0")/.."
export PHAROS_GAP_KG_DATA="$PWD/data"
mkdir -p logs
( for m in "$@"; do
    echo "=== $STAGE $m  $(date -u +%FT%TZ)"
    .venv/bin/python -m pharos_gap_kg.model "$STAGE" pharos319 --cutoff 2018 --model "$m" || echo "FAILED: $STAGE $m"
  done
  echo "=== done  $(date -u +%FT%TZ)" ) > "logs/$STAGE.log" 2>&1 &
echo "running detached (pid $!); follow with:  tail -f logs/$STAGE.log"
