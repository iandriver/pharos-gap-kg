#!/usr/bin/env bash
# Bring back what the instance produced: each model's chosen setting (a small file to commit), its
# tuning record, and after `final` its score matrix.
#
#   scripts/brev_pull.sh <ssh-host>
set -euo pipefail
HOST="${1:?usage: scripts/brev_pull.sh <ssh-host>}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p data/eval/gene_trait_2018/models logs/brev
rsync -az --include 'model_*_2018.json' --exclude '*' "$HOST:pharos-gap-kg/docs/data/" docs/data/
rsync -az "$HOST:pharos-gap-kg/data/eval/gene_trait_2018/models/" data/eval/gene_trait_2018/models/
rsync -az "$HOST:pharos-gap-kg/logs/" logs/brev/ || true
git status --short docs/data | sed 's/^/  /'
echo "pulled. Read logs/brev/, then commit docs/data/model_*_2018.json before running 'final'."
