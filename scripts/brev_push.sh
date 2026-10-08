#!/usr/bin/env bash
# Copy the code and the small files a model needs to a GPU instance you launched, over SSH.
# Nothing is published. The raw Pharos dumps and the extracted tables stay here: fitting needs
# only the frozen split (2 MB) and the context tables (9 MB).
#
#   scripts/brev_push.sh <ssh-host>        # the instance name, once `ssh <name>` connects
set -euo pipefail
HOST="${1:?usage: scripts/brev_push.sh <ssh-host>}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
for f in data/eval/gene_trait_2018/train.parquet data/eval/gene_trait_2018/test.parquet data/eval/gene_trait_2018/unscored.parquet \
         data/context/pharos319/partners.parquet data/context/pharos319/pathways.parquet data/context/pharos319/go_terms.parquet \
         data/context/pharos319/compounds.parquet data/context/pharos319/disease_links.parquet data/context/pharos319/disease_links_pruned.parquet; do
  [ -f "$f" ] || { echo "missing $f: build it here first (README, 'Run it')"; exit 1; }
done
ssh "$HOST" 'mkdir -p ~/pharos-gap-kg'
# code and the committed summaries; data/ is sent separately and never deleted on the instance
rsync -az --delete --exclude .git --exclude .venv --exclude '/data' --exclude '/logs' --exclude __pycache__ --exclude .pytest_cache ./ "$HOST:pharos-gap-kg/"
rsync -azR data/eval/gene_trait_2018/train.parquet data/eval/gene_trait_2018/test.parquet data/eval/gene_trait_2018/unscored.parquet \
      data/context/pharos319/partners.parquet data/context/pharos319/pathways.parquet data/context/pharos319/go_terms.parquet \
      data/context/pharos319/compounds.parquet data/context/pharos319/disease_links.parquet data/context/pharos319/disease_links_pruned.parquet \
      "$HOST:pharos-gap-kg/"
echo "pushed to $HOST:pharos-gap-kg; next:  ssh $HOST 'cd pharos-gap-kg && bash scripts/brev_setup.sh'"
