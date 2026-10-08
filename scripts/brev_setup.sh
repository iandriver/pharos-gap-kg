#!/usr/bin/env bash
# Run ON the instance, once, from the pushed checkout:   cd ~/pharos-gap-kg && bash scripts/brev_setup.sh
# Makes a Python 3.12 environment with the package and torch, then checks the GPU and the code.
set -euo pipefail
cd "$(dirname "$0")/.."
export PHAROS_GAP_KG_DATA="$PWD/data"
python3 -m pip install --user --quiet uv            # from PyPI; uv then fetches Python 3.12 if the image has an older one
export PATH="$HOME/.local/bin:$PATH"
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[test,eval,model]"
echo "--- GPU"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || echo "no GPU visible: the models will fit on the CPU, slowly"
.venv/bin/python -c "import torch; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
echo "--- the split is the frozen one"
.venv/bin/python -c "from pharos_gap_kg import split; p = split.load(2018); print({k: len(v) for k, v in p.items()})"
echo "--- self-tests"
.venv/bin/python -m pharos_gap_kg.model --selftest
echo "ready; next:  bash scripts/brev_run.sh tune mf lightgcn lightgcn_ctx lightgcn_ctx_pruned"
