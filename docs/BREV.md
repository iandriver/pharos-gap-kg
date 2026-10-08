# Fitting the models on a rented GPU (NVIDIA Brev)

The four models of `docs/PROTOCOL.md` ("Fitted models") are fitted by `pharos_gap_kg.model`. The two
without context run on a laptop. The two with the context graph, and the whole grid of settings, are
what a GPU is for. This page is the route from a Brev instance to a committed result.

Nothing here publishes anything. Code and data go to a machine you launched, over SSH.

## What the instance needs

| | size | where it comes from |
|---|---|---|
| the code and `docs/data` | 2 MB | `scripts/brev_push.sh` |
| the frozen split | 2 MB | same |
| the context tables, pruned and unpruned | 9 MB | same |

The Pharos dumps and the extracted tables stay on this machine. The last stage, `model score`, needs
them (the ontology, gene positions) and is run here on the score matrix the instance sends back.

Any single GPU with 16 GB of memory is more than enough: the largest graph has under 100,000 nodes and
about a million edges. Pick the cheapest one offered, with 8 CPUs or more, because validation is ranked on
the CPU.

## Once

Installing the CLI, logging in and launching are yours to do.

```bash
brew install brevdev/homebrew-brev/brev
```

```bash
brev login
```

Create the instance in the Brev console, or with `brev create <name> --gpu <type>` using a type from
the console's list. Then:

```bash
brev refresh
```

```bash
ssh <name> nvidia-smi
```

`brev refresh` writes the SSH entry that the last command and the scripts below use. If `ssh <name>`
does not connect, `brev shell <name>` opens a terminal on the instance, and the scripts accept any SSH
host (`user@address`) in place of the name. The Brev quickstart documents `login`, `create`, `shell`,
`refresh` and `list`; for stopping and deleting an instance use the console or `brev --help`.

## Each run

1. **Send the code and the small data, and set the instance up.**

   ```bash
   scripts/brev_push.sh <name>
   ```

   ```bash
   ssh <name> 'cd pharos-gap-kg && bash scripts/brev_setup.sh'
   ```

   The setup makes a Python 3.12 environment, installs the package with torch, prints the GPU, checks
   that the split on the instance is the frozen one, and runs the model self-test.

2. **Tune.** Validation only. No test link is read.

   ```bash
   ssh <name> 'cd pharos-gap-kg && bash scripts/brev_run.sh tune mf lightgcn lightgcn_ctx lightgcn_ctx_pruned'
   ```

   It runs detached and writes `logs/tune.log` on the instance.

3. **Bring the settings back and commit them.**

   ```bash
   scripts/brev_pull.sh <name>
   ```

   Each model's chosen setting is `docs/data/model_<name>_2018.json`. The protocol requires it to be
   committed before the model is scored on the test links.

4. **Score on the test links, once per model.** Push again so the instance has the committed
   settings, then:

   ```bash
   ssh <name> 'cd pharos-gap-kg && bash scripts/brev_run.sh final mf lightgcn lightgcn_ctx lightgcn_ctx_pruned'
   ```

   `final` refuses a model it has already scored, and a setting from a partial grid.

5. **Pull the score matrices and run the controls here.**

   ```bash
   scripts/brev_pull.sh <name>
   ```

   ```bash
   python -m pharos_gap_kg.model score pharos319 --cutoff 2018 --model lightgcn
   ```

   `score` recomputes the neighbours and degree baselines beside the model and refuses to report
   unless they equal the committed ones.

6. **Stop the instance.** It bills while it runs.

## What to expect

Measured on this laptop's CPU while it was busy with other work, one setting each, validation only
(the neighbours baseline scores 0.0712 on the same validation links):

| model | settings in the grid | one setting, CPU | first setting, validation |
|---|---|---|---|
| `mf` | 8 | 62 s | 0.0543 |
| `lightgcn` | 24 | 85 s | 0.0658 |
| `lightgcn_ctx`, `lightgcn_ctx_pruned` | 24 each | 573 s | 0.0618 |

Those are single settings from a dry run that chose nothing, not results. They say the first two
models can be tuned on a laptop in under an hour, and that the GPU's job is the two context models
(about four hours each on this CPU) and the five-seed refits. None of the three first settings beats
the neighbours baseline on validation.

`lightgcn_ctx_pruned` reads `disease_links_pruned.parquet`, so push after the pruning rules have been
applied (`python -m pharos_gap_kg.applycut write pharos319`).
