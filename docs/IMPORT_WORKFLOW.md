# Inspect, check and safely add a table

The MCP import workflow now explains what happens to source rows and records how a joined table was made. It does not change the scientific acceptance rules, remove repeated rows, repair gene names automatically or join a source without the person's instruction.

## 1. Inspect

Call `add_card(file, about)` to inspect sheets. Select a sheet and header row and call it again. The response includes a full-table `profile`: row count, missing/distinct values per column and repeated rows beyond the first. Repeats may be separate measurements; they are reported, not removed.

The original bounded model card stays unchanged. This profile is additional information for the person or client, not a change to recorded model experiments.

## 2. Check and repair

Call `add_check` with stage `mapping`. Beside the verdict and existing link-count preview, `report` now provides:

| Field | Meaning |
|---|---|
| `source_rows`, `accounted_source_rows` | These must match; every input row has one outcome |
| `outcomes` | Rows contributing graph pairs, conditional gaps, explicit skips, unresolved endpoints, pairs outside the graph, or unhandled split values |
| `outcome_examples` | Up to ten source-row references per outcome |
| `issue_counts`, `issue_examples` | All issue counts, with up to ten examples total, including unmatched and multiply matched symbols |
| `exported_edges`, `exported_gaps` | Actual exporter totals when validation reached export; otherwise null |

A row holding `APOE/UNKNOWN` can produce an edge for APOE and still report UNKNOWN as unmatched. A symbol that resolves to multiple graph targets is explicitly reported, while retaining the existing expansion rule. No target is selected by guessing.

Source-row numbers are one-based data rows below the selected header. Sheet-row numbers use the same numbering as the sheet card, after blank rows were skipped by the reader; they are not guaranteed to be physical Excel row numbers. Counts are before exporter merging, so ten rows need not produce ten edges. A row expanded into several outcomes has a combined label, and is still counted once.

Repair suggestions:

- Unmatched symbol: inspect the source value and identifier vocabulary; do not silently replace it with a similar name.
- Multiple symbol matches: inspect all mappings and the source's intended entity; several protein products may be legitimate.
- Explicit skip: confirm the exclusion reason matches the scientific question.
- Unhandled split value: add an explicit rule or skip for that value, then check again.
- Conditional gap: confirm the condition and restoration rule; a failed condition is not the same as a missing identifier.
- Mapping outside the graph: inspect category/identifier choice and whether the needed vocabulary is available.

A refused mapping can still have a diagnostic report. Its predicted contributions are not authorization to import it. A schema-valid mapping can also be scientifically inappropriate.

## 3. Join once, retain evidence

After the person has reviewed the accepted mapping and report and instructed the join, `add_join` validates again and writes:

- an immutable content-addressed Parquet table;
- a JSON receipt binding the source-file digest, selected sheet/header, mapping, study metadata, bundle manifest, resolver state, Biolink schema and validator-code identities;
- an atomically replaced manifest referencing both.

The source file itself is not copied into the receipt or release. Retain it in an appropriate licensed location if later rechecking is needed.

Identical repeat imports verify the stored receipt and table and return `already_present: true`. A conflicting study key, including a DOI-normalisation collision, is refused. Legacy imports without receipts remain readable but cannot be replaced through this API. No existing study IDs or tables are migrated automatically.

An interruption before manifest replacement leaves the previous manifest readable. It may leave an unreferenced immutable file; this is not an imported study and is not automatically deleted. The writer lock serialises imports using this API. Older CLI/build writers do not use that lock and must not run concurrently. Explicit versioned updates and user-directed rollback remain future work.

## 4. Recheck a saved recipe offline

With the source file retained and the same bundle available:

```sh
python -m pharos_gap_kg.add replay-import /absolute/path/to/bundle/pharos319 /absolute/path/to/receipt.json /absolute/path/to/source.xlsx
```

This reads the saved mapping and runs the current validator; it calls no model and writes no import. Changed source or bundle bytes are refused. `same_validator_code` says whether the code identity matches the original. A successful replay establishes a current validation verdict, not byte-for-byte equality of the output or scientific correctness. Changed code/reference state requires review of the new verdict.

## Scope and checks

These capabilities are exposed by the MCP add tools and the offline replay command. The browser's existing walkthrough remains a recorded demonstration, not a new arbitrary-file upload/import UI. Tests cover missing/partial/multiple symbol resolution, skips and conditional gaps, bounded examples, repeat imports, legacy/colliding identities, tampering, concurrent writers, interrupted commit and read-only replay. No new hosted-model evaluation, source import or prediction scoring was run to build this workflow.

Verification on November 9: full `make check` passed all module selftests and 180 pytest tests, including HTTP tests with localhost permission. Two existing kgx-gaps warnings concern duplicate columns in refusal diagnostics. `git diff --check` passed. Tests used temporary import stores; no real source was joined or migrated.
