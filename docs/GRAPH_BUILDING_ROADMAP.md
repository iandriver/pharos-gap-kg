# Making graph building easier

November 9, 2026. Human evaluation remains useful, but it is not a dependency for this engineering work. Do not claim usability improvement until it is measured.

The current system can propose a sheet and mapping, validate it, preview link counts and join a later source. The browser walkthrough makes these steps accessible. The next goal is to make every refusal repairable and every accepted import inspectable, reproducible and reversible.

## First increment: inspect before mapping

The MCP selected-sheet `add_card` response now includes a read-only `profile`: full-table row count, missing and distinct values by column, and repeated rows beyond the first. It uses the same missing-value convention as validation and does not remove anything. The existing bounded model card is unchanged; scientific interpretation and acceptance rules are unchanged. The profile helps identify an empty identifier column, mostly missing measurements, or repeated records before asking a model to map them.

## Work that does not need a participant

| Priority | Capability | How to test without a human session |
|---|---|---|
| 1 | Explain each unmatched identifier and affected source row; show exact, ambiguous and unmatched cases separately | Fixed alias/obsolete/unknown/many-to-many fixtures. Suggest aliases with provenance; never silently replace an ambiguous identifier |
| 1 | A complete dry-run change report: source rows, exported edges, gaps, dropped rows and reasons; sample each outcome | Account for every input row, including intentional gaps, footnotes and rows expanded into several edges; counts must reconcile. Schema acceptance must remain separate from scientific appropriateness |
| 1 | Immutable import receipt and safe updates | Bind source bytes, selected sheet/header, mapping, validator version, bundle identity and outputs. Test DOI-key collisions, repeat imports, interrupted writes and rollback; refuse destructive replacement by default |
| 2 | Save and replay mapping recipes | Rebuild an import offline from source+recipe, with identical identifiers and output; changed source/schema yields an explicit difference, not stale cache reuse |
| 2 | Guided repair in the browser | Each refusal names the field, offending rows, rule and suggested next action; test known failures end-to-end. Let the person select another sheet without restarting the workflow |
| 2 | Schema-change preview for a later source version | Report added/removed columns, changed units/types and evidence dates, affected outputs and expected gap changes before accepting the update |
| 3 | Portable source adapters beyond Pharos | Define a minimal adapter contract for identifiers, evidence, dates, units and provenance; prove it on a small independent public source without assuming Pharos target IDs everywhere |

## Evaluation contract

Use small offline regression fixtures for parser/validator faults and a separately frozen unseen-source set for model generalisation. Measure correct sheet selection, identifier coverage, ambiguity handling, row accounting, mapping agreement, reproducibility, repair attempts and API cost separately. A known study replay tests reproducibility, not generalisation. A code-valid mapping can still be scientifically inappropriate.

Current broader opportunities include historical symbols missed by exact lookup and discordance between paper-named genes and positional Catalog assignments. The right response is an explicit mapping method and provenance, not choosing whichever list yields a higher overlap.

Keep study-specific publication licences and source evidence semantics intact. Do not infer that a gene named near a GWAS variant is causal. Preserve spent prediction tests; this roadmap does not authorize new scientific scoring or external API spending.

## Suggested next build

Extend the existing dry run with an unresolved-row report, then add import receipts and atomic updates. These steps reduce debugging effort and protect existing graphs regardless of whether a person has tested the interface. The human task packet can later measure whether the same reports help someone make a correct decision.

## Implemented November 9

The next vertical slice is built: source-row accounting and unresolved/multiply matched symbol diagnostics in MCP `add_check`; immutable content-addressed table/receipt commits and verified idempotence in MCP `add_join`; and a read-only saved-recipe replay command. See [IMPORT_WORKFLOW.md](IMPORT_WORKFLOW.md) for field definitions, commands and limitations. Conflicting/legacy imports are refused rather than silently migrated. Versioned updates, rollback UI, browser repair controls and portable adapters remain future work; no human-usability gain has yet been measured.
