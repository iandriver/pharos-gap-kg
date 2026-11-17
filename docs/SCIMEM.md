# Project memory and build coordination

Hypothesary/scimem project: `pharos-gap-kg`. Root: `/Users/iandriver/Downloads/pharos-gap-kg`. Home: `/Users/iandriver/.scimem`.

The project indexes README, docs/*.md and findings/*.md. Input files and secrets are not included by this recipe. The October 8 user request authorized project documentation and agent communication.

With MCP, pass `project="pharos-gap-kg"` on every call. Start with `start(q=...)`, then `brief(q=...)` before proposing work and `state()` for the card. Read and handle `inbox()` notes, then acknowledge only those handled. Open a task and pass its id to writes. Check draft conclusions with `check(text=...)`. Post build updates to `project:pharos-gap-kg`, naming exact code/input identities, artifacts, checks, limits and next owner. Close only completed tasks. Preserve existing owners and do not infer process liveness from an idle/active card.

CLI on this machine:

```sh
/Users/iandriver/Downloads/scimem/.venv/bin/scimem --home /Users/iandriver/.scimem --project pharos-gap-kg --writer codex start 'current build question'
```

Use the same prefix for `brief`, `state`, `task`, `inbox`, `check` and `post`. CLI executable is named scimem here; its help calls it hypothesary. Refresh documents through project `start` or `rebuild` after changes. Do not use VCC for this entry.

Plans and findings are indexed notes, not automatically validated results. Attach measured result claims and explicit dependencies when a test completes; a document count is not a competition score.
