---
status: superseded by ADR-0009
---

# Pre-flight existence check, with `workflow.output.overwrite = true`

The run must never overwrite a CRAM. Nextflow's `workflow.output.overwrite = false` doesn't fail: it
skips an existing target silently (no warning, exit 0, and the file is still listed as published). With
that setting, registration would then record a file this run didn't write. So instead, `main.nf` checks
every CRAM, `.crai` and metrics path the run would publish before any alignment starts, and fails with the
full list. `overwrite` is left `true`, so a file that appears mid-run is replaced (or refused by IAM)
loudly, rather than skipped.

## Consequences

- If a CRAM was published but never registered (for example, its registration failed), the next run for
  that SG stops at the pre-flight check until someone registers or deletes the file.
- Two runs started together for the same SGs can both pass the check. Nothing guards against that.
- The check is plain Groovy, not a task, so it runs again on `-resume`, where the earlier attempt's
  outputs exist. To keep `-resume` usable, a run that passes the check writes the paths it claims to
  `<workDir>/nf-align-claims/<sessionId>.txt`, and claimed paths are allowed to exist. `-resume` reuses
  both the session ID and the work dir, so only the same session's earlier attempts count; a fresh run
  (a new session) checks everything.
