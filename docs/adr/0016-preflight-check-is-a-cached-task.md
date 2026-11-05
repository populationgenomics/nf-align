# The pre-flight check is a cached task, so `-resume` skips it

The output pre-flight check ([ADR 0009](0009-never-overwrite-preflight-and-publishdir-no-overwrite.md)) runs
as `PREFLIGHT`, a native (`exec:`) process. Its input is every sample sheet row with every URI the run would
publish, and it passes the rows on only if none of those URIs exists. Nothing downstream starts until it has
passed.

The check has to let a resumed run through. Its first attempt may already have published some SGs' files,
and they are not clashes, because this session wrote them. Nextflow already records which tasks of a session
succeeded, and `-resume` skips them. So once `PREFLIGHT` has passed, a resumed run takes it from the cache and
doesn't check again. If it failed, it isn't cached, and a resumed run checks again.

This replaces the claims file of ADR 0009, `<workDir>/nf-align-claims/<sessionId>.txt`. A check written as a
function in `main.nf` runs again on every resume, and to tell this session's own files from anyone else's it
had to keep that list itself.

## Considered Options

- **A claims file (ADR 0009)**: the passed check wrote down every URI it covered, and a resumed run allowed
  those to exist. It works, but it duplicates what Nextflow's cache already knows.
- **Skip the check whenever `workflow.resume` is set**: simplest, but a run whose first attempt failed the
  check would then be let through on resume, and register a file it didn't write
  ([ADR 0013](0013-only-register-files-this-run-wrote.md)).

## Consequences

- If the task's input changes between attempts (an edited `--samplesheet`, a different `output_root`, or a
  rebuilt sample sheet), the cache misses and the check runs again. It then also finds the files the earlier
  attempt published and fails. That fails safe: launch a fresh run.
- `PREFLIGHT` must stop the run, so its `errorStrategy` is `terminate`, not the per-SG `ignore`.
- It runs inside the Nextflow process, so it needs no VM, but it checks every URI one at a time, as before.
