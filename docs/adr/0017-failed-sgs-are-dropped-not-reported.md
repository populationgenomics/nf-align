# A failed SG is dropped, the run exits 0, and recovery is `-resume`

Per-SG tasks still use `errorStrategy 'ignore'` once OOM retries are used up, so one corrupt input doesn't
throw away the rest of a cohort run ([ADR 0007](0007-ignore-failed-sgs-and-report.md)). But the run no
longer reports its dropped SGs: there is no `failed.csv`, and no warning that lists them.

Working out the first step each SG didn't get through took a dedicated part of `main.nf`: every stage
emitted an event, and the events were grouped per SG and compared with a list of steps that changed with
`--register` and `--cram_qc`. The same information is in Seqera Platform, which lists every failed task by SG.
And the run is resumable: `-resume` re-runs only the tasks that failed, and those that come after them. The
pre-flight check is cached ([ADR 0016](0016-preflight-check-is-a-cached-task.md)), so the files the first
attempt published don't stop it.

This supersedes ADR 0007, keeping its `ignore` strategy and dropping its report.

## Considered Options

- **Report each SG's failed step in `failed.csv` (ADR 0007)**: machine-readable, but costly to keep right,
  and duplicates what Seqera shows.
- **List the SGs that didn't finish, without the step**: still needs the per-SG bookkeeping in `main.nf`.

## Consequences

- The run's exit code is 0 even when SGs failed, and nothing in the run records says so. Check the run in
  Seqera, or compare `crams.csv` with the sample sheet. `crams.csv` lists the CRAMs that were aligned and
  published, whether or not their registration or QC then succeeded.
- **Resume, don't relaunch.** A resumed run picks up an SG that failed at fingerprinting or registration, and
  registers the CRAM it published. A new run (a new session) doesn't: it finds that published, unregistered
  CRAM and stops at the pre-flight check, as before
  ([ADR 0013](0013-only-register-files-this-run-wrote.md)).
