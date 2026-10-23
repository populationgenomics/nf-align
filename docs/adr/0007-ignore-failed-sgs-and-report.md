# A failed SG is dropped, and the run carries on and exits 0

Per-SG tasks use `errorStrategy 'ignore'` once OOM retries are used up. In a cohort-sized run, one
corrupt input shouldn't throw away hundreds of finished alignments. Registration happens per SG, so a
partial run still leaves Metamist consistent. The dropped SGs and the step they failed at are written to
`alignment_samplesheets/<arGuid>/failed.csv` and logged as a warning. A later run on the same cohorts
picks them up again, because they have no CRAM analysis.

## Considered Options

- **`terminate`**: kills tasks that are still running, and their work is lost.
- **`finish`**: lets running tasks complete but starts nothing new, so the rest of the cohort waits for
  another run.

## Consequences

The run's exit code is 0 even when SGs failed. Anything automated has to read `failed.csv` (or the
warning), not the exit status.
