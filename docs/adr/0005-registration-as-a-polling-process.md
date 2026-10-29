# Registration is a per-SG process that polls for the published file

`REGISTER_ANALYSIS` runs once per new CRAM, inside the workflow. Workflow outputs are published
asynchronously: the registration task can start before the copy has begun, and locally it was even seen
reading a partly copied file. So the task polls until both the CRAM and the `.crai` exist (up to
`register_wait_minutes`), then calls `cpg_utils.metamist_registration`, which re-checks that the files
exist. GCS objects only become visible once complete, which makes "exists" a safe signal there. The
registered path is built by the same `cramPath()` function that the output block publishes with, so the
two can't drift apart.

## Considered Options

- **A post-run script that reads the output index file.** It doesn't wait for anything, but it puts
  registration outside the pipeline, where it can be forgotten or run against the wrong run's index.
- **Registering from `onComplete`.** It isn't a process, so it can't be retried or containerised.
