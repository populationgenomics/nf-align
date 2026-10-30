# Never overwrite: a pre-flight check, plus `publishDir overwrite: false`

Our buckets don't allow existing files to be overwritten, so the pipeline never tries to. Every
`publishDir` sets `overwrite: false`. On Nextflow 26.10.4 that setting leaves an existing target in place
silently: there's no warning and the exit code is 0. A clash would therefore go unnoticed, and
registration would then record a file this run didn't write. So, before any alignment, `main.nf` checks
every CRAM, `.crai` and metrics URI the run would publish, in each SG's own dataset bucket. If any exist,
it fails and lists them all.

Before the sample sheet is built, there's also a cheap check of the run-records directory
(`alignment_samplesheets/<arGuid>/` in the analysis dataset bucket). It fails the run if that directory
already holds a sample sheet, `crams.csv` or `failed.csv`. `arGuid` is unique in production, but it's set
by hand while testing.

This supersedes [ADR 0002](0002-preflight-check-instead-of-overwrite-false.md), which left overwrite on so
that a clash in the middle of a run would fail loudly.

## Consequences

- **`-resume`:** a run that passes the output check records the paths it claims in
  `<workDir>/nf-align-claims/<sessionId>.txt`. On `-resume`, those paths are allowed to exist, and
  `publishDir` skips them when it republishes cached tasks. The run-records check is skipped when
  `workflow.resume` is set.
- **Unregistered CRAM:** a CRAM that was published but never registered (for example, because
  registration failed) stops the next run for that SG at the pre-flight check, until someone registers or
  deletes the file.
- **Accepted gap:** two runs launched close together for the same SG can both pass the check. The second
  run's publish is then skipped silently, and its registration records the first run's CRAM under the
  second run's `ar-guid`. The CRAM itself is valid, and the sample sheet builder takes the latest CRAM
  analysis.
