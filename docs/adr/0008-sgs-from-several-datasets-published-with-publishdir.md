# A run covers SGs from several datasets, each published and registered in its own dataset

A run is launched under one analysis dataset (`params.dataset`), but its cohorts or SG IDs can span many
datasets. Each SG's dataset is its sample's Metamist project, with any `-test` removed. The cohort it came
through doesn't matter. Its CRAM, `.crai` and metrics are published with `publishDir` to that dataset's
bucket, in the usual per-dataset layout (`[exome/]cram/`, `[exome/]qc/...`), and the CRAM is registered in
that dataset's Metamist project. Per-SG tasks carry that dataset in their `dataset` resource label. Only
the run records (`alignment_samplesheets/<arGuid>/`) go to the analysis dataset's bucket.

The reason is operational. Workflow outputs can only publish under a single `outputDir`
([clarification](../clarifications/workflow-outputs-single-bucket.md)), which would mean one launch per
dataset: about 35 separate launches for a typical batch. `publishDir` can choose a destination for each
task, so a single launch can cover them all.

This supersedes [ADR 0001](0001-one-dataset-one-sequencing-type-per-run.md). The rule of one sequencing
type per run stays.

## Considered Options

- **One run per dataset, using workflow outputs** (ADR 0001). This was the original design, chosen to use
  Nextflow's newest publishing mechanism. It was dropped because of the number of launches.
- **Workflow outputs with absolute `gs://` paths.** This isn't possible. ADR 0001 said such paths were
  accepted, but a `gs://` string is resolved as a subdirectory of `outputDir`, and nothing reaches the
  other bucket.
- **An allow-list of datasets a run may write to.** Rejected. The SG's own project picks the bucket, so a
  wrong cohort can only bring in unintended SGs. It can't send data to the wrong place.

## Consequences

- The analysis dataset's service account must be able to read every dataset its SGs belong to, and write
  to their buckets and Metamist projects.
- `accessLevel` still applies to the whole run. Any SG whose project doesn't match it (a production
  project in a test run, or the other way round) fails the sample sheet step.
- Output paths are worked out once in `main.nf` and carried in each SG's `meta`. `publishDir`, the
  pre-flight check and registration all read them from there.
- A prepared `--samplesheet` must have a `dataset` column.
