# Workflow outputs can't publish to more than one bucket

The `output {}` block publishes everything under a single `outputDir`. A `path` closure that returns a
`gs://` string does **not** publish to that bucket. Instead the string becomes a subdirectory of
`outputDir`. This corrects the "Considered Options" in
[ADR 0001](../adr/0001-one-dataset-one-sequencing-type-per-run.md), which says that path closures accept
absolute `gs://` paths.

## What happens

- **Relative string**: resolved under `outputDir`, as documented.
- **Absolute local path** (`/x/y`): escapes `outputDir`. This is a side effect of Java path joining and
  isn't documented.
- **`gs://bucket/x/`**: becomes `<outputDir>/gs:/bucket/x/`. Nothing is written to `gs://bucket`.
- **A `Path` object**, such as `file('gs://...')`: the run fails with "Invalid `path` directive".

## Evidence

- **Docs.** [Workflow outputs](https://docs.seqera.io/nextflow/workflow#outputs) say only "Specify the
  publish path relative to the output directory". Every closure example returns a relative string. The
  page doesn't mention absolute or remote targets.
- **Source.** In `PublishOp.groovy`, the closure result has to be a string. That string goes into
  `outputDir.resolve(String)` in v26.04.0, v26.04.9 and v26.10.4. Nextflow doesn't override `resolve`.
  Java treats a string that doesn't start with `/` as relative, and `gs://...` doesn't. On a GCS
  `outputDir`, google-cloud-nio's `CloudStoragePath.resolve(String)` also parses the string as a relative
  path. Unlike `getPath()`, it doesn't reject a `gs:` scheme.
- **Experiment.** ✅ (Nextflow 26.10.4, local `outputDir`) A closure returning `'gs://some-bucket/abs/'`
  published to `results/gs:/some-bucket/abs/a.txt`. The "Outputs" summary listed it as
  `gs:/some-bucket/abs/a.txt` and raised no error. Behaviour on a GCS `outputDir` is inferred from the
  library source and hasn't been run against a real bucket.

## Consequences

- One run can publish to one bucket only, which is `outputDir`. That supports the single-dataset rule in
  ADR 0001 for a stronger reason than the one the ADR gives: the multi-dataset option isn't possible with
  workflow outputs at all.
- A mistaken `gs://` target fails silently. Files land under an odd `gs:` prefix inside the dataset
  bucket, and the run still succeeds.
- Publishing to several buckets would need `publishDir`, a separate copy process, or one run per bucket.
