# One dataset and one sequencing type per run

A run takes `params.dataset` and `params.sequencing_type`. Its Metamist query is scoped to that dataset's
project, and `outputDir` is the dataset bucket, so the output block can use cpg-flow's fixed relative
layout (`cram/`, `exome/cram/`, `qc/...`). An SG from another project fails the sample sheet step, and
SGs of the other sequencing type are skipped (or, when named with `--sgids`, fail it). A cohort that spans datasets or sequencing types needs one
run per (dataset, type).

## Considered Options

- **Several datasets per run, with absolute `gs://` paths per SG.** Testing on Nextflow 26.04 showed that
  `output {}` path *closures* do accept absolute paths (static strings are rejected), so this is possible.
  It was rejected because registration, `accessLevel` and the pre-flight check would all have to become
  per-SG, and a single dataset per run keeps permissions and billing in one place.
