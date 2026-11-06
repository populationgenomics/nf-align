# Workflow outputs can't publish to more than one bucket

The initial design proposal for the CPG-Flow -> NextFlow transition used [Workflow Outputs](https://docs.seqera.io/nextflow/tutorials/workflow-outputs), 
a relatively new NF feature where each workflow has a single section documenting all workflow products
and output locations. Based on NF documentation this is the future direction, with the previous version
(publishDir) directives potentially being phased out.

Using Workflow outputs, the `output {}` block publishes everything under a single `outputDir`. In cloud
environments, an `outputDir` is a bare bucket, not targeting any directory structure within the bucket.

An implication of this is that each workflow writing to a bucket must write all its outputs to a single 
bucket, which does not match the way CPG workflows are typically run. We often use cross-dataset aggregation 
steps, later sharding each dataset's portion of the aggregate into its own bucket for sharing with that 
specific collaborator. See GATK-SV, gCNV, RNA-Seq, and small-variant joint-calling. 

Even workflows which operate on each sample in parallel (Mito, STR) can be operationally simplified by
running across a range of datasets in one run, instead of requiring one run per dataset. 

## Incorrect Assumptions

Claude at some point became convinced that passing absolute paths into the workflow outputs `output` block
could be used to direct outputs to an exact bucket, independent of the `workflow.outputDir`. This was proven
incorrect:

A `path` closure that returns a `gs://` string does **not** publish to that bucket. Instead, the string 
becomes a subdirectory of `outputDir`. This corrects an earlier design statement, which says that path 
closures accept absolute `gs://` paths.

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

- One run can publish to one bucket only, which is `outputDir`.
- A mistaken `gs://` target fails silently. Files land under an odd `gs:` prefix inside the dataset
  bucket, and the run still succeeds.
- Publishing to several buckets would need `publishDir`, a separate copy process, or one run per bucket.
