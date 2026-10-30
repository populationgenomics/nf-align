# nf-align: design brief

Aligns the short-read sequencing groups of a set of Metamist cohorts, or a list of sequencing groups,
which can come from several datasets. The output for each is a deduplicated, coordinate-sorted CRAM,
published to and registered in the SG's own dataset. It's a Nextflow (≥ 25.10, strict syntax) pipeline that runs on Google
Batch against data in GCS, and it can be tested offline with stub runs.

The pipeline adapts the alignment half of the Nextflow port in `cpg-flow-align-genotype/nextflow/`.
It aligns with **bwa-mem3** (that port used DRAGMAP), and adds an in-pipeline sample sheet step,
Metamist registration and a pre-flight output check. Genotyping is left out for now, but it is the next
addition, and the shapes below are chosen so it can be wired in from `main.nf` alone.

Terms are defined in [GLOSSARY.md](GLOSSARY.md). Why things are the way they are is recorded in
[adr/](adr/). [RISKS.md](RISKS.md) lists what may bite later, and
[AMBIGUOUS-CHOICES.md](AMBIGUOUS-CHOICES.md) lists the calls that could reasonably have gone the other way.
The full design conversation is in [design-session-transcript.md](design-session-transcript.md).

## Scope of a run

One run is launched under one **analysis dataset** and covers **one sequencing type**. Its SGs can come
from any number of datasets. Each SG is published to its own dataset's bucket and registered in its own
dataset's Metamist project
([ADR 0008](adr/0008-sgs-from-several-datasets-published-with-publishdir.md)).

An SG's dataset is its sample's Metamist project, with any `-test` removed. The cohort it came through
doesn't matter: a cohort in `alpha` that holds an SG from `beta` sends that SG's CRAM to `beta`'s bucket.

| param | meaning |
|---|---|
| `dataset` | the analysis dataset (no `-test` suffix): the run is launched under it, it pays for the run, and it holds the run records. Its service account must be able to see every SG's dataset and write to their buckets and Metamist projects |
| `sequencing_type` | `genome` or `exome` |
| `accessLevel` | `full` / `standard` (equivalent), or `test`. Applies to the whole run: every SG's Metamist project must be `<ds>` (full) or `<ds>-test` (test), and its bucket is `gs://cpg-<ds>-main` or `-test` |
| `arGuid` | the analysis-runner GUID, used in resource labels, the per-run record directory and the registered analyses. Mandatory |
| `cohorts` **or** `sgids` **or** `samplesheet` | exactly one: Metamist cohort IDs, Metamist SG IDs (queried directly), or a prepared sheet that skips the Metamist step |
| `output_root` | null by default, so each dataset's outputs go to its bucket. If set (the `test` profile sets `results`), they go to `<output_root>/<ds>/` instead, with the same layout below that |
| `ref_fasta` | the alignment reference. Its `.fai`, bwa-mem3 index and `.alt` (`<ref_fasta>.{fai,bwt.2bit.64,amb,ann,pac,alt}`) must sit beside it, and `ref_fasta_dict` is its alt-annotated `.dict`. Aligned alt-aware ([ADR 0010](adr/0010-gatk-unmasked-reference-alt-aware.md)) |
| `align_container`, `cpg_container` | the images, with defaults in `nextflow.config`. The run fails at launch if either is unset, in every profile |

There is no `outputDir`: every destination is worked out per SG (see Outputs). The launcher supplies
`-work-dir` and `google.project`; the config never sets them.

## Flow

```
run-records check ─ --cohorts | --sgids ─ BUILD_SAMPLESHEET ─┐   (or --samplesheet)
                                                             ▼
          sample sheet ─ parse, per-SG destinations into meta ─┬─ aligned ────────────────────────────▶ (all_crams, for genotyping)
                                                               ├─ aligned, no somalier ─┐
                                                               └─ to align ─────────────┴─ pre-flight ─┬─ ALIGN ─ publishDir: CRAM, .crai, metrics → SG's dataset bucket
                                                                                                       └─ SOMALIER_EXTRACT (new and gapfilled CRAMs) ─ publishDir: .somalier → SG's dataset bucket
                                                       then REGISTER_ANALYSIS (cram and somalier, per SG, SG's dataset)
                                                       and  WRITE_RUN_RECORDS (crams.csv, failed.csv → analysis dataset bucket)
```

### 1. Sample sheet (`BUILD_SAMPLESHEET`, `bin/build_samplesheet.py`)

- One GraphQL query, using the top-level `cohorts(id:)` or `sequencingGroups(id:)` lookup, so it isn't
  limited to one project. It fetches cohort status and the cohorts' SGs (`--cohorts`), or the named SGs
  directly (`--sgids`), with each SG's project, assay meta and completed, active `cram` and `somalier`
  analyses. The
  script takes `--access-level`, not a project.
- Cohorts that are missing (or not visible to the run's credentials), inactive or empty fail the run. So
  does an SG whose project doesn't match the access level: a production project in a test run, or a
  `-test` project in a full run. Every offender is reported in one error.
- `--cohorts`: SGs are de-duplicated across cohorts. SGs that aren't short-read, or aren't
  `sequencing_type`, are skipped with a warning.
- `--sgids`: each SG was asked for by name, so one that's missing (or not visible), or isn't short-read
  `sequencing_type`, fails the run instead of being skipped.
- Read locations are used **exactly as Metamist stores them**.
- Each SG gets one row per alignment unit, with a `kind` (glossary): `fastq` and `fastq_ora` (one row
  per lane), `cram` or `bam` (one row), or `aligned` (it already has a CRAM analysis; one row pointing at
  that CRAM). An `aligned` row's `somalier` is the latest `somalier` analysis's output, or empty if it has
  none, and then the CRAM is fingerprinted (**gapfilled**).
- Each row records the SG's `dataset`: its project, with any `-test` removed.
- A run with no in-scope SG fails.
- The sheet is published to `alignment_samplesheets/<arGuid>/samplesheet.tsv` in the analysis dataset
  bucket. A `--samplesheet` input is not republished, and must have a `dataset` column.

Columns: `sg_id, dataset, kind, n_units, unit, r1, r2, reads, reads_index, reads_ref, somalier`. For a CRAM,
`reads_ref` is the reference it was written against, and is required: a CRAM is never decoded with `ref_fasta`.
For ORA, it's the ORA reference tarball. For an `aligned` row, it's the `ref_fasta` its analysis's meta records
(nf-align records it); if empty, `aligned_ref_fasta` (default: the DRAGEN-masked FASTA the CRAMs made before
nf-align used, ADR 0006) is assumed. A `--samplesheet` input without a `somalier` column has every `aligned`
row gapfilled.

### 2. Pre-flight

The pipeline never overwrites a file ([ADR 0009](adr/0009-never-overwrite-preflight-and-publishdir-no-overwrite.md)).
`publishDir` uses `overwrite: false`, which skips an existing target silently, so clashes are caught
up front:

- **Run records**: before the sample sheet is built, the run fails if
  `alignment_samplesheets/<arGuid>/` in the analysis dataset bucket already holds a `samplesheet.tsv`,
  `crams.csv` or `failed.csv`. This guards against an `arGuid` reused by hand while testing. It's skipped
  on `-resume`.
- **Outputs**: before anything is aligned, every CRAM, `.crai`, `.somalier` and metrics URI that the run would
  publish is checked in its SG's dataset bucket. If any exist, the run fails and lists them all. For a gapfilled
  `aligned` row only the `.somalier` is checked, and other `aligned` rows are not checked: for those, an existing
  CRAM is the expected state.

A run that passes the check records the outputs it claims in `<workDir>/nf-align-claims/<sessionId>.txt`.
`-resume` keeps the session ID, so a resumed run lets those outputs exist and republishes them from the
cached tasks (or, where they're already there, skips them). Anything else that already exists still
fails the run.

Two runs launched close together for the same SG can both pass the check. That gap is accepted
(ADR 0009).

### 3. ALIGN subworkflow

| input | path |
|---|---|
| CRAM/BAM without an index | `SORT_INDEX`, then sharded like an indexed input |
| CRAM/BAM | `realignment_shards` bazam shards ([ADR 0003](adr/0003-read-name-hash-shards-not-regions.md)), defaulting to 6 for genomes and 1 for exomes |
| FASTQ / FASTQ.ora | one unit per lane; ORA files are streamed through `orad` into bwa-mem3 via FIFOs |

- **One unit** goes to `ALIGN_TO_CRAM`: `align_stream.sh | dedup_to_cram.sh` in a single task.
- **Several units** go to `ALIGN_SHARD` (aligned, then name-sorted to BAM). They are then grouped per SG
  with `groupKey`, and `MERGE_TO_CRAM` runs `samtools merge -n | dedup_to_cram.sh`. Duplicates are
  marked once per SG, after the merge ([ADR 0004](adr/0004-mark-duplicates-once-per-sg-after-merge.md)).
- **Alignment:** `bwa-mem3 mem -t <⅔ cpus> -K 100000000 --bam -R '@RG\tID:<sg>\tSM:<sg>\tPL:ILLUMINA'`.
  CRAM/BAM input is piped in interleaved from bazam (`-p -`).
- **`dedup_to_cram.sh`:** `dupblaster | samtools sort -u | samtools view -O cram,version=3.0 --write-index`,
  then the dupblaster duplicate metrics are renamed to `<sg>.markduplicates-metrics`.
- **CPU split:** `bin/threads.sh` is the only place the per-tool ratios live:

  | CPUs | bwa-mem3 | bazam / merge input | sort while aligning | final CRAM write |
  |---|---|---|---|---|
  | 16 (align) | 10 | 3 | 2 | 15 |
  | 8 (merge, `--no-aligner`) | – | 4 | 3 | 7 |

- **Guards:** every script runs under `set -euo pipefail`. `check_alignments.sh` runs `samtools quickcheck`
  and requires at least one record, on every shard BAM and every final CRAM.
- **Emits:** `cram` (newly aligned only), `markdup`, and `all_crams` (new crams plus `aligned` rows, with
  no tasks of their own), which is ready for genotyping.

### 4. Outputs

Everything is published with `publishDir`, `mode: 'copy'` and `overwrite: false`. There is no
`output {}` block: it can only publish under a single `outputDir`
([clarification](clarifications/workflow-outputs-single-bucket.md)).

When `main.nf` parses a sample sheet row, it works out the SG's destinations and puts them in `meta`
next to the rest of its inputs:

- `dataset`, and its Metamist project (`<ds>` or `<ds>-test`)
- its output root: `gs://cpg-<ds>-{main|test}`, or `<output_root>/<ds>`
- the relative paths of its CRAM, `.crai`, somalier fingerprint and metrics

The `publishDir` settings in config only apply these values, and the pre-flight check and registration
read the same `meta`, so publishing and registration can't drift apart. `main.nf` stays the one place the
layout is defined:

| target | path under the SG's output root (`exome/` prefix for exomes) |
|---|---|
| CRAM + index | `[exome/]cram/<sg>.cram`, plus `.crai` |
| somalier fingerprint | `[exome/]cram/<sg>.cram.somalier` (cpg-flow's `CramPath.somalier_path`) |
| duplicate metrics | `[exome/]qc/markduplicates_metrics/<sg>.markduplicates-metrics` |

Run records go under the analysis dataset's output root, at
`alignment_samplesheets/<arGuid>/`:

| file | written by |
|---|---|
| `samplesheet.tsv` | `BUILD_SAMPLESHEET` |
| `crams.csv` (`sg_id, dataset, cram, crai`, the published URIs) | `WRITE_RUN_RECORDS`, a native (`exec:`) process that gets the collected records, sorted by SG |
| `failed.csv` (`sg_id, dataset, failed_at`) | `WRITE_RUN_RECORDS` |

Only newly aligned SGs are published, plus the fingerprints of gapfilled SGs. A gapfilled fingerprint goes
to the standard layout, `[exome/]cram/<sg>.cram.somalier`, even when the existing CRAM is somewhere else.

Each per-SG task's `dataset` resource label is `meta.dataset`, so its cost goes to the SG's dataset.
`BUILD_SAMPLESHEET` and `WRITE_RUN_RECORDS` are labelled with the analysis dataset.

### 5. Registration (`REGISTER_ANALYSIS`)

When `register` is true, this runs once per new CRAM and once per fingerprint (new or gapfilled), in the SG's
own Metamist project. It gets the output URI from `meta`, and polls for up to `register_wait_minutes` until the
output (and a CRAM's `.crai`) exist, because publishing is asynchronous
([ADR 0005](adr/0005-registration-as-a-polling-process.md)). It then runs:

```
python3 -m cpg_utils.metamist_registration --project <SG's dataset>[-test] --output <cram> --type cram \
    --secondary index=<cram>.crai --sgs <sg> --meta sequencing_type=<type> ar-guid=<arGuid> ref_fasta=<ref>
python3 -m cpg_utils.metamist_registration --project <SG's dataset>[-test] --output <cram>.somalier \
    --type somalier --sgs <sg> --meta sequencing_type=<type> ar-guid=<arGuid> ref_fasta=<CRAM's ref>
```

Metrics are not registered.

### 6. Failures

- **Per-SG tasks:** an exit code of 137 or 140 (OOM) is retried twice, with memory scaled by the attempt
  number. Any other failure is **ignored**: that SG is dropped and the rest carry on
  ([ADR 0007](adr/0007-ignore-failed-sgs-and-report.md)).
- **Registration** is retried once and then ignored.
- **Spot preemption** is retried by `google.batch.maxSpotAttempts = 3`.
- **`BUILD_SAMPLESHEET`**, **`WRITE_RUN_RECORDS`**, the run-records check and the pre-flight check fail the
  whole run.
- **Reporting:** dropped SGs are written to `alignment_samplesheets/<arGuid>/failed.csv` with the step
  they failed at (`align`, `somalier`, `register` for the CRAM, or `register_somalier`), and logged as a warning. **The run still exits 0.**

## Resources (`conf/base.config`)

| process | CPUs | memory | disk (genome / exome) | machine type, by attempt |
|---|---|---|---|---|
| `ALIGN_TO_CRAM`, `ALIGN_SHARD` | 16 | 48 GB × attempt | 400 / 100 GB | `n2-standard-16`, `n2-highmem-16`, `n2-highmem-32` |
| `MERGE_TO_CRAM` | 8 | 24 GB × attempt | 400 / 100 GB | `n2-standard-8`, `n2-highmem-8`, `n2-highmem-16` |
| `SORT_INDEX` | 8 | 32 GB × attempt | 550 / 250 GB | `n2-highmem-8`, `n2-highmem-16`, `n2-highmem-16` |
| `SOMALIER_EXTRACT` | 1 | 4 GB × attempt | 100 / 20 GB | `e2-standard-2`, `e2-standard-4`, `e2-standard-8` |
| `BUILD_SAMPLESHEET`, `REGISTER_ANALYSIS` | 1 | 4 / 2 GB | default | `e2-standard-2` |
| `WRITE_RUN_RECORDS` | native `exec:`, in the Nextflow process | | | |

**Every process names one exact `machineType`**, chosen by `task.attempt`, so an OOM retry lands on a
bigger VM. Don't replace these with a list of types or a wildcard (`n2-*`). With
`NXF_CLOUDINFO_ENABLED=false`, Nextflow records a task's machine type in the trace only when it was given
an exact type. Seqera's cost estimates come from that record, so with a list or wildcard Google Batch
still picks a VM, but the type isn't recorded and the task gets no cost estimate ([RISKS.md](RISKS.md)).
Each type leaves memory headroom over the request, because a request equal to the VM's memory may not
fit. The align processes use N2 because bwa-mem3 needs AVX2.

## Files

```
main.nf                       params validation, run-records check, sample sheet, per-SG destinations, pre-flight, ALIGN, somalier, registration, run records
nextflow.config               params, profiles (standard: Google Batch, docker, test)
nextflow_schema.json          launch-form schema (main.nf does its own validation)
conf/base.config              containers, resources, error strategy, per-SG resource labels, publishDir
subworkflows/local/align/     ALIGN
modules/local/<process>/      build_samplesheet, sort_index, align_to_cram, align_shard, merge_to_cram, somalier_extract, register_analysis, write_run_records
bin/build_samplesheet.py      cohorts or SG IDs, from any datasets -> sample sheet (needs only `metamist`)
bin/align_stream.sh           per-kind input -> bwa-mem3 BAM stream
bin/dedup_to_cram.sh          dupblaster -> sort -> CRAM 3.0 + metrics
bin/threads.sh                CPU split between the tools in a pipe
bin/check_alignments.sh       quickcheck + at least one record
test/                         pytest: sample sheet unit tests, stub runs of main.nf
```

## Running

```bash
# production (the launcher adds -work-dir and google.project)
nextflow run main.nf \
    --dataset <analysis ds> --sequencing_type genome --arGuid <guid> --cohorts COH1,COH2 \
    --align_container <image> --cpg_container <image>
  
# or realistically... 
analysis-runner seqera ...

# offline wiring check
nextflow run main.nf -profile test -stub --dataset ds --sequencing_type genome \
    --samplesheet sheet.tsv --ref_fasta refs/ref.fasta

# tests (sample sheet unit tests + stub runs; the stub runs are skipped without nextflow)
uv run pytest
```

`dupblaster_extra_args` is a string that is passed into `ALIGN_TO_CRAM` and `MERGE_TO_CRAM` and appended
verbatim to the dupblaster command, so quote it as you would for bash:

```bash
--dupblaster_extra_args "--read-name-format 'regex:^(?<su>\w+:\d+):(?<tile>\d+):\d+:\d+'"
```

## Status and next steps

**Done:** the multi-dataset design described above (ADR 0008, ADR 0009), implemented from
[PLAN.md](PLAN.md): `build_samplesheet.py` with `--access-level`, top-level queries and a `dataset` column;
per-SG destinations in `meta`; `publishDir` in place of the `output {}` block; `WRITE_RUN_RECORDS`; the
run-records check; per-SG resource labels; per-SG registration projects. `nextflow lint` is clean, and 40
tests pass. The stub runs span two datasets and cover every input kind, the exome layout, registration
paths and projects (including `accessLevel test`), injected failures, the run-records and pre-flight
clashes, -resume, cohorts and SG ID modes, a sheet row without a dataset, per-SG resource labels and
parameter validation.

**Not done yet:**
1. **A real run on a small test dataset.** Only stubs have run so far, so the process scripts, images,
   resources, Metamist access from Batch and registration are all untested.
2. **Tiny real-data fixtures** (a small contig with its index, and a few thousand pairs as FQ/BAM/CRAM)
   for a local docker run.
3. **Genotyping:** consume `ALIGN.out.all_crams`, and bring back the "a `gvcf` analysis means the SG is
   done" rule in the sample sheet.
