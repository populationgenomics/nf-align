# nf-align: design brief

Aligns the short-read sequencing groups of a set of Metamist cohorts, or a list of sequencing groups. The output for each is a
deduplicated, coordinate-sorted CRAM. It's a Nextflow (≥ 25.10, strict syntax) pipeline that runs on Google
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

One run covers **one dataset and one sequencing type** ([ADR 0001](adr/0001-one-dataset-one-sequencing-type-per-run.md)).

| param | meaning |
|---|---|
| `dataset` | Metamist project name, without the `-test` suffix |
| `sequencing_type` | `genome` or `exome` |
| `accessLevel` | `full` / `standard` (equivalent), or `test`. Selects the Metamist project (`<dataset>` or `<dataset>-test`) and the default bucket (`gs://cpg-<dataset>-main` or `-test`) |
| `arGuid` | the analysis-runner GUID, used in resource labels, the per-run record directory and the registered analyses. Mandatory |
| `cohorts` **or** `sgids` **or** `samplesheet` | exactly one: Metamist cohort IDs, Metamist SG IDs (queried directly), or a prepared sheet that skips the Metamist step |
| `ref_fasta` | the alignment reference. Its `.fai` and bwa-mem3 index (`<ref_fasta>.{fai,bwt.2bit.64,amb,ann,pac}`) must sit beside it ([ADR 0006](adr/0006-dragen-masked-reference-alt-unaware.md)) |
| `align_container`, `cpg_container` | null until the images are built; the `gcb` profile refuses to run without them |

`outputDir` defaults to the dataset bucket. The launcher supplies `-work-dir` and `google.project`;
the config never sets them.

## Flow

```
--cohorts | --sgids ─ BUILD_SAMPLESHEET ─┐          (or --samplesheet)
                               ▼
                         sample sheet ─ parse ─┬─ aligned ───────────────────────────────▶ (all_crams, for genotyping)
                                               └─ to align ─ pre-flight ─ ALIGN ─┬─ publish CRAM, .crai, metrics
                                                                                 └─ REGISTER_ANALYSIS (per SG)
                                                                    failed-SG report ◀─┘
```

### 1. Sample sheet (`BUILD_SAMPLESHEET`, `bin/build_samplesheet.py`)

- One GraphQL query, scoped to `project(name: <dataset>[-test])`. It fetches cohort status and the
  cohorts' SGs (`--cohorts`), or the named SGs directly (`--sgids`), with assay meta and completed, active
  `cram` analyses.
- Cohorts that are missing from the project, inactive or empty fail the run. So does an SG from another
  project. Every offender is reported in one error.
- `--cohorts`: SGs are de-duplicated across cohorts. SGs that aren't short-read, or aren't
  `sequencing_type`, are skipped with a warning.
- `--sgids`: each SG was asked for by name, so one that's missing from the project, or isn't short-read
  `sequencing_type`, fails the run instead of being skipped.
- Read locations are used **exactly as Metamist stores them**.
- Each SG gets one row per alignment unit, with a `kind` (glossary): `fastq` and `fastq_ora` (one row
  per lane), `cram` or `bam` (one row), or `aligned` (it already has a CRAM analysis; one row pointing at
  that CRAM).
- A run with no in-scope SG fails.
- The sheet is published to `alignment_samplesheets/<arGuid>/samplesheet.tsv`. A `--samplesheet` input
  is not republished.

Columns: `sg_id, kind, n_units, unit, r1, r2, reads, reads_index, reads_ref`. For a CRAM, `reads_ref`
is the reference it was written against; for ORA, it's the ORA reference tarball.

### 2. Pre-flight

Before anything is aligned, every CRAM, `.crai` and metrics path that the run would publish is checked
under `outputDir`. If any exist, the run fails and lists them all
([ADR 0002](adr/0002-preflight-check-instead-of-overwrite-false.md)). `aligned` rows are not checked:
for those, an existing CRAM is the expected state.

A run that passes the check records the outputs it claims in `<workDir>/nf-align-claims/<sessionId>.txt`.
`-resume` keeps the session ID, so a resumed run lets those outputs exist and republishes them from the
cached tasks. Anything else that already exists still fails the run.

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

Workflow outputs, with `mode 'copy'`. Paths come from three functions in `main.nf`, which are the
single source of truth for both publishing and registration:

| target | path (`exome/` prefix for exomes) |
|---|---|
| CRAM + index | `cramPath(meta)` = `[exome/]cram/<sg>.cram`, plus `.crai` |
| duplicate metrics | `metricsPath(meta)` = `[exome/]qc/markduplicates_metrics/<sg>.markduplicates-metrics` |
| per-run records | `runPath(name)` = `alignment_samplesheets/<arGuid>/{samplesheet.tsv, crams.csv, failed.csv}` |

Only newly aligned SGs are published.

### 5. Registration (`REGISTER_ANALYSIS`)

When `register` is true, this runs once per new CRAM. It gets `uriOf(outputDir.resolve(cramPath(meta)))`,
and polls for up to `register_wait_minutes` until the CRAM and `.crai` exist, because publishing is
asynchronous ([ADR 0005](adr/0005-registration-as-a-polling-process.md)). It then runs:

```
python3 -m cpg_utils.metamist_registration --project <dataset>[-test] --output <cram> --type cram \
    --secondary index=<cram>.crai --sgs <sg> --meta sequencing_type=<type> ar-guid=<arGuid>
```

Metrics are not registered.

### 6. Failures

- **Per-SG tasks:** an exit code of 137 or 140 (OOM) is retried twice, with memory scaled by the attempt
  number. Any other failure is **ignored**: that SG is dropped and the rest carry on
  ([ADR 0007](adr/0007-ignore-failed-sgs-and-report.md)).
- **Registration** is retried once and then ignored.
- **Spot preemption** is retried by `google.batch.maxSpotAttempts = 3`.
- **`BUILD_SAMPLESHEET`** and the pre-flight check fail the whole run.
- **Reporting:** dropped SGs are written to `alignment_samplesheets/<arGuid>/failed.csv` with the step
  they failed at (`align` or `register`), and logged as a warning. **The run still exits 0.**

## Resources (`conf/base.config`)

| process | CPUs | memory | disk (genome / exome) | machines |
|---|---|---|---|---|
| `ALIGN_TO_CRAM`, `ALIGN_SHARD` | 16 | 48 GB × attempt | 400 / 100 GB | `n2-*,n2d-*,c2d-*,c3d-*` (AVX2) |
| `MERGE_TO_CRAM` | 8 | 24 GB × attempt | 400 / 100 GB | any |
| `SORT_INDEX` | 8 | 32 GB × attempt | 550 / 250 GB | any |
| `BUILD_SAMPLESHEET`, `REGISTER_ANALYSIS` | 1 | 4 / 2 GB | default | any |

## Files

```
main.nf                       params validation, sample sheet, pre-flight, ALIGN, registration, failed-SG report, outputs
nextflow.config               params, outputDir, workflow.output settings, profiles (gcb, docker, test)
nextflow_schema.json          launch-form schema (main.nf does its own validation)
conf/base.config              containers, resources, error strategy, resource labels
subworkflows/local/align/     ALIGN
modules/local/<process>/      build_samplesheet, sort_index, align_to_cram, align_shard, merge_to_cram, register_analysis
bin/build_samplesheet.py      cohorts or SG IDs -> sample sheet (needs only `metamist`)
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
    --dataset <ds> --sequencing_type genome --arGuid <guid> --cohorts COH1,COH2 \
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

**Done:** everything above. `nextflow lint` is clean, and 24 tests pass, covering every input kind,
exome layout, registration paths, injected failures, the pre-flight clash, -resume, cohorts and SG ID modes and
parameter validation.

**Not done yet:**
1. **Build the images.** `align_container` needs bazam, bwa-mem3, dupblaster, samtools and orad;
   `cpg_container` needs metamist and cpg_utils.
2. **Build the bwa-mem3 index** for `Homo_sapiens_assembly38_masked.fasta`. That needs about 76 GB of RAM
   ([RISKS.md](RISKS.md)). Write it next to the FASTA (`bwa-mem3 index <fasta>`), then set the container
   defaults in `nextflow.config`.
3. **A real run on a small test dataset.** Only stubs have run so far, so the process scripts, images,
   resources, Metamist access from Batch and registration are all untested.
4. **Tiny real-data fixtures** (a small contig with its index, and a few thousand pairs as FQ/BAM/CRAM)
   for a local docker run.
5. **Genotyping:** consume `ALIGN.out.all_crams`, and bring back the "a `gvcf` analysis means the SG is
   done" rule in the sample sheet.
